function [sys] = mlinearize(model,low_bnd,up_bnd,level,max_order)
% Obtains a multilinear model from a Simulink model, using sparse 
% grid interpolation methods.
%
%  MLINEARIZE takes a Simulink model, a domain with lower and upper bounds, 
%  the sparse grid level, and the maximum multilinear order as inputs. It 
%  outputs a multilinear time-invariant state-space model as an mss object.
%
%  input parameters:
%     -model
%     -low_bnd & up_bnd: bounded region intervals
%               low_bnd = [x1_l,...,xn_l,u1_l,...,um_l,y1_l,...,yp_l],
%                up_bnd = [x1_u,...,xn_u,u1_u,...,um_u,y1_u,...,yp_u].
%
%     -level: number of points on the sparse grid, O((n+m)^level).
%
%     -max_order: maximum multilinear order allowed.
%
%  output parameters:
%      - sys: continous or discrete mss object.
%
%  For detailed documentation see <a href="matlab:open(which('mlinearizeDoc.html'))">here</a>
%  See also mss
 
%  Leandro Samaniego, Kai Kruppa - 11.06.2024


% Check if Sparse toolkit is live
res = environmentChecker(true,"Sparse_Toolkit");
if res == 0 
 error(['Sparse Grids MATLAB kit not found in the MATLAB path. ', ...
        'Ensure that the Sparse Grids MATLAB kit (release 23-5 onwards) is installed ', ...
        'and added to the MATLAB path./']);
end


% Compile Simulink Model
[sys,~,str,ts] = feval(model,[],[],[],'compile');

settings = getActiveConfigSet(model);
disp("Order of states taken by  Simulink:")
disp(str)

% Get sample time
tsim = ts(end,1);
if tsim == 0
    tsim = [];
end

% Number of continuous or discrete states x, inputs u and outputs y
if sys(1)==0 && sys(2)~=0
%    only discrete states
    nx = sys(2);
    disc_flag = 1;
    noStates_flag = false;
elseif sys(1)~=0 && sys(2)==0
%     only continuous states
    nx = sys(1);
    disc_flag = 0;
    noStates_flag = false;
elseif sys(1)==0 && sys(2)==0
%     no states
    nx = 0;
    noStates_flag = true;
    disc_flag = 0;
    warning("System has no discrete or continuous states.");
else
    error(['Models with discrete AND continuous states are not supported.' ...
        'Avoid using Unit Delay blocks to break algebraic loops.'])
end
nu = sys(4);
ny = sys(3);
if nx+nu+ny ~= length(up_bnd)
    error('Wrong size for the argument up_bnd')
end
if nx+nu+ny ~= length(low_bnd)
    error('Wrong size for the argument low_bnd')
end

% Sparse Grid is created on [-1,1]^(d) to the specifed level of interpolation.
d = nx + nu;
Sr = createScaledGird(d,level);
xu_scal = Sr.knots;

% Create scaling factors (a,b) for xu_scal = a*xu+b
[a,b] = scaleFactors(up_bnd,low_bnd);

% Get sampling points xu from xu_scal
xu = scaleGrid(xu_scal,a,b,nx,nu);
% Number of sampled points
m = size(xu,2);      

% Evaluation of dx/dt = f(x,u) and y = g(x,u) on the sparse grid

% Initialize state derivative and output vectors.
dx = zeros(nx,m);
y = zeros(ny,m);
 
% Change 'Initial state is array' Diagnostic setting to 'none'
% before sampling. 
set_param(settings,'InitInArrayFormatMsg','none');

% Sampling Model on sparse grid xu
for i=1:m
    y(:,i) = feval(model,tsim,xu(1:nx,i),xu(nx+1:d,i),'outputs');
    if noStates_flag == false
        if disc_flag == 1
            dx(:,i) = feval(model,tsim,xu(1:nx,i),xu(nx+1:d,i),'update'); 
        elseif disc_flag == 0
            dx(:,i) = feval(model,tsim,xu(1:nx,i),xu(nx+1:d,i),'derivs');
        else
            error('Models with discrete AND continuous states are not supported yet')
        end
    end
end

% Finish the simulation of the model. 
set_param(settings,'InitInArrayFormatMsg','warning');
feval(model,[],[],[],'term');

% Scale state derivtive dx/dt and output y vectors into the domain [-1,1]^d
[dx_scal,y_scal] = scaleModel(dx,y,a,b,nx,nu,ny);

% Reduced index by maximum multilinear order
red_index = buildIndex(d,max_order);

% Calculation of the scaled F and G matrices
[F_scal,G_scal] = scaledFG(xu_scal,dx_scal,y_scal,red_index,Sr,nx,nu,ny,m,disc_flag);


% Scale the coefficents back into the original domain
[F_trafo,G_trafo] = deScaledFG(F_scal,G_scal,a,b,nx,nu,ny,red_index);

% Create the entire F and G Matrices
F = zeros(nx,2^(nx+nu));
G = zeros(ny,2^(nx+nu));
F(:,red_index) = F_trafo;
G(:,red_index) = G_trafo;

% Create mss object
if (disc_flag==1)
    sys = mss(F,G,tsim);
    sys.stateName=str;
else
    sys = mss(F,G);
    sys.stateName=str;
end

end

%% Sub-functions
function [Sr] = createScaledGird(d,level)
    knots=@(n) knots_CC(n,-1,1,'nonprob'); % Clenshaw-Curtis knots
    w = level;
    S = create_sparse_grid(d,w,knots,@lev2knots_doubling);
    Sr = reduce_sparse_grid(S);
end

% Scaling factors for states inputs and outputs
function [a,b] = scaleFactors(up_bnd,low_bnd)  
    a = 2./(up_bnd - low_bnd);
    b = -(up_bnd + low_bnd)./(up_bnd - low_bnd);

    % Check for same bounds in the domain
    logic1 = isinf(a)|isnan(a)|isinf(b)|isnan(b);
    if any(logic1)
        warning("Same bound(s) in position(s): " + mat2str(find(logic1)));
        a(isinf(a)|isnan(a)) = 0;
        b(isinf(b)|isnan(b)) = 0;
    end

    % Check for upper bound smaller than lower bounds
    logic2 = up_bnd<low_bnd;
    if any(logic2)
        warning("Lower bound(s) greater than upper bound at position(s): " + mat2str(find(logic2)));
    end
end

% Get xu from xu_scal = a*xu+b
function xu = scaleGrid(xu_scal,a,b,nx,nu)
    xu = (xu_scal-b(1:nx+nu)').*(a(1:nx+nu)'.^-1);
    xu(isinf(xu)|isnan(xu)) = 0;
end

% Get scaled version of dx and y
function [dx_scal,y_scal] = scaleModel(dx,y,a,b,nx,nu,ny)
    dx_scal = (a(1:nx)').*(dx);
    y_scal = (a(nx+nu+1:nx+nu+ny)').*y + b(nx+nu+1:nx+nu+ny)';
end

function [red_index] = buildIndex(nxu, max_order) 
% Obtain the indices of the monomial which order is less or equal to max_order
    if max_order > nxu
        max_order = nxu;
        warning('Maximal multilinear order is higher than the number of variables. => Maximal multilinear order is set to the number of variables.')
    end
    ind_mx = 0;
    for ind1 = 1:nxu
        ind_old = ind_mx;
        ind_mx = zeros(ind1+1,2^ind1);
        ind_mx(1) = 1;
        for ind2 = 1:ind1
            ind_mx(1,ind2+1) = 2^(ind2-1)+1;  
        end
        for ind2 = 2:ind1
            ind_mx_temp = [ind_old(ind2,1:length(find(ind_old(ind2,:)))) 2^(ind1-1)+ind_old(ind2-1,1:length(find(ind_old(ind2-1,:))))];
            ind_mx(ind2,1:length(ind_mx_temp)) = ind_mx_temp;
        end
        ind_mx(ind1+1,:) = 1:2^ind1;
    end
    red_index_aux = ind_mx(max_order,:);
    red_index = red_index_aux(red_index_aux~=0);
end

function [F_scal,G_scal] = scaledFG(xu_scal,dx_scal,y_scal,red_index,Sr,nx,nu,ny,m,disc_flag)

    % Size of monomial vector
    n_xu = length(red_index);
    
    % Initialize monomial vector
    mXU = zeros(n_xu,m);
    
    % Calculate scaled basis u(x_scal,u_scal)
    for i=1:m
    mXU(:,i) = reducedMXU(xu_scal(:,i),red_index)';
    end

    if disc_flag == 1 
         % Orthonormalization of the basis u(x_scal,u_scal) into
         % u_hat(x_scal,u_scal) in the discrete case
         [Q,R] = gson(mXU'); % Q-R matrix decomposition
         mXU = Q';
    
         % Obtain the orthonormal-scaled phi and gamma coefficients
         phi_hat_scal = dx_scal*mXU';
         gamma_hat_scal = y_scal*mXU';
    
         % Transform coefficents into regular basis u(x_scal,u_scal)
         F_scal = phi_hat_scal/(R');
         G_scal = gamma_hat_scal/(R');
         
    elseif disc_flag == 0
         % Orthonormalization of the basis u(x_scal,u_scal) into
         % u_hat(x_scal,u_scal) in the continouos case
         
         % Orthonormalization coefficents
         on_coeff = orthoNormalCoeff(nx + nu,red_index); 
         for i=1:m
            mXU(:,i) = mXU(:,i)'.*on_coeff; % u_hat(x_scal,u_scal)
         end
    
         % Multiply dx_scal with each element of u_hat(x_scal,u_scal).
         dx_u_hat = zeros(size(mXU,1),size(mXU,2),nx);
         for k=1:nx
            for i=1:n_xu
                dx_u_hat(i,:,k) = dx_scal(k,:).*mXU(i,:);  % f(x_scal,u_scal)*u_hat(x_scal,u_scal)
            end
         end
     
        % Multiply y_scal with each element of u_hat(x_scal,u_scal).
         y_u_hat = zeros(size(mXU,1),size(mXU,2),ny);
         for k=1:ny
            for i=1:n_xu
                y_u_hat(i,:,k) = y_scal(k,:).*mXU(i,:); % g(x_scal,u_scal)*u_hat(x_scal,u_scal)
            end
         end
    
         % Obtain the normalized scaled phi and gamma coefficients
         phi_hat_scal = zeros(nx,n_xu);
         gamma_hat_scal = zeros(ny,n_xu);
         for i = 1:nx
            phi_hat_scal(i,:) = quadrature_on_sparse_grid(dx_u_hat(:,:,i),Sr)'; % integrates -> f(x_scal,u_scal)*u_hat(x_scal,u_scal) over Sr
         end
         for i = 1:ny
            gamma_hat_scal(i,:) = quadrature_on_sparse_grid(y_u_hat(:,:,i),Sr)'; % integrates -> g(x_scal,u_scal)*u_hat(x_scal,u_scal) over Sr
         end
         
         % Transform coefficents into regular basis u(x_scal,u_scal)
         F_scal = phi_hat_scal.*on_coeff;
         G_scal = gamma_hat_scal.*on_coeff;
    end
    end
    
    % Obtain partial F and G matrices from scaled version
    function [F_trafo,G_trafo]=deScaledFG(F_scal,G_scal,a,b,nx,nu,ny,red_index)
        AB_n = 1;
        for ind1 = 1:nx+nu
            AB_n = [1*AB_n 0*AB_n;b(ind1)*AB_n a(ind1)*AB_n];
        end
        AB_n = AB_n(red_index,red_index);
        a_diag_F = diag(a(1:nx).^-1);
        F_trafo = a_diag_F*(F_scal*AB_n);
        a_diag_G = diag(a(nx+nu+1:nx+nu+ny).^-1);
        offset_G = zeros(ny,length(red_index));
        offset_G(:,1) = b(nx+nu+1:nx+nu+ny)';
        G_trafo = a_diag_G*(G_scal*AB_n - offset_G);
end

function [mXU] = reducedMXU(xu, red_index)
    % Creates the monomial vector, while excluding the monomials with
    % degree hihger than max_order.
    
    % xu = [x;u];
    nxu = length(xu);
    beta = 2^nxu;
    switch nxu
        case 1
            mXU = [1 xu(1)]';
        case 2
            mXU = [ 1, xu(1), xu(2), xu(1)*xu(2)]';
        otherwise
            mXU = zeros(1,beta); 
            mXU(:,1:4) = [ 1, xu(1), xu(2), xu(1)*xu(2)]; 
            for i = 2:nxu-1
                mXU(:,1:2^(i+1)) = [mXU(:,1:2^i) mXU(:,1:2^i).*xu(i+1)]; 
            end
            mXU = mXU';
    end
    mXU = mXU(red_index);
end

function [norm_coeff] = orthoNormalCoeff(n,red_index)
    % Computes a set of coefficients q, to transfrom the monomial m(x) to
    % the normalized monomial m_hat(x), such that m_hat(x) = q.*m(x).
    beta = 2^n;
    norm_coeff = zeros(1,beta);
    q = zeros(1,n+1);

    % Calculate the coefficients.
    %  The coefficients arranged with respect to the degree of the components in m(x).
    %  For example,  m(x) = [1  x1 x2 x1*x2 x3 x1*x3 x2*x3 x1*x2*x3]',
    %                   q = [q0 q1     q2                     q3], 
    %  where the subindex in q_j is the degree of the corresponding monomial. 
    %  Then
    %                m_hat(x) = [q0 q1*x1 q1*x2 q2*x1*x2 q1*x3 q2*x1*x3 q2*x2*x3 q3*x1*x2*x3]'.
    %                m_hat(x) = q.*m(x);
    %
    % Coefficients values
    for i=0:n
       q(:,i+1) = sqrt(3^i*2^(-n));
    end

    %Arrange the coeffcients with respect to the monomial vector.
    % The implementation uses the fact that,
    %   for n = 3, m(x)_3 = [[1 x1 x2 x1x2], x3.*[x1 x2 x1x2]]',
    %   for n = 4, m(x)_4 = [m(x)_3', x4.*m(x)_3']'.
    %   and so on for othervalues of n.

    % A varible indx is constructed as a row vector with the degree of the monomials in m(x).
    indx = zeros(1,beta); 
    indx(:,1:4) = [0 1 1 2]; % initial monomial degrees
    for i = 2:n-1
        indx(:,1:2^(i+1)) = [indx(:,1:2^i) indx(:,1:2^i)+1]; 
    end

    % Adding one to indx to match the positions of the coefficients in the normalized monomial vector.
    indx = indx + 1; 
     
    % Order the coefficients in the correct position.
    for i=1:beta
        norm_coeff(:,i) = q(indx(:,i));
    end

    if nargin>1
        % Reduced index, tells the function which coefficents to drop.
        norm_coeff = norm_coeff(red_index);
    end
end

function [Q, R] = gson(X)
    % Gram-Schmidt orthonormalization which produces the same result as [Q,R]=qr(X,0)
    % Written by Mo Chen (sth4nth@gmail.com).
    [d,n] = size(X);
    m = min(d,n);
    R = zeros(m,n);
    Q = zeros(d,m);
    for i = 1:m
        R(1:i-1,i) = Q(:,1:i-1)'*X(:,i);
        v = X(:,i)-Q(:,1:i-1)*R(1:i-1,i);
        R(i,i) = norm(v);
        Q(:,i) = v/R(i,i);
    end
    R(:,m+1:n) = Q'*X(:,m+1:n);
end
