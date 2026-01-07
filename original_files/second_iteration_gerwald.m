%% Linearization and analysis for iMTI models in CPN representation 
% basic & functional MATLAB code
% gerwald.de 16.-19.12.2025

clear
check = 0;     % 0 = sparse = scalable = fast | 1 = full = only for checks
singular = 0;  % 0 = random example | 1 = singular worst case  
krylov = 0;    % 0 = no computation | 1 = compute Krylov subspace rank
boole = 0;     % 0 = ones are double | 1 = ones are logical

%% Dimensions
n = 1000;       % number of states
m = 10;         % number of inputs
p = 0;          % number of outputs  (currently: only 0 possible)
r = 15000;      % tensor rank

if singular         % Overwrite dimensions
    n = 2;          % number of states
    m = 1;          % number of inputs
    p = 0;          % number of outputs  (currently: only 0 possible)
    r = 6;          % rank
end

if check & ~singular & n>10
    warning('consider checking a small example or turn check = 0')
end

q = n;          % number of equations
N = 2*n+m+p;    % total number of signals

%% iMTI Model 
if singular
    % Simple example 
    % dx1 = -x1 + u
    % dx2 = -x2 +x1*x2
    S = [1 0 0 0 0 0;0 1 0 0 0 0;0 0 1 0 0 1;0 0 0 1 0 1;0 0 0 0 1 0];
    S = sparse(S);                              % Structure matrix
    P = sparse([-1 0 -1 0 1 0;0 -1 0 -1 0 1]);  % Parameter matrix
    % singular linearization point
    dx = ones(n,1);         % state derivative      
    x = zeros(n,1);         % state
    u = ones(m,1);          % input
    y = ones(p,1);          % output
else
    % random iMTI model in CPN representation  
    P = sprand(q,r,1/n);    % Random Parameter matrix
    S = sprand(N,r,1/n);    % Random Structure matrix 
    if boole
        S = S>0;            % Project to Boolean {0,1}
    end
    % random linearization point
    dx = rand(n,1);         % state derivative
    x = rand(n,1);          % state
    u = rand(m,1);          % input
    y = rand(p,1);          % output
end

v = [dx;x;u;y];             % signal vector of linearzation point

%% Linearization  (Scalable = sparse & low rank)
tic                             % start the cputime clock
X = S.*v-abs(S);                % sparse matrix of all Factors - 1
izero = sparse(S.*v==abs(S)&S); % indices of nonzeros of S where X==0
imone = sparse(X==-1);          % indices where X==-1 
ctwo = sum(imone)>1;            % get colums where more than one X==-1
X = spfun(@(x) x+1, X);         % add one to all nonzero elements
%X(izero)=1 corrects indices in S but not X and X(imone)=1 for one X=-1
X(izero|imone) = 1;             % do both corrections in one
[rowi,coli,val] = find(X);      % get the indices and values of X
Y = accumarray(coli,val,[r 1],@prod)'; % all products at operation point
Y(ctwo) = 0;                    % Columns with 2 zero factors will be zero
X = spfun(@(x) 1./x, X);        % invert only nonzero elements
F = S.*Y.*X;                    % compute factor matrix -> Paper
EABC = P*F';                    % Compute combined LTI model matrix
disp('Linearization: done')     % display when linearization is done
t1 = toc;                       % get cputime for linearization

%% Full Linearization  (NOT scalable ! Only for didactics & checks)
if check 
    FC = S;                 % initialize product matrix
    X = 1-abs(S)+S.*v;      % evaluate factors 1-|S|+Sv (norm-1) 
    %X = sqrt(1-S.^2)+S.*v; % alternative (norm-2)
    for k=1:N               % Compute products for each variable & summand
        Y = X;              % init with all factors 
        Y(k,:) = S(k,:);    % only the linear coefficient for this variable 
        FC(k,:) = prod(Y);  % compute the product over columns (lot of 1s) 
    end
    disp('Linearization (full): done')   % display when linearization is done
    t1 = toc; 
end
%% Extract LTI matrices
E = -EABC(:,1:n);               % Extract "mass" matrix E       
A = EABC(:,(n+1):2*n);          % Extract system matrix A
B = EABC(:,(2*n+1):(2*n+m));    % Extract input matrix B
   
%% Stability Analysis 
try   
    lambda = eigs(A,E,4,'largestreal');  % dominant generalized eigenvalues
    if max(real(lambda)>0)
        disp('Eigenvalues: unstable'),
    else 
        disp('Eigenvalues: stable');
    end
catch
    disp('Eigenvalues: not computable')
end
t2 = toc;                       % cputime after stability analysis
disp(['Lin: ' num2str(t1), '  Eig: ' num2str(t2-t1) ])

%% Krylov subspace investigation
if krylov  
    try
        % only for performance here, has to be adapted to implicit!!      
        rc = rank(ctrb(A,B));   % rank check 
        if (rc<n)
            disp('Controllability matrix: rank < n'),
        else 
            disp('Controllability matrix: rank n');
        end
    catch
        disp('Controllability matrix: not computable')
    end
    t3 = toc;                   % cputime after controllability analysis
    disp(['Krylov: ' num2str(t3-t2)])
end

%% Debug
if check
    if n<10                     % only makes sense for small problems
        disp(full(P*FC'))       % display concatination of LTI matrices
        disp(full(EABC))        % display concatination of LTI matrices
    end
    try
        disp(lambda)
    end
    if any(abs(F-FC)>10*eps)
        warning('results for factors differ')
    end
end 