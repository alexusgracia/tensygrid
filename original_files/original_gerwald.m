%% Start file for TenSyGrid Hackaton 16.12.2025
% linearization for iMTI models in CPN representation 
% very basic, functional, MATLAB
%
% gerwald.de 16.2.2025

clear
check = 0;   % 0 = sparse = scalable = fast | 1 = full = only for checks

%% Define problem

% Dimensions
n = 2;          % number of states
m = 1;          % number of inputs
p = 0;          % number of outputs
r = 7;          % rank

q = n;          % number of equations
N = 2*n+m+p;    % total number of signals

% iMTI model in CPN representation  
S = rand(N,r);          % Structure matrix
S = round(2*S-1);       % Project to {-1,0,1} - un/comment for tests
S = sparse(S);          % Convert to sparse - un/comment for tests
P = rand(q,r);          % Parameter matrix
P = P.*(P<0.3);         % Set some zeros - un/comment for tests
P = sparse(P);          % Convert to sparse - un/comment for tests

% linearization point
dx = rand(n,1);         % state derivative
disp(dx);
x = rand(n,1);          % state
u = rand(m,1);          % input
y = rand(p,1);          % output
v = [dx;x;u;y];         % signal vector

tic     % start the cputime clock
%% Linearization  (Scalable = sparse & low rank)
X = S.*v-abs(S);        % Important term, better than S.*(v*ones(1,r)) 
icrit = (X==-1);        % get numerical critical indices
if any(icrit(:))
    error('specials to be implemented ?')
end
X = spfun(@(x) x+1 ,X);     % add one only to the nonzero elements
[rowi,coli,val] = find(X);  % get the indices and values of X
Y = accumarray(coli,val,[r 1],@prod)'; % product over nonzero elements
X = spfun(@(x) 1./x, X);    % invert only nonzero elements
F = S.*Y.*X;                % compute factor matrix -> Paper
toc                         % get first cputime 

%% Full version (NOT scalable ! Only for didactics & checks)
if check 
    FC = S;                 % initialize product matrix
    X = 1-abs(S)+S.*v;      % evaluate factors 1-|S|+Sv (norm-1) 
    %X = sqrt(1-S.^2)+S.*v;  % alternative (norm-2)
    for k=1:N        % Compute products for each variable & summand
        Y = X;              % init with all factors 
        Y(k,:) = S(k,:);    % only the linear coefficient for this variable 
        FC(k,:) = prod(Y);  % compute the product over columns (lot of 1s) 
    end
end

%% Extract linear system matrices
EABC = P*F';                    % Compute combined matrix
E = -EABC(:,1:n);               % dx matrix        
A = EABC(:,(n+1):2*n);          % system matrix
B = EABC(:,(2*n+1):(2*n+m));    % input matrix

%% Do analysis MATLAB based

% Local stability: dominant (geralized) eigenvalues
eigs(A,E)
% Controllability: Krylov subspace is full rank
try
    rc = rank(ctrb(A,B));
    if (rc<n)
        disp('not controllable'),
    else 
        disp('controllable');
    end
catch
    warning('not computable')
end
toc         % get first cputime 

%% Debug
if check
    disp(full([F;FC]));
end