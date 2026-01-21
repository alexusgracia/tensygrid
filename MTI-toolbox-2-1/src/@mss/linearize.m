function [lsys] = linearize(msys,x,u)
% <linearize> Linearizes an explicit MTI model in CPN format
%   
% Input arguments: 
%  - msys: mss object
%  - x  : scalar or vector state operating point  (or op struct / vector [x;u])
%  - u  : scalar or vector input operating point (optional if op used)
%
% Output arguments: 
%  - linSs: ss object 
%
% Remarks:
%  - If an output equation of the mss object is not specified
%    all states are outputs and feedthrough is set to zero.
%
% Example
% linsys=linearize(sys,x_op,u_op)
%
% For detailed documentation see <a href="matlab:open((which('linearizeDoc.html')))">here</a>

% Enrico Uhlenberg, Marah Engels, Torben Warnecke, Leandro Samaniego,
% Carlos Cateriano Yáñez, Christoph Kaufmann, Leona Schnelle, Georg
% Pangalos, and Gerwald Lichtenberg - 12/06/2024

% Accept unified op argument (struct / vector / cell) or legacy (x,u)
if nargin >= 2 && (isstruct(x) || iscell(x) || (isnumeric(x) && isvector(x) && numel(x)~=msys.n))
    % Treat x as 'op' and parse it
    [x,u] = i_parseOpMSS(msys, x);
else
    % Legacy path: x (required), u optional
    if nargin < 2
        error('No operating point for state variables x or input variables u specified. Cannot execute linearization.');
    end
    if nargin < 3
        if msys.m~=0
            error('No operating point for input variables u specified even though the system has inputs. Cannot execute linearization.');
        else
            u = [];
        end
    end
end

% Dimension consistency checks 
if msys.n ~= length(x)
    error("Number of states of the operating point argument x do not match with number of states of the mss object.")
    
end

% ===== Local helper  =====
function [x,u] = i_parseOpMSS(sys, op)
% Returns column vectors x (n x 1) and u (m x 1 or [])

n = sys.n; m = sys.m;
x = []; u = [];

if isstruct(op)
    if isfield(op,'x'), x = op.x; end
    if isfield(op,'u'), u = op.u; end
elseif isnumeric(op)
    v = op(:);
    if numel(v) == n+m
        x = v(1:n);
        u = v(n+1:end);
    elseif numel(v) == n && m==0
        x = v;
        u = [];
    else
        error('linearize:BadOpMSS', 'Numeric op length must be n+m (or n if m=0).');
    end
elseif iscell(op)
    if numel(op) < 1
        error('linearize:BadOpMSS', 'Cell op must contain at least {x}.');
    end
    x = op{1};
    if numel(op) >= 2, u = op{2}; end
else
    error('linearize:BadOpMSS', 'Unsupported op type for mss. Use struct with fields x,u or vector [x;u].');
end

% Shape & presence checks
if isempty(x)
    error('linearize:MissingX', 'State operating point op.x (x) is required for mss.');
end
x = x(:);
if numel(x) ~= n
    error('linearize:DimX', 'Length of x (%d) must equal number of states n=%d.', numel(x), n);
end

if isempty(u)
    if m>0
        error('linearize:MissingU', 'Input operating point op.u (u) is required for systems with m>0.');
    else
        u = [];
    end
else
    u = u(:);
    if numel(u) ~= m
        error('linearize:DimU', 'Length of u (%d) must equal number of inputs m=%d.', numel(u), m);
    end
end
end

if nargin==3
   if msys.m ~= length(u)
        error("Number of inputs of the operating point argument u do not match with the number of inputs of the mss object.")
   end    
end

% Get linear state-space model of type ss

% Construct vector containing operating point of x and u
xu=[x(:);u(:)];

if issparse(msys.F.U)
    
    AB=msys.F.cpn2LinSparse(xu);
    CD=msys.G.cpn2LinSparse(xu);

    % generate the output as sparse linear state space model object
    lsys = sparss(AB(:,1:msys.n),AB(:,(msys.n+1):end),CD(:,1:msys.n),CD(:,(msys.n+1):end));

else % Linearization for non-sparse model
        
        % Jacobian of states    
        AB = msys.F.cpn2Lin(xu);
        
        % Jacobian of outputs
        if isempty(msys.G)
            warning('No Outputs defined, assuming states as outputs.')
            CD = [eye(msys.n),zeros(msys.n,msys.m)];
        else 
            CD = msys.G.cpn2Lin(xu);
        end,  
        
        % Extract linear state space model
        lsys = ss(AB(:,1:msys.n),AB(:,(msys.n+1):end),CD(:,1:msys.n),CD(:,(msys.n+1):end));
end



% Handle discrete-time models
    lsys.Ts = msys.ts;



end

