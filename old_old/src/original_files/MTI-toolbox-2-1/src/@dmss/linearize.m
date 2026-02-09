function [sparse_dss_sys,dss_sys] = linearize(msys,xpop,xop,uop,yop,zop)
%<linearize> Linearize an dmss-object around an operating point (op) to get a
% linear descriptor state-space model (dss)
%
% Input arguments:  
%  - msys : (hybrid) implicit multilinear state-space model as dmss object
%  - op.x  : scalar or vector state operating point (or op struct / vector [x;u])
%  - op.u  : scalar or vector input operating point
%  - op.xp : scalar or vector state derivative operating point 
%  - op.y  : scalar or vector algebraic operating point 
%  - op.z  : scalar or vector boolean operating point 
%
% Output arguments: 
%  - sparse_dss: sparse descriptor state-space model
%
% Remarks:
%  - If an output equation of the mss object is not specified
%    all states are outputs and feedthrough is set to zero.
%
% Example
% ldss=linearize(dsys,op)
%
% For detailed documentation see <a href="matlab:open((which('linearizeDoc.html')))">here</a>
%
% 15.12.2024, Christoph Kaufmann


% check for the number of arguments and input types

% Unified op argument support (struct / vector / cell): allow calling as
%   linearize(msys, op)
% while retaining legacy signature (xpop,xop,uop,yop,zop)
if nargin == 2 && (isstruct(xpop) || iscell(xpop) || (isnumeric(xpop) && isvector(xpop)))
    [xpop,xop,uop,yop,zop] = parseOpDMSS(msys, xpop);
end


if nargin < 2 || isempty(xop)
    error ('No data for the state equilibrium provided. Please provide a scalar, or vector with the state values of the operating point.')
end

if msys.n ~= length(xop)
    error("Number of states of the operating point argument x do not match with number of states of the dmss object.")
end


% Dimension consistency checks 
if msys.n ~= length(xop)
    error("Number of states of the operating point argument x do not match with number of states of the dmss object.")
elseif ~isempty(uop) && msys.m ~= length(uop)
    error("Number of inputs of the operating point argument u do not match with the number of inputs of the dmss object.")
end 

% 


%--------------------------------------------------------------------------
% Linearize model
%--------------------------------------------------------------------------

J=jacobian(msys,xpop,xop,uop,yop,zop); 

%--------------------------------------------------------------------------
% Extract descriptor state-space model
%--------------------------------------------------------------------------

E=[J.equality.stateDerivative,zeros(msys.nEq,msys.nEq-msys.n)]; % E needs to be extended

A=-[J.equality.state, J.equality.algebraic]; % both must be added  
B=-J.equality.input;

C=eye(msys.n+msys.p);

D=zeros(msys.n+msys.p,size(B,2));

%--------------------------------------------------------------------------
% Create descriptor state-space model
%--------------------------------------------------------------------------


% storing as sparse descriptor state-space model

sparse_dss_sys=sparss(sparse(A),sparse(B),sparse(C),sparse(D),sparse(E),msys.ts);
sparse_dss_sys.InputName=msys.inputName;
sparse_dss_sys.OutputName=[msys.stateName,msys.algebraicName];


% storing as full descriptor state-space model

    if nargout>1
    dss_sys = dss(A,B,C,D,E,msys.ts);
    dss_sys.stateName=[msys.stateName,msys.algebraicName];
    dss_sys.inputName=msys.inputName;
    dss_sys.outputName=dss_sys.stateName;

    end

end




% ===== Local parser   =====
function [xp,x,u,y,z] = parseOpDMSS(sys, op)
% Returns column vectors xp (n x 1), x (n x 1), u (m x 1), y (p x 1), z (q x 1)

n = sys.n; m = sys.m; p = sys.p; q = sys.q;
xp = []; x = []; u = []; y = []; z = [];

if isstruct(op)
    if isfield(op,'xp'), xp = op.xp; end
    if isfield(op,'x'),  x  = op.x;  end
    if isfield(op,'u'),  u  = op.u;  end
    if isfield(op,'y'),  y  = op.y;  end
    if isfield(op,'z'),  z  = op.z;  end
elseif isnumeric(op)
    v = op(:);
    expected = n+n+m+p+q;
    if numel(v) ~= expected
        error('linearize:BadOpDMSS', 'Numeric op length must be n+n+m+p+q = %d.', expected);
    end
    pos = 0;
    xp = v(pos+1:pos+n); pos = pos + n;
    x  = v(pos+1:pos+n); pos = pos + n;
    u  = v(pos+1:pos+m); pos = pos + m;
    y  = v(pos+1:pos+p); pos = pos + p;
    z  = v(pos+1:pos+q);
elseif iscell(op)
    if numel(op) >= 1, xp = op{1}; end
    if numel(op) >= 2, x  = op{2}; end
    if numel(op) >= 3, u  = op{3}; end
    if numel(op) >= 4, y  = op{4}; end
    if numel(op) >= 5, z  = op{5}; end
else
    error('linearize:BadOpDMSS', 'Unsupported op type for dmss. Use struct with fields xp,x,u,y,z.');
end

% col shape; [] allowed
toCol = @(v) reshape(v,[],1);
if ~isempty(xp), xp = toCol(xp); end
if ~isempty(x),  x  = toCol(x);  end
if ~isempty(u),  u  = toCol(u);  end
if ~isempty(y),  y  = toCol(y);  end
if ~isempty(z),  z  = toCol(z);  end

% Basic dimension checks for provided fields
if ~isempty(xp) && numel(xp) ~= n
    error('linearize:DimXP', 'Length of xp (%d) must equal n=%d.', numel(xp), n);
end
if ~isempty(x) && numel(x) ~= n
    error('linearize:DimX', 'Length of x (%d) must equal n=%d.', numel(x), n);
end
if ~isempty(u) && numel(u) ~= m
    error('linearize:DimU', 'Length of u (%d) must equal m=%d.', numel(u), m);
end
if ~isempty(y) && numel(y) ~= p
    error('linearize:DimY', 'Length of y (%d) must equal p=%d.', numel(y), p);
end
if ~isempty(z) && numel(z) ~= q
    error('linearize:DimZ', 'Length of z (%d) must equal q=%d.', numel(z), q);
end
end

