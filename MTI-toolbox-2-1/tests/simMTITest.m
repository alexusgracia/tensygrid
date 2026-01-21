% Test header
function tests = simMTITest
    tests = functiontests(localfunctions);
end

function BackWardsCompatibilityTest(testCase)
    BackWardsCompatibility(testCase,0);
    
end

function BackWardsCompatibilityBooleanTest(testCase)
    BackWardsCompatibility(testCase,1);
end

function ContinuousSimulationTest(testCase)
    A = [-1.5 -3; 3 -1];
    B = [1.3; 0];
    C = [1.15 2.3];
    D = [0];
    lsys = ss(A, B, C, D);

    x0 = [-0.2 0.3];
    t = 0:0.05:8;
    u = zeros(length(1),1);
    u(t>=2) = 1;

    [y_lin, t_lin, x_lin] = lsim(lsys, u, t, x0, 'foh');

    F = CPN1(diag(ones([3,1])), [A, B]);
    G = CPN1(diag(ones([3,1])), [C, D]);
    
    msys = mss(F,G);
    
    [y_mlin,t_mlin,x_mlin] = msim(msys, u', t, x0);
    
    y_lin_interp = interp1(t_lin, y_lin, t_mlin);
    [R,p] = corrcoef(y_mlin, y_lin_interp);

    verifyGreaterThan(testCase, R, 0.99*ones(2));

end




function BackWardsCompatibility(testCase,boolean)
    sys = drmss(2,2,2,5,0.2,boolean);
    ts = 1;
    msys = mss(CPN1(sys.F.U, sys.F.phi),CPN1(sys.G.U, sys.G.phi), ts);
    x0 = [ 0;0];
    u = [0 1 2 3 4; 0 1 2 3 4]';
    t = 5;

    [y_sim_legacy, x_sim_legacy] = LegacyMsim(sys,u',t,x0);
    [y_sim, x_sim] = LegacyMsim(msys,u',t,x0);
    
    osim = simMTI(msys);
    [y_result,t_result,x_result]= osim.simulate(u,t,x0);

    verifyEqual(testCase,{y_sim_legacy, x_sim_legacy},{y_sim, x_sim});
    verifyEqual(testCase,{x_result'},{x_sim},"AbsTol", 1.0e-14);
    verifyEqual(testCase,{y_result'},{y_sim},"AbsTol", 1.0e-14);
 
    
end

% Accesory Functions
function [y,x] = LegacyMsim(msys,u,t,x0)
%function [y,x] = msim(msys,u,t,x0)
% discrete-time simulation for multilinear CPN-1 models
% focus large scale, high effiency
% best if factor matrices are sparse logical 
% 
% gerwald.de
% 3.11.2022

% get dimensions
[n rf] = size(msys.F.phi);         % number of states
[p rg] = size(msys.G.phi);         % number of outputs
m = size(msys.F.U,1)-n;         % number of inputs 
if size(msys.G.U,1)-n ~= m 
    error('size mismatch') 
end
% Initialize result matrices
if length(t)==1
    T = t; 
    t = 1:T; 
else
    T = t(end);
end
if nargout>1 
    x = zeros(n,T);
else
    x = zeros(n,1);
end
x(:,1)= x0(:);
y = zeros(p,T);

% get factor matrix indices and values (only if double)
if islogical(msys.F.U) 
    [rfi,cfi] = find(msys.F.U);
    lfi = true;
else
    [rfi,cfi,vfi] = find(msys.F.U);
    sfi = -sign(vfi);  % sign vector of non-zero elements 
    lfi = false;
end
if islogical(msys.G.U)
    [rgi,cgi] = find(msys.G.U);
    lgi = true;
else
    [rgi,cgi,vgi] = find(msys.G.U);
    sgi = -sign(vgi);  % sign vector of non-zero elements
    lgi = false;
end

kx = 1;
kxp = 1;

% discrete time simulation
for k = t
    % index only needed for 2nd output arg (state trajectory)
    if nargout>1
        kx = k;
        kxp = k+1;
    end  
    % stack state-input vector
    xu = [x(:,kx);u(:,k)];
    % compute output 
    if lgi   % Boolean
        y(:,k) = msys.G.phi*accumarray(cgi,xu(rgi),[rg 1],@prod,1);
    else     % Double
        y(:,k) = msys.G.phi*accumarray(cgi,1+vgi.*(xu(rgi)+sgi),[rg 1],@prod,1);
    end
    % compute next state
    if lfi   % Boolean
        x(:,kxp) = msys.F.phi*accumarray(cfi,xu(rfi),[rf 1],@prod,1);
    else    % Double
        x(:,kxp) = msys.F.phi*accumarray(cfi,1+vfi.*(xu(rfi)+sfi),[rf 1],@prod,1);
    end
end
end

function sys = drmss(n,p,m,r,s,b)
%DRMSS  Generate random discrete-time multilinear state-space models.

if nargin < 1
   n=max([1,round(abs(10*randn(1,1)))]);
end
if nargin < 2
   p=1;
end
if nargin < 3   
   m=1;
end
if nargin < 4
   r=round(n+m)/2;
end
if nargin < 5  
    s = 0.2;
end
if nargin < 6
    b = false
end

sys.n = n;
sys.F.ntype = 1;
sys.F.U = sprand(n+m,r,s);
if b
    sys.F.U = sys.F.U>0;
end
sys.F.phi = sprand(n,r,0.7);

sys.G.ntype = 1;
sys.G.U = sprand(n+m,r,s);
if b
    sys.G.U = sys.G.U>0;
end
sys.G.phi = sprand(p,r,0.7);


end

