function msys = drmss(n,m,p,r,Ts)
% <rmss> - Generates a random discrete-time multilinear state-space model
%
% Input parameter:
%   n - number of states, defaults to random number between 1 and
%       a positive random integer drawn from a standard normal distribution
%       with standard deviation 10
%   m - number of inputs, , defaults to 1
%   p - number of outputs, defaults to 1
%   r - rank of the resulting decomposed tensor, defaults to
%          the nearest integer of (n+m)/2
%   Ts - discretization time  
%
% Output parameter:
%  msys - mss object  
%
% Example:  
% msys=drmss(4,1,1,6,1e-3)
% For a system with n=4 states, m=1 input, p=1 output, a rank of r=6, 
% and discretization time of Ts=1e-3 seconds.           

if nargin < 1
    warning ('No number of state variables specified. Random number of state variables will be taken.')
    n=max([1,round(abs(10*randn(1,1)))]);
end

if nargin < 2
    warning ('No number of input variables m specified. Default is set to 1.')
    m=1;
end

if nargin < 3
    warning ('No number of output variables p specified. Default is set to 1.')
    p=1;
end

if nargin < 4
    warning ('No rank is specified. Default is set to the nearest integer of (n+m)/2.')
    r=round((n+m)/2);
end

if nargin < 5
    warning ('No discretization time is specified. Default is set to 1e-4 seconds.')
    Ts=1e-4;
end

[F_U,F_phi,G_U,G_phi]=cpnTens.randCpn(n,p,m,r,false,'1');
msys=mss(CPN1(F_U,F_phi),[CPN1(G_U,G_phi)],Ts);
end

