function x = otvlsimA(otvlTensObj, ~, t, x0)
%OTVLSIMA Function performing discrete time state simulation of an MTI 
%system whose equation structure and parameters are given as an
%otvlTens-object. 
%
% inputs: 
%   otvlTensObj: otvlTens object containing multilinear state equations
%   t: time Vector
%   x0: initial state condition
%   
% output:
%   x: simulated state trajectory
%
% See also msim otvl otvlTens

% Marah Engels - 31/05/2024

% second function argument is a placeholder for inputs

if ~isa(otvlTensObj,'otvlTens')
    error(['Wrong class for input argument at first position. ' ...
        'Object has to be of class otvl'])
end

if ~(length(x0) == max(otvlTensObj.NVariables))
    error(['Dimension of initial condition not consistent with ' ...
        'number of variables'])
end

if size(t,1) > 1
    t = t'; 
end

FPhi = otvlTensObj.FPhi;
FaInter = otvlTensObj.FaInter;
Fc = otvlTensObj.Fc;

N_EQUATIONS = otvlTensObj.NEquations;
nRowsTotal = otvlTensObj.NRowsTotal; % sum of rows in all TVLs of the system
dc = otvlTensObj.DontCareArray; % merged Dont-Care-arrays of all TVLs
b = otvlTensObj.BooleanArray; % merged Boolean-arrays of all TVLs
timesteps = length(t);
xPre = zeros(timesteps,nRowsTotal); 
x = zeros(timesteps, N_EQUATIONS); % states
index = otvlTensObj.IndexVector; % index array defines groups for xPre -> x
x(1,:) = x0;

for ts = t(:,2:end) 
    xPre(ts,:) = (FPhi.*prod(dc + (~dc).*(b.*(x(ts-1,:) + FaInter(1,:)) ...
        +(~b).*(FaInter(2,:) - x(ts-1,:))),2))';
    x(ts,:) = accumarray(index, xPre(ts,:)')' + Fc';
end
