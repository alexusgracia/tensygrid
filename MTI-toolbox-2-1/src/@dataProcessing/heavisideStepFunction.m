function [H] = heavisideStepFunction(x)
% <heavisideStepFunction> - The Heaviside step function H(x)
% is a step function which is zero for negative arguments
% of x and otherwise one. 
%   
%   input parameter:
%   - x: [input parameter vector].
%
%   output parameter:
%   - H: [Heaviside function evaluation vector]
%
%   Example:
%   H = heavisideStepFunction([-3 0 2]);
%   [Here H=[0 1 1]]
H = double(x >= 0);
end