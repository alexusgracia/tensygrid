function sys = sym2dmss(symfun, Ts, symbols)
% SYM2DMSS - Converts symbolic or string equations into dmss object.
%   
%   input parameter:
%   - symfun: symbolic or string array/cell with polynomial equations 
%   only containing integer exponents (auxilary variables and equations
%   will be used to represent integer exponents not equal to 1).
%   Math operators must be used in between each variable and parameter. 
%   Possible operators are "*", "+", "-", "/", "^", as well as brackets "("
%   and ")". For different variables of the same type use integers 
%   after the variable type symbol, e.g. for states "x1" and "x2". 
%   If an integer is missing, the leftover numbers will be counted as
%   parameters. "e+" and "e-" is understood as scientific format, e.g. in "3.0e+07".
%   - Ts: stepsize of the model (ts = 0: continuous-time,
%   ts > 0: discrete-time)
%   - symbol (optional): Symbols used for the different variable types.
%   Must be string- or symbolic-array.
%   Default is ["xp", "x", "u", "y", "z"]
%
%   output parameter:
%   - sys: the model as a dmss-object
%
%   Example:
%   eq = {"xp1 - 0.13*u1*u2 + 0.13*x1*u3 = 0";...
%   "y1 - u1*u2 = 0";...
%   "y2 - u3*x1 = 0"};
%
%   ts = 0;
%
%   sys = sym2dmss(eq, ts);
%
%   Note:
%   Also equations like "u1*f*(xp1-x1)/(u2*y1-k*x1*y2^2)  = u1*u3^(-1)"
%   are possible, it is only important that they can be converted into an
%   implicit multilinear equation.
%
% For detailed documentation see <a href="matlab:open((which('sym2dmssDoc.html')))">here</a>

% Torben Warnecke - 11/06/2024

arguments
    symfun (:,1) string
    Ts (1,1) double
    symbols (1,:) string = ["xp", "x", "u", "y", "z"] 
    % Can also also a sym-array e.g. [sym("dx"), sym("x"), sym("in"), sym("out"), sym("binary")] 
    % and will be converted into a string array by the input parser.
    % Also just parsing 4 symbols is possible, if dmss-model is not hybrid.
end

%str_array = string(symfun);

str_cell = cell(length(symfun),1);
for k = 1:length(symfun)%k = 1:length(str_array)
    str_cell{k,1} = symfun(k);%str_cell{k,1} = str_array(k);
end

%try
    sys = dmss.poly2dmss(str_cell, Ts, symbols);
% catch
%     error('An error occured: Maybe multiplication signs are missing, they must be written.')
% end

end