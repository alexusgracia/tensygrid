function J = cpn2Lin(obj,x)
% <cpn2Lin> - Computes Jacobian around operating point
%
% Input parameter: 
%  - obj: CPN1 object 
%  - x  : operating point 
% 
% Output parameter:
%  - J  : Jacobian matrix
%
% Example:
% J = cpn2Lin(obj,x_op)
%
% For detailed documentation see <a href="matlab:open((which('cpn2LinDoc.html')))">here</a>

% Enrico Uhlenberg, Marah Engels, Torben Warnecke, Leandro Samaniego,
% Carlos Cateriano Yáñez, Christoph Kaufmann, Leona Schnelle, Georg
% Pangalos, and Gerwald Lichtenberg - 12/06/2024

%% Dimensionality check 
if size(obj.U,1) ~= length(x)
    error("Dimensions of the operating point argument x do not match with the CPN1 object dimensions.")
end
  
%% Linearization procedure 

temp = 1-abs(obj.U)+obj.U.*x;   % temporary matrix to check for zeros 
E = ~temp;                      % zero entries
temp(E) = 1;                    % zeplace zero entries by ones
fp = prod(temp);                % compute products of all nonzeros factors
Fd = obj.U./temp;               % compute individual factor  
sE = sum(E);                    % count the number of zero entries 
% Compute Jacobian (which is only nonzero iff at most one factor is zero)
J = obj.phi*((Fd.*fp).*~(~E.*(sE==1)).*(sE<=1)).';  

end % of cpn2Lin

