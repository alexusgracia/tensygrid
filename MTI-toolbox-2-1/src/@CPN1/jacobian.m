function [J] = jacobian(obj,op)
% <jacobian> - Computes Jacobian around operating point
%
% Input parameter: 
%  - obj: CPN1 object 
%  - x  : operating point 
% 
% Output parameter:
%  - J  : Jacobian matrix
%
% Example:
% J = jacobian(obj,x_op)
%
% For detailed documentation see <a href="matlab:open((which('cpn2LinDoc.html')))">here</a>

% Enrico Uhlenberg, Marah Engels, Torben Warnecke, Leandro Samaniego,
% Carlos Cateriano Yáñez, Christoph Kaufmann, Leona Schnelle, Georg
% Pangalos, and Gerwald Lichtenberg - 12/06/2024


    if issparse(obj.U)
         J=cpn2LinSparse(obj,op);
    else 
         % Jacobian of states    
            J = cpn2Lin(obj,op);
    end

end

