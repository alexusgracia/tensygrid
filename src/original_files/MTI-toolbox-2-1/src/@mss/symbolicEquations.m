function [xp,y] = symbolicEquations(sys)
% <symbolicEquations> Creates a symbolic equation of the mss-object
%
% Note that this command requires the Symbolic Math Toolbox from Matlab. 
%
% Input arguments
%   - sys: mss object
% Output arguments:
%   - xp: symbolic equations of the state derivatives
%   - y : symbolic equations of the output equations
%
% Example:
% [xp,y] = symbolicEquations(sys)
%
% For detailed documentation see <a href="matlab:open((which('symbolicEquationsDoc.html')))">here</a>

% Enrico Uhlenberg, Marah Engels, Torben Warnecke, Leandro Samaniego,
% Carlos Cateriano Yáñez, Christoph Kaufmann, Leona Schnelle, Georg
% Pangalos, and Gerwald Lichtenberg - 12/06/2024

x = sym('x', [sys.n,1]);
u = sym('u', [sys.m,1]);
var=[x; u];

F_S=sys.F.U;
F_phi=sys.F.phi;

% state equations

eq=sym(zeros(length(F_phi(:,1)),1));
    for i=1:1:length(F_phi(:,1))
        for r=1:1:length(F_phi(1,:))
            monomial=1;
            for j=1:1:length(F_S(:,1))
                monomial=monomial*(1-abs(F_S(j,r))+F_S(j,r)*var(j));
            end
            eq(i)=eq(i)+F_phi(i,r)*monomial;
        end 
    end 
xp=eq;


% output equations

if sys.p~=0

G_S=sys.G.U;
G_phi=sys.G.phi;

eq=sym(zeros(length(G_phi(:,1)),1));
    for i=1:1:length(G_phi(:,1))
        for r=1:1:length(G_phi(1,:))
            monomial=1;
            for j=1:1:length(G_S(:,1))
                monomial=monomial*(1-abs(G_S(j,r))+G_S(j,r)*var(j));
            end
            eq(i)=eq(i)+G_phi(i,r)*monomial;
        end 
    end 
   y=eq;    
else
    y=[];
end 
 



end