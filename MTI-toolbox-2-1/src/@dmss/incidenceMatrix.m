function I = incidenceMatrix(sys)
%INCIDENCEMATRIX Calculate incidence matrix of the set of equations of a
%   dmss model.
%
%   The incidence matrix will be stored as a struct with entries for each
%   variable type. The rows of the sub matrices correspond to the different
%   equations.
%
% For detailed documentation see <a href="matlab:open((which('incidenceMatrixDoc.html')))">here</a>

% Torben Warnecke - 11/06/2024

% HFc = spones([sys.H.F.stateDerivative.c; sys.H.F.state.c; sys.H.F.input.c; sys.H.F.algebraic.c; sys.H.F.boolean.c]);
% HFt = [sys.H.F.stateDerivative.t; sys.H.F.state.t; sys.H.F.input.t; sys.H.F.algebraic.t; sys.H.F.boolean.t];
% HFf = [sys.H.F.stateDerivative.f; sys.H.F.state.f; sys.H.F.input.f; sys.H.F.algebraic.f; sys.H.F.boolean.f];

phiEq = spones(sys.H.phi.equality.c+sys.H.phi.equality.t-sys.H.phi.equality.f);
phiIneq = spones(sys.H.phi.inequality.c+sys.H.phi.inequality.t-sys.H.phi.inequality.f);

if sys.nEq>0
    I.equality.stateDerivative = spones(phiEq * spones([sys.H.F.stateDerivative.c])' + phiEq * spones([sys.H.F.stateDerivative.t])' + phiEq * spones([sys.H.F.stateDerivative.f])');
    I.equality.state = spones(phiEq * spones([sys.H.F.state.c])' + phiEq * spones([sys.H.F.state.t])' + phiEq * spones([sys.H.F.state.f])');
    I.equality.input = spones(phiEq * spones([sys.H.F.input.c])' + phiEq * spones([sys.H.F.input.t])' + phiEq * spones([sys.H.F.input.f])');
    I.equality.algebraic = spones(phiEq * spones([sys.H.F.algebraic.c])' + phiEq * spones([sys.H.F.algebraic.t])' + phiEq * spones([sys.H.F.algebraic.f])');
    I.equality.boolean = spones(phiEq * spones([sys.H.F.boolean.c])' + phiEq * spones([sys.H.F.boolean.t])' + phiEq * spones([sys.H.F.boolean.f])');
else
    I.equality.stateDerivative = double.empty([], sys.n);
    I.equality.state = double.empty([], sys.n);
    I.equality.input = double.empty([], sys.m);
    I.equality.algebraic = double.empty([], sys.p);
    I.equality.boolean = double.empty([], sys.q);
end
    
if sys.nIneq>0
    I.inequality.stateDerivative = spones(phiIneq * spones([sys.H.F.stateDerivative.c])' + phiIneq * spones([sys.H.F.stateDerivative.t])' + phiIneq * spones([sys.H.F.stateDerivative.f])');
    I.inequality.state = spones(phiIneq * spones([sys.H.F.state.c])' + phiIneq * spones([sys.H.F.state.t])' + phiIneq * spones([sys.H.F.state.f])');
    I.inequality.input = spones(phiIneq * spones([sys.H.F.input.c])' + phiIneq * spones([sys.H.F.input.t])' + phiIneq * spones([sys.H.F.input.f])');
    I.inequality.algebraic = spones(phiIneq * spones([sys.H.F.algebraic.c])' + phiIneq * spones([sys.H.F.algebraic.t])' + phiIneq * spones([sys.H.F.algebraic.f])');
    I.inequality.boolean = spones(phiIneq * spones([sys.H.F.boolean.c])' + phiIneq * spones([sys.H.F.boolean.t])' + phiIneq * spones([sys.H.F.boolean.f])');
else
    I.inequality.stateDerivative = double.empty([], sys.n);
    I.inequality.state = double.empty([], sys.n);
    I.inequality.input = double.empty([], sys.m);
    I.inequality.algebraic = double.empty([], sys.p);
    I.inequality.boolean = double.empty([], sys.q);
end

end

