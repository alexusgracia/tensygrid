% Test header
function tests = dmsimTest
    tests = functiontests(localfunctions);
end

function sys = setupOnce(testCase)
H = hyCPN1();

H.F.state = [0 0 1 0 0 0 0 0 0 1 0 0 0;...
    0 1 0 0 0 1 0 0 0 0 0 0 0;...
    0 0 0 0 0 0 0 0 0 0 1 0 0];
H.F.stateDerivative = [1 0 0 0 0 0 0 0 0 0 0 0 0;...
    0 0 0 1 0 0 0 0 0 0 0 0 0;...
    0 0 0 0 0 0 0 0 0 0 0 1 1];
H.F.algebraic = [0 0 0 0 0 0 1 0 0 0 0 0 0;...
    0 0 0 0 0 0 0 1 0 0 0 0 0];
H.F.boolean = [0 0 0 0 1 0 0 0 2 2 0 0 0];
H.F.input = [0 1 0 0 0 0 0 0 0 0 0 0 0;...
    0 0 1 0 0 0 0 0 0 0 0 0 0;...
    0 0 0 0 0 0 0 0 1 0 0 0 0];

H.phi.equality = [1 -0.13  0.13 0 0 0 0 0 0 0 0 0 0;...
    0 0 0 1 0 0 0 0 0 0 -1 0 0;...
    0 0 0 0 -1/1 1/1 0 0 0 0 15 0.5 0.5;...
    0 -1 0 0 0 0 1 0 0 0 0 0 0;...
    0 -1 0 0 0 0 1 0 0 0 0 0 0;...
    0 0 -1 0 0 0 0 1 0 0 0 0 0];
H.phi.inequality = [0 0 0 0 0 0 0 0 -1 1 0 0 0;...
    0 0 -1 0 0 0 0 1 0 0 0 0 0];

sys = dmss(H, 0);

testCase.TestData.dmss = sys;
end

function dmsimContinuousTimeTest(testCase)
sys = testCase.TestData.dmss;

ut = 0:150;
x0 = [0.01; 0; 0];
u = [1.2, 0.5, 0.8] .* ones(size(ut))';
u(50:end,3) = u(50:end,3) + 0.4;
opt = odeset('RelTol', 10.0^(-3),'AbsTol',10.0^(-7));

sim = dmsim(sys, x0, [], [], ut, u, opt);
verifyEqual(testCase, sim.tsim(end), ut(end));
verifyGreaterThan(testCase, sim.x(end,1), 1.2*(1-0.05));
verifyLessThan(testCase, sim.x(end,1), 1.2*(1+0.05));
end

function dmsimDiscontinuousTest(testCase)
sys = testCase.TestData.dmss;
sys = sys.c2d(0.1);

ut = 0:0.1:150;
x0 = [0.01; 0; 0];
u = [1.2, 0.5, 0.8] .* ones(size(ut))';
u(50:end,3) = u(50:end,3) + 0.4;

sim = dmsim(sys, x0, [], [], ut, u);
verifyEqual(testCase, sim.tsim(end), ut(end));
verifyLength(testCase, sim.tsim, length(ut));
verifyGreaterThan(testCase, sim.x(end,1), 1.2*(1-0.05));
verifyLessThan(testCase, sim.x(end,1), 1.2*(1+0.05));
end