% Test header
function tests = dmssTest
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

function dmssConstructorTest(testCase)
emptySys= dmss();

sys = testCase.TestData.dmss;
verifyClass(testCase, sys.H, 'hyCPN1');
verifyEqual(testCase, sys.H.F.algebraic.t(1,7), sparse(true));
end

function mss2dmssConstructorTest(testCase)
F = CPN1([1 0 1 1; 0 1 0 1; 0 0 1 1], [1 -1.2 3 2]);
G = CPN1([1 1 0; 0 1 1; 0 0 1], [1 0 0; 0 1 -0.5]);
eSys = mss(F,G, 2);

iSys = dmss(eSys);

%disp(iSys.symbolicEquations)

verifyEqual(testCase, iSys.n, eSys.n);
verifyEqual(testCase, iSys.m, eSys.m);
verifyEqual(testCase, iSys.p, eSys.p);
verifyEqual(testCase, iSys.ts, eSys.ts);
verifyEqual(testCase, full(iSys.H.phi.equality.c+iSys.H.phi.equality.t-iSys.H.phi.equality.f), [-1 0 0 1 -1.2 3 2 0 0; 0 -1 0 1 0 0 0 0 0; 0 0 -1 0 0 0 0 1 -0.5]);
end

function mss2dmssMethodOfMssTest(testCase)
F = CPN1([1 0 1 1; 0 1 0 1; 0 0 1 1], [1 -1.2 3 2]);
G = CPN1([1 1 0; 0 1 1; 0 0 1], [1 0 0; 0 1 -0.5]);
eSys = mss(F,G, 2);

iSys = eSys.mss2dmss;

%disp(iSys.symbolicEquations)

verifyEqual(testCase, iSys.n, eSys.n);
verifyEqual(testCase, iSys.m, eSys.m);
verifyEqual(testCase, iSys.p, eSys.p);
verifyEqual(testCase, iSys.ts, eSys.ts);
verifyEqual(testCase, full(iSys.H.phi.equality.c+iSys.H.phi.equality.t-iSys.H.phi.equality.f), [-1 0 0 1 -1.2 3 2 0 0; 0 -1 0 1 0 0 0 0 0; 0 0 -1 0 0 0 0 1 -0.5]);

iSys2 = mss2dmss(eSys);

verifyEqual(testCase, iSys2.n, eSys.n);
verifyEqual(testCase, iSys2.m, eSys.m);
verifyEqual(testCase, iSys2.p, eSys.p);
verifyEqual(testCase, iSys2.ts, eSys.ts);
verifyEqual(testCase, full(iSys.H.phi.equality.c+iSys.H.phi.equality.t-iSys.H.phi.equality.f), [-1 0 0 1 -1.2 3 2 0 0; 0 -1 0 1 0 0 0 0 0; 0 0 -1 0 0 0 0 1 -0.5]);
end

