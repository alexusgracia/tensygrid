% Test header
function tests = ktensor2CPN1Test
    tests = functiontests(localfunctions);
end

function Cp2CpnShapeTest(testCase)
    r=6;
    n=2;
    m=1;
    F=ktensor(ones(r,1),ones(2,r),ones(2,r),ones(2,r),ones(2,r));
    [FU,Fphi]=CPN1.ktensor2CPN1(F);
    verifyEqual(testCase,size(FU),[n+m r]);
    verifyEqual(testCase,size(Fphi),[n r]);
end

function fromNegativeCPTest(testCase)
    load kTensTest;
    [FU,Fphi]=CPN1.ktensor2CPN1(kTens);
    verifyEqual(testCase,Fphi,FCpn.phi);
    verifyEqual(testCase,FU,FCpn.U);
end
function fromrandomNegativTest(testCase)
    load kTens_neg;
    [FU,Fphi]=CPN1.ktensor2CPN1(kTens);
    verifyEqual(testCase,Fphi,CpnTens.Phi);
    verifyEqual(testCase,FU,CpnTens.U);
end
function fromZeroFactorsCPTest(testCase)
    CPNexp.F=[1 0.5;
              0.5 0.5;
              0.5 0.5];
    CPNexp.Phi=[4 8;
                4 8];
               %lambda  x1              x2              u              phi 
    F=ktensor([1;1;1],[0 1 1;1 0 1],[1 0 1;1 0 1],[1 1 1;1 0 1],[1 1 1;1 1 1]);
    [FU,Fphi]=CPN1.ktensor2CPN1(F);
    verifyEqual(testCase,FU,CPNexp.F);
    verifyEqual(testCase,Fphi,CPNexp.Phi);
end