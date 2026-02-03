% Test header
function tests = cpnTensTest
    tests = functiontests(localfunctions);
end

function createRandomCpnTensTest(testCase)
n = 5;
p = 3;
m = 4;
expectedBooleanU = false;
expectednormtype = class(CPN1());
[FU,FPHI,GU,GPHI] = cpnTens.randCpn(n,p,m);
msys = mss(CPN1(FU,FPHI),CPN1(GU,GPHI));
verifyEqual(testCase,[n m p expectedBooleanU],[msys.n msys.m msys.p isa(msys.F.U,'logical')] );
verifyClass(testCase,msys.F,expectednormtype);
end
