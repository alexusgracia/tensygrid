% Test header
function tests = mssTest
    tests = functiontests(localfunctions);
end

function constructFromNormOneOnlyFTest(testCase)
    F = 1;
    F_Phi = 2;
    sys = mss(CPN1(F,F_Phi));
    
    verifyInstanceOf(testCase,sys,'mss');
    verifyEqual(testCase,sys.F.U,F);
    verifyEqual(testCase,sys.F.phi,F_Phi);

end

function nmChecktest(testCase)
    dataFile = "simMTITest" + "Fultensor.mat";
    F = load(dataFile).S;
    n = 2; % dataFile contains 2x2x2x2 tensor
    m = 1;

    msys = mss(CPN1(F));

    verifyEqual(testCase,msys.n,n);
    verifyEqual(testCase,msys.m,m);
    
end

function  constructFromMatrixTest(testCase)
    F = [0.5 0.32 1 -0.27];
    G = [0 1 0 0; 0 0 -2 0];
    sys = mss(F, G);

    verifyEqual(testCase,sys.n,1);
    verifyEqual(testCase,sys.m,1);
    verifyEqual(testCase,sys.p,2);
    verifyEqual(testCase,sys.G.U,[0 1 0 1; 0 0 1 1]);
    verifyEqual(testCase,sys.G.phi(2,3),-2);
    
    F = [1 -1];
    G = [1 2; 3 4];
    verifyWarning(testCase, @() mss(F, G), 'CPN1:matrix2x2Warning');
end