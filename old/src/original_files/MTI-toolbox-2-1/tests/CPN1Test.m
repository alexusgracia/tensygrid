% Test header
function tests = CPN1Test
    tests = functiontests(localfunctions);
end

% emptytest prevents warning. Mathlab need at least 1 test function
% If you implement your own tests simply remove it
function fromFullTest(testCase)
    dataFile = "simMTITest" + "Fultensor.mat";
    F = load(dataFile).S;
    n = 2; % dataFile contains 2x2x2x2 tensor
    m = 1;
    cpn = CPN1(F);

    verifyEqual(testCase,ndims(F)-1,size(cpn.U,1));
    verifyEqual(testCase,n,size(cpn.phi,1));
    verifyEqual(testCase,m,size(cpn.U,1) - size(cpn.phi,1)); 
end

function fromNorm1CPNTest(testCase)
    dataFile = "simMTITest" + "Data.mat";
    sys = load(dataFile).sys;

    cpn = CPN1(sys.F.U,sys.F.phi);

    verifyEqual(testCase,sys.F.U,cpn.U);
    verifyEqual(testCase,sys.F.phi,cpn.phi);
end