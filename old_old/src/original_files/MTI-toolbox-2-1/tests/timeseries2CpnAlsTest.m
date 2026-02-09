% Test header
function tests = timeseries2CpnAlsTest
    tests = functiontests(localfunctions);
end

% emptytest prevents warning. Mathlab need at least 1 test function
% If you implement your own tests simply remove it
function timeseries2CpnAlsnoErrorTest(testCase)
    dataFile = "testdata.mat";
    load(dataFile);
    n=2;
    m=1;
    r=1;
    sys=mlgreyest(testdata,r,"AlsTolerance",10);
    verifyNotEmpty(testCase,sys.F.U,sys.F.phi)
end

function timeseries2CpnAlsSizeTest(testCase)
    dataFile = "testdata.mat";
    load(dataFile);
    n=2;
    m=1;
    r=2;
    [sys]=mlgreyest(testdata,r,"AlsIterations",10);
    verifySize(testCase,[sys.F.U;sys.F.phi],[n+m+n,r])
end
