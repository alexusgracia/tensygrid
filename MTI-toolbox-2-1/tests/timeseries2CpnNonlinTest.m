% Test header
function tests = timeseries2CpnNonlinTest
    tests = functiontests(localfunctions);
end

% emptytest prevents warning. Mathlab need at least 1 test function
% If you implement your own tests simply remove it
function ParameterReshapingTest(testCase)
    dataFile = "testdata.mat";
    load(dataFile);
    n=2;
    m=1;
    r=3;
    [sys]=mlgreyest(testdata,r,Method="fmincon",Initialize="zero");
    verifySize(testCase,sys.F.U,[n+m r])
    verifySize(testCase,sys.F.phi,[n r])
end
function ParameterOptionTest(testCase)
    dataFile = "testdata.mat";
    load(dataFile);
    r=1;
    [sys]=mlgreyest(testdata,r,upperBound=0,Method="fmincon",Initialize="ga");
    verifyLessThanOrEqual(testCase,sys.F.U,1)
end

function ParameterIdentificationAutonomousTest(testCase)
end