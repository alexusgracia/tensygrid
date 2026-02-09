% Test header
function tests = OtiMTITest
    tests = functiontests(localfunctions);
end

% emptytest prevents warning. Mathlab need at least 1 test function
% If you implement your own tests simply remove it
function OptiMTInoErrorTest(testCase)
    dataFile = "testdata.mat";
    load(dataFile);
    n=2;
    m=1;
    r=3;
    W=nan(1,(2*n+m)*r);
    cost=cpnTens.OptiMTI(W,testdata.y,testdata.u,testdata.SamplingInstants,n,m,r,0,1,"simulation",0,"1-norm");
    verifyNotEmpty(testCase,cost)
end

function OptiMTIAutonomousTest(testCase)
end