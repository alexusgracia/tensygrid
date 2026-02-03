% Test header
function tests = drmssTest
    tests = functiontests(localfunctions);
end

function firstTest(testCase)
% Simple test if drmss is executable
dmsys=mss.drmss(4,2,1,5,1e-3);
[n,r]=size(dmsys.F.phi);
verifyEqual(testCase,[n r],[4 5]);
end

