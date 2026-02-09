% Test header
function tests = rmssTest
    tests = functiontests(localfunctions);
end

function firstTest(testCase)
% Simple test if rmss is executable
msys=mss.rmss(4,2,1,5);
[n,r]=size(msys.F.phi);
verifyEqual(testCase,[n r],[4 5]);
end

