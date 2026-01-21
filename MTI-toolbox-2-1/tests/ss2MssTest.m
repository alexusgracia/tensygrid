% Test header
function tests = ss2MssTest
    tests = functiontests(localfunctions);
end

function firstTest(testCase)
lsys=rss(6,2);
msys=mss.ss2Mss(lsys,'1');
verifyEqual(testCase,msys.F.phi(1:6,1:6),lsys.A);
end

