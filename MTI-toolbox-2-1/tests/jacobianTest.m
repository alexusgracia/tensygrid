% Test header
function tests = jacobianTest
    tests = functiontests(localfunctions);
end


function firstTest(testCase)
% Simple test for the cpn2lin function for 1normCpn
U=[0.5 -0.5; 1 0; 0 1];
phi=[2 0; 4 6; 0 8];
obj=CPN1(U,phi);
x_op=[-1;1;3]; % to get a zero entry

Jtrue=[1 0 0; -7 0 6; -12 0 8];
J=jacobian(obj,x_op);
verifyEqual(testCase,J,Jtrue);
end

