% Test header
function tests = linearizeTest
    tests = functiontests(localfunctions);
end



function linearizeMSStest(testCase)
% Simple test for the cpn2lin function for 1normCpn
U=[0.5 -0.5; 1 0; 0 1];
phi=[2 0; 4 6];
x_op=[-1;1]; % to get a zero entry
u_op=3;

sys=mss(CPN1(U,phi));
linsys=linearize(sys,x_op,u_op);


Atrue=[1 0; -7 0]; 
Btrue=[0;6]; 
verifyEqual(testCase,[linsys.A,linsys.B],[Atrue, Btrue]);
end

function linearizeDMSStest(testCase)
% Test for linearization of a dmss object in CPN1 format

% Create iMTI model as dmss object
eqs={'xp1==x1*x2^2-4';
    'xp2==x2-2';};

iMTImodel=sym2dmss(eqs,0);

% Operating point
xpeq=[0 0].';
xeq=[1 2]';

% Linearizing iMTI model
descriptorModel=linearize(iMTImodel,xpeq,xeq,[],xeq(2),[]);


% Compute generalized eigenvalues of the obtained descriptor model
 if isMATLABReleaseOlderThan('R2025a')
        generalizedEigenvaluesDescriptorModel=eig(full(descriptorModel));
 else 
        generalizedEigenvaluesDescriptorModel=eig(descriptorModel);
 end

%Nonlinear reference
A=[xeq(2)^2  2*xeq(1)*xeq(2);
    0                   1];
eigenvaluesLinearModel=eig(A);


verifyEqual(testCase,generalizedEigenvaluesDescriptorModel,eigenvaluesLinearModel)

end


function thirdTest(testCase)
    % check the operation with op 
	% Construct small mss example
	U=[0.5 -0.5; 1 0; 0 1];
	phi=[2 0; 4 6];
	x_op=[-1;1];
	u_op=3;
	sys=mss(CPN1(U,phi));

	op = struct('x', x_op, 'u', u_op);
	linsys = linearize(sys, op); % struct-based wrapper

	% Compare with legacy call for consistency
	linsysLegacy = linearize(sys, x_op, u_op);
	verifyEqual(testCase, linsys.A, linsysLegacy.A);
	verifyEqual(testCase, linsys.B, linsysLegacy.B);
	verifyEqual(testCase, linsys.C, linsysLegacy.C);
	verifyEqual(testCase, linsys.D, linsysLegacy.D);
end

function fourthTest(testCase)
    % check the operation with op 
    
	% Simple implicit system from existing test
	eqs={'xp1==x1*x2^2-4';
		 'xp2==x2-2'};
	dsys=sym2dmss(eqs,0);

	xpeq=[0 0].';
	xeq=[1 2].';

	% struct with only xp and x
	op = struct('xp', xpeq, 'x', xeq, 'y',xeq(2));
	dss = linearize(dsys, op);

	% Compare generalized eigenvalues with legacy call
	dssLegacy = linearize(dsys, xpeq, xeq, [], xeq(2), []);
    if isMATLABReleaseOlderThan('R2025a')
	    verifyEqual(testCase, eig(full(dss)), eig(full(dssLegacy)));
    else 
        verifyEqual(testCase, eig(dss), eig(dssLegacy));
    end 
end
