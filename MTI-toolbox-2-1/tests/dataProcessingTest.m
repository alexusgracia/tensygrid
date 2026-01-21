% Test header
function tests = dataProcessingTest
    tests = functiontests(localfunctions);
end

function calculateScaleOffsetTest(testCase)
% this borrows from mss2mss doc example, but is very much hard-coded
u = [zeros(20, 1); ones(20, 1); -ones(20, 1)]; % input signal
xsim = zeros(61, 2); % corresponding states
xsim(22:61, 1) = [0.1714, 0.3147, 0.4347, 0.5350, 0.6190, 0.6892, ...
    0.7480, 0.7972, 0.8383, 0.8727, 0.9015, 0.9256, 0.9458, 0.9626, ...
    0.9767, 0.9885, 0.9984, 1.0067, 1.0136, 1.0194, -1.0242, 0.6855, ...
    -0.7449, 0.4518, -0.5494, 0.2883, -0.4126, 0.1738, -0.3168, 0.0937, ...
    -0.2497, 0.0376, -0.2028, -0.0017, -0.1699, -0.0292, -0.1469, ... 
    -0.0484, -0.1308, -0.0619]; 
xsim(1:61, 2) = [0.0000, 0.5710, 0.9388, 1.1756, 1.3281, 1.4263, ...
    1.4895, 1.5302, 1.5564, 1.5733, 1.5842, 1.5912, 1.5957, 1.5986, ...
    1.6005, 1.6017, 1.6024, 1.6029, 1.6033, 1.6035, 1.6036, 1.4255, ...
    1.2984, 1.2070, 1.1406, 1.0915, 1.0545, 1.0261, 1.0039, 0.9864, ...
    0.9723, 0.9610, 0.9517, 0.9441, 0.9379, 0.9327, 0.9285, 0.9249, ...
    0.9220, 0.9195, 0.9175, 1.1448, 1.6190, 1.6563, 1.9702, 1.9415, ...
    2.1472, 2.1013, 2.2360, 2.1912, 2.2798, 2.2421, 2.3009, 2.2712, ...
    2.3105, 2.2881, 2.3145, 2.2980, 2.3159, 2.3039, 2.3161];
xsim2 = zeros(61, 2); % another input signal with all equal values
xsim2(:, 1) = 1.2916; 
xsim2(:, 2) = -0.0183;
xsim3 = zeros(0); % another (empty) input trajectory
lbx = [0, 0];                                  % lower limit state
ubx = [1, 1];                                  % upper limit state
lbu = 0;                                       % lower limit input
ubu = 1;                                       % upper limit input
lbwrong = [0, -1, -2]; % another lower bound that has the wrong length

% Test case 1: matrix with 2 trajectories, separate bounds
[scaleFactor1, offset1] = dataProcessing.calculateScaleOffset( ... 
    xsim, lbx, ubx);
scaleFactor1Expected = [2.0436, 2.3161];
offset1Expected = [-1.0242, 0]; 
verifyEqual(testCase, scaleFactor1, scaleFactor1Expected);
verifyEqual(testCase, offset1, offset1Expected); 

% Test case 2: one trajectory only
[scaleFactor2, offset2] = dataProcessing.calculateScaleOffset( ... 
    u, lbu, ubu);
scaleFactor2Expected = 2;
offset2Expected = -1; 
verifyEqual(testCase, scaleFactor2, scaleFactor2Expected);
verifyEqual(testCase, offset2, offset2Expected);

% Test case 3: matrix with 2 trajectories, same bound for all (scalar)
[scaleFactor3, offset3] = dataProcessing.calculateScaleOffset( ... 
    xsim, lbu, ubu);
scaleFactor3Expected = [2.0436, 2.3161];
offset3Expected = [-1.0242, 0]; 
verifyEqual(testCase, scaleFactor3, scaleFactor3Expected);
verifyEqual(testCase, offset3, offset3Expected);

% Test case 4: 2 trajectories, upper bound vector, lower bound scalar
[scaleFactor4, offset4] = dataProcessing.calculateScaleOffset( ... 
    xsim, lbx, ubu);
scaleFactor4Expected = [2.0436, 2.3161];
offset4Expected = [-1.0242, 0]; 
verifyEqual(testCase, scaleFactor4, scaleFactor4Expected);
verifyEqual(testCase, offset4, offset4Expected);

% Test case 5: only trajectory provided, no bounds
[scaleFactor5, offset5] = dataProcessing.calculateScaleOffset(xsim);
scaleFactor5Expected = [1, 1];
offset5Expected = [0, 0]; 
verifyEqual(testCase, scaleFactor5, scaleFactor5Expected);
verifyEqual(testCase, offset5, offset5Expected);

% Test case 6: trajectory contains all the same values
[scaleFactor6, offset6] = dataProcessing.calculateScaleOffset(...
    xsim2, lbx, ubu);
scaleFactor6Expected = [1, 1];
offset6Expected = [1.2916, -0.0183]; 
verifyEqual(testCase, scaleFactor6, scaleFactor6Expected);
verifyEqual(testCase, offset6, offset6Expected);

% Test case 7: trajectory contains all the same values and no bounds
[scaleFactor7, offset7] = dataProcessing.calculateScaleOffset(xsim2);
scaleFactor7Expected = [1, 1];
offset7Expected = [0, 0]; 
verifyEqual(testCase, scaleFactor7, scaleFactor7Expected);
verifyEqual(testCase, offset7, offset7Expected);

% Note that providing only bounds, but no trajectory, is not possible
% easily due to the limitations of MATLABs name-argument syntax which is 
% different to other languages (see comment in calculateScaleOffset). The
% only way to do it is to provide an empty trajectory
% Test case 8: only bounds provided, but no trajectory (empty)
[scaleFactor8, offset8] = dataProcessing.calculateScaleOffset(xsim3, ...
    lbu, ubx);
scaleFactor8Expected = zeros(0);
offset8Expected = zeros(0); 
verifyEqual(testCase, scaleFactor8, scaleFactor8Expected);
verifyEqual(testCase, offset8, offset8Expected);

% Test case 9: neither bounds nor trajectory provided
[scaleFactor9, offset9] = dataProcessing.calculateScaleOffset();
scaleFactor9Expected = zeros(0);
offset9Expected = zeros(0); 
verifyEqual(testCase, scaleFactor9, scaleFactor9Expected);
verifyEqual(testCase, offset9, offset9Expected);

% Test case 10: flipped lower and upper bounds
[scaleFactor10, offset10] = dataProcessing.calculateScaleOffset(xsim, ...
    ubx, lbu);
scaleFactor10Expected = [2.0436, 2.3161];
offset10Expected = [-1.0242, 0];  
verifyEqual(testCase, scaleFactor10, scaleFactor10Expected);
verifyEqual(testCase, offset10, offset10Expected);

% Test case 11: lower bound of wrong length 
[scaleFactor11, offset11] = dataProcessing.calculateScaleOffset(xsim, ...
    ubx, lbwrong);
scaleFactor11Expected = [2.0436, 2.3161];
offset11Expected = [-1.0242, 0];  
verifyEqual(testCase, scaleFactor11, scaleFactor11Expected);
verifyEqual(testCase, offset11, offset11Expected);

end 



function scaleDataTest(testCase)
% this borrows from mss2mss doc example, but is very much hard-coded.
% We do not test everything in as much detail because this is already done 
% in calculateScaleOffsetTest (scaleData calls calculateScaleOffset). 
u = [zeros(20, 1); ones(20, 1); -ones(20, 1)]; % input signal
xsim = zeros(61, 2); % corresponding states
xsim(22:61, 1) = [0.1714, 0.3147, 0.4347, 0.5350, 0.6190, 0.6892, ...
    0.7480, 0.7972, 0.8383, 0.8727, 0.9015, 0.9256, 0.9458, 0.9626, ...
    0.9767, 0.9885, 0.9984, 1.0067, 1.0136, 1.0194, -1.0242, 0.6855, ...
    -0.7449, 0.4518, -0.5494, 0.2883, -0.4126, 0.1738, -0.3168, 0.0937, ...
    -0.2497, 0.0376, -0.2028, -0.0017, -0.1699, -0.0292, -0.1469, ... 
    -0.0484, -0.1308, -0.0619]; 
xsim(1:61, 2) = [0.0000, 0.5710, 0.9388, 1.1756, 1.3281, 1.4263, ...
    1.4895, 1.5302, 1.5564, 1.5733, 1.5842, 1.5912, 1.5957, 1.5986, ...
    1.6005, 1.6017, 1.6024, 1.6029, 1.6033, 1.6035, 1.6036, 1.4255, ...
    1.2984, 1.2070, 1.1406, 1.0915, 1.0545, 1.0261, 1.0039, 0.9864, ...
    0.9723, 0.9610, 0.9517, 0.9441, 0.9379, 0.9327, 0.9285, 0.9249, ...
    0.9220, 0.9195, 0.9175, 1.1448, 1.6190, 1.6563, 1.9702, 1.9415, ...
    2.1472, 2.1013, 2.2360, 2.1912, 2.2798, 2.2421, 2.3009, 2.2712, ...
    2.3105, 2.2881, 2.3145, 2.2980, 2.3159, 2.3039, 2.3161];
lbx = [0, 0];                                  % lower limit state
ubx = [1, 1];                                  % upper limit state
lbu = 0;                                       % lower limit input
ubu = 1;                                       % upper limit input

% scale state only
[xsc, ~, T, c] = dataProcessing.scaleData(xsim, u, lbx, ubx); 
TExpected = [2.0436, 0, 0; 0, 2.3161, 0; 0, 0, 1];
cExpected = [-1.0242, 0, 0];
xscExpected = (xsim - c(1:2)) / T(1:2, 1:2);

verifyEqual(testCase, T, TExpected);
verifyEqual(testCase, c, cExpected);
verifyEqual(testCase, xsc, xscExpected, AbsTol = 1e-15);

% scale state and input
[xsc, usc, T, c] = dataProcessing.scaleData(xsim, u, lbx, ubx, lbu, ubu);
TExpected(3, 3) = 2;
cExpected(3) = -1;
uscExpected = (u - c(3)) / T(3, 3);

verifyEqual(testCase, T, TExpected);
verifyEqual(testCase, c, cExpected);
verifyEqual(testCase, xsc, xscExpected, AbsTol = 1e-15);
verifyEqual(testCase, usc, uscExpected, AbsTol = 1e-15);

end



function heavisideStepFunctionTest(testCase)
testInput = [-3 0 2];
testOutput = dataProcessing.heavisideStepFunction(testInput);
expectedOutput = [0 1 1];
verifyEqual(testCase, testOutput, expectedOutput);
end 
