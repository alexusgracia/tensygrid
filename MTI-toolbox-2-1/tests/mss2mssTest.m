% Test header
function tests = mss2mssTest
    tests = functiontests(localfunctions);
end

function transformStateTest(testCase)
% this just basically goes through the example listed in the docs

% define expected output - hard-coded! (do this at the start so I can have 
% both cases with and without scaling input in one piece of code)
% not very nice code and does not follow the naming conventions etc, but 
% hopefully it should do the job... these tests are a mess anyway. 
u_exp = zeros(60, 1); 
u_exp(21:40) = 1; 
u_exp(41:60) = -1;

usc_exp = zeros(60, 1); 
usc_exp(1:20) = 0.5; 
usc_exp(21:40) = 1;

ysc_exp = zeros(60, 2);
ysc_exp(22:60, 1) = [0.1714, 0.3147, 0.4347, 0.5350, 0.6190, 0.6892, ...
    0.7480, 0.7972, 0.8383, 0.8727, 0.9015, 0.9256, 0.9458, 0.9626, ... 
    0.9767, 0.9885, 0.9984, 1.0067, 1.0136, 1.0194, -1.0242, 0.6855, ...
    -0.7449, 0.4518, -0.5494, 0.2883, -0.4126, 0.1738, -0.3168, 0.0937, ...
    -0.2497, 0.0376, -0.2028, -0.0017, -0.1699, -0.0292, -0.1469, ...
    -0.0484, -0.1308];
ysc_exp(2:60, 2) = [0.5711, 0.9388, 1.1756, 1.3281, 1.4263, 1.4895, ...
    1.5302, 1.5564, 1.5733, 1.5842, 1.5912, 1.5957, 1.5986, 1.6005, ...
    1.6017, 1.6024, 1.6029, 1.6033, 1.6035, 1.6036, 1.4255, 1.2984, ...
    1.2070, 1.1406, 1.0915, 1.0545, 1.0261, 1.0039, 0.9864, 0.9723, ... 
    0.9610, 0.9517, 0.9441, 0.9379, 0.9327, 0.9285, 0.9249, 0.9220, ... 
    0.9195, 0.9175, 1.1448, 1.6190, 1.6563, 1.9702, 1.9415, 2.1472, ...
    2.1013, 2.2360, 2.1912, 2.2798, 2.2421, 2.3009, 2.2712, 2.3105, ... 
    2.2881, 2.3145, 2.2980, 2.3159, 2.3039];
ysc_tol = zeros(60, 2) + 5e-5; % because I rounded values to 4 digits

xsimsc_exp = zeros(61, 2);
xsimsc_exp(1:61, 1) = [ones(1, 21)*0.5012, 0.5850, 0.6552, 0.7139, ...
    0.7630, 0.8041, 0.8385, 0.8672, 0.8913, 0.9114, 0.9282, 0.9423, ...
    0.9541, 0.9640, 0.9722, 0.9791, 0.9849, 0.9897, 0.9938, 0.9972, ...
    1.0000, -0.0000, 0.8366, 0.1367, 0.7223, 0.2323, 0.6423, 0.2993, ...
    0.5862, 0.3462, 0.5470, 0.3790, 0.5196, 0.4020, 0.5004, 0.4180, ...
    0.4869, 0.4293, 0.4775, 0.4372, 0.4709];
xsimsc_exp(2:61, 2) = [0.2466, 0.4053, 0.5076, 0.5734, 0.6158, 0.6431, ...
    0.6607, 0.6720, 0.6793, 0.6840, 0.6870, 0.6890, 0.6902, 0.6910, ...
    0.6915, 0.6919, 0.6921, 0.6922, 0.6923, 0.6924, 0.6155, 0.5606, ...
    0.5212, 0.4925, 0.4713, 0.4553, 0.4430, 0.4335, 0.4259, 0.4198, ...
    0.4149, 0.4109, 0.4076, 0.4049, 0.4027, 0.4009, 0.3993, 0.3981, ...
    0.3970, 0.3961, 0.4943, 0.6990, 0.7151, 0.8507, 0.8383, 0.9271, ...
    0.9073, 0.9654, 0.9461, 0.9843, 0.9681, 0.9934, 0.9806, 0.9976, ...
    0.9879, 0.9993, 0.9922, 0.9999, 0.9947, 1.0000];
xsimsc_tol = zeros(61, 2) + 5e-5; % different dimensions, so need again...

% Build multilinear model (MTI), 2.order, 1 input
F.U = [[0.83, -0.1]; [0, 0.53]; [1, -0.1]];    % state transition matrix
F.phi = [[1.0080, 0]; [0, 1.50]];              % state weigting matrix
u = [zeros(20, 1); ones(20, 1); -ones(20, 1)]; % input signal
t = 1:1:60;                                    % time vector
Ts = 1;                                        % sampling time (discrete)
obj = mss(CPN1(F.U, F.phi), Ts);               % build mti model object

% Simulate MTI Model
x0 = [0; 0];                                   % initial state
[y, ~, xsim] = msim(obj, u, t, x0);

% Set limits for states (inputs not scaled in this example)
lbx = [0, 0];                                  % lower limit state
ubx = [1, 1];                                  % upper limit state

% Case 1: Linear state transformation of mti model without scaling input
[xsc, ~, T, c] = dataProcessing.scaleData(xsim, u, lbx, ubx); 

% scale data within limits and obtain transformation matrix and offset
msys = mss2mss(obj, T, c);                 % transform (scale) system 
% Simulation of transformed mti model
[ysc, ~, xsimsc] = msim(msys, u, t, xsc(1,:)); 
% simulate scaled mti model with scaled initial state

verifyEqual(testCase, u, u_exp);
verifyEqual(testCase, ysc, ysc_exp, AbsTol = ysc_tol);
verifyEqual(testCase, xsimsc, xsimsc_exp, AbsTol = xsimsc_tol);

% Case 2: Linear state transformation of mti model with scaling input 
lbu = 0;                                       % lower limit input
ubu = 1;                                       % upper limit input
[xsc, usc, T, c] = dataProcessing.scaleData(xsim, u, lbx, ubx, lbu, ubu);
% scale data within limits and obtain transformation matrix and offset
msys = mss2mss(obj, T, c);                 % transform (scale) system
% Simulation of transformed mti model
[ysc, tOut, xsimsc] = msim(msys, usc, t, xsc(1,:)); 
% simulate scaled mti model with scaled initial state

verifyEqual(testCase, usc, usc_exp);
verifyEqual(testCase, ysc, ysc_exp, AbsTol = ysc_tol);
verifyEqual(testCase, xsimsc, xsimsc_exp, AbsTol = xsimsc_tol);


end

