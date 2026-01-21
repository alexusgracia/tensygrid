%% append
% Appends multiple dmss models (diagonally).
%% Syntax
% sys = sys1.append(sys2, sys3, ..., sysn) Appends multiple dmss models in
% a diagonal fashion.
%% Description
% The decomposed tensor matrices (structure matrix sys.H.F and parameter
% matrix sys.H.phi) will be appended diagonally, as example:
%
% sys.H.F.state = [sys1.H.F.state, 0; 0, sys2.H.F.state] -> similar with
% the other sub matrices.
%
% sys.H.phi.equality = [sys1.H.phi.equality, 0; 0, sys2.H.phi.equality]
%
% The number of variables of the resulting system is the sum of the number
% of variables of the appended systems.
% To connect systems via variable index matrices, they must be appended
% first.
%% Example
% *Appending component models of a (simplified) HVAC system*
% 
% The supplied room of an HVAC system can be modelled by its power and
% transport balance, as well as the enthalpy equation for the humid air.
powerBalance = ["0 = -m*dh1 + p1 + p2"];
transportBalance = ["0 = -m*dh2 + p3"];
enthalpy = ["h1 = c_d * T1 + c_v * T1 * h2 + r0 * h2"];
room = sym2dmss([powerBalance; transportBalance; enthalpy], 0, ["dh", "h", "p", "T"]);
%%
% where h1 is the room air enthalpy, h2 is the absolute room air humidity, dh1 and dh2 are their derivatives, T1 is the room temperature, p1 is the heat
% transported by the air flow of the air conditioning system (either heating or cooling the room air), p2 is the
% heat loss to the environment and p3 is the water vapor mass flow
% transported by the air flow of the air conditioning system (either humidifying
% or dehumidifying the room air).
%
% The energy and humidity transfer by the air conditioning system can be
% modelled by the following equations
airConditioning = ["y1 = c_d*u1*u2 + c_v*u1*u2*u3 + r0*u1*u3 - u1*u4";...
                    "y2 = u1*u3 - u1*u5"];
airConditioning = sym2dmss(airConditioning, 0, ["dx", "x", "u", "y"]);
%%
% where y1 corresponds to the power transported by the air flow of the air
% conditioning system, y2 to the water vapor mass flow transported by the air
% flow of the air conditioning system, u1 to the supplied mass flow of air, u2
% to the supply temperature, u3 to the supply humidity, u4 to the rooms
% enthalpy and u5 to the rooms absolute humidity.
%
% Additionally the heating loss to the environment can be calculated by
loss = ["y1 = k*A*(u1-u2)"];
loss = sym2dmss(loss, 0, ["dx", "x", "u", "y"]);
%%
% where y1 is the heat loss, u1 is the ambient temperature and u2 is the
% rooms temperature.
%
% The component models can be appended to obtain an overall model of the
% system
sys = append(room, airConditioning, loss);
disp(sys.symbolicEquations)
%%
% Due to the appending the number of variables as well
% as their index is increasing accordingly. 
% 
% The input and output signals can be connected by using connection matrices, where the first column corresponds to the input
% index and the second column to the state/algebraic signal index.
stateInputConnection = [7 1; 8 2];
algebraicInputConnection = [9 1; 1 2; 3 3; 2 4];
similarInputs = [];
sys = sys.connect(stateInputConnection, algebraicInputConnection, [], similarInputs);

disp(sys.symbolicEquations)

%% References
%
%% See also
% <dmssDoc.html dmss>,
% <connectDoc.html connect>,
% <hyCPN1Doc.html hyCPN1>,
% <sym2dmssDoc.html sym2dmss>
%
%
% Author(s): Torben Warnecke