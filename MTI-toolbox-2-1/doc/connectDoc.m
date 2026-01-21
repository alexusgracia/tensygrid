%% connect
% Connect signals of dmss model(s)
%% Syntax
% sys = connect(sys)
%   Connect an appended system based on signal names
%
% sys = connect(sys, stateConnections, algebraicConnections, booleanConnections, similarInputs)
%   Connect an appended system based on indexes (state-input-/
%   algebraic-inputs-/ boolean-input-connections and similar
%   inputs). Indexes are parsed as connection matrices with input indexes
%   in the first column and the variable indices in the following columns [inputIdx, variableIdx, ...; ...].
%   If more then 2 values are in a row a additional (linear) connection equation will be added: u = var(1) + var(2) + ...
%   If a variable index is parsed as negative value it is subtracted in the
%   connection equation: u = -var(1) ...
%
% sys = connect(sys1, sys2, ..., sysN):
%   Connect multiple systems based on signal names
%
%% Examples
% *Connect via index matrices: simplified HVAC system*
% 
% The supplied room of an HVAC system can be modelled by its power and
% transport balance, as well as the enthalpy equation for the humid air.
powerBalance = ["0 = -m*dx1 + u1 + u2"];
transportBalance = ["0 = -m*dx2 + u3"];
enthalpy = ["x1 = c_d * y1 + c_v * y1 * x2 + r0 * x2"];
room = sym2dmss([powerBalance; transportBalance; enthalpy], 0, ["dx", "x", "u", "y"]);
%%
% where x1 is the room air enthalpy, x2 is the absolute room air humidity,
% dx1 and dx2 are their derivatives, y1 is the room temperature, u1 is the
% heat transfered by the air flow of the air conditioning system (either heating or cooling the room air), u2 is the
% heat loss to the environment and u3 is the water vapor mass flow
% transfered by the air flow of the air conditioning system (either humidifying
% or dehumidifying the room air).
%
% The energy and humidity transfer by the air conditioning system can be
% modelled by the following equations
airConditioning = ["y1 = c_d*u1*u2 + c_v*u1*u2*u3 + r0*u1*u3 - u1*u4";...
                    "y2 = u1*u3 - u1*u5"];
airConditioning = sym2dmss(airConditioning, 0, ["dx", "x", "u", "y"]);
%%
% where y1 corresponds to the heat transfered by the air flow of the air
% conditioning system, y2 to the water vapor mass flow transfered by the air
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
%%
% *Connect via signal names: simplified HVAC system*
%
% Anopther way of connecting the sub-models is via signal names. 
% Therefor the signal names need to be set accordingly.
% For the naming it is important that there are no duplicats in the
% states and algebraic variable names. If there are duplicats in the input
% names, they will be set as similar inputs.
room.stateName = ["room enthalpy", "room humidity"];
room.algebraicName = ["room temperature"];
room.inputName = ["heat transfer by air conditioning", "heat loss", "humidity transfer by air conditioning"];

airConditioning.algebraicName = ["heat transfer by air conditioning", "humidity transfer by air conditioning"];
airConditioning.inputName = ["air flow", "supply temperature", "supply humidity", "room enthalpy", "room humidity"];

loss.algebraicName = ["heat loss"];
loss.inputName = ["ambient temperature", "room temperature"];

sys = connect(room, airConditioning, loss);

disp(sys.symbolicEquations)

%% References
%
%% See also
% <dmssDoc.html dmss>,
% <appendDoc.html append>,
% <sym2dmssDoc.html sym2dmss>,
% <replaceSymbolicParametersDoc.html replaceSymbolicParameters>
%
%
% Author(s): Torben Warnecke