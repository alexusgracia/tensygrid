function sfunc_mss(block)
% msstest level-2 s-func for test to add an mss object
% gerwald.de
% 22.8.25 

setup(block);
  
function setup(block)
  % Parameters
  block.NumDialogPrms = 2;   % 1st: model, 2nd: initial state
  block.DialogPrmsTunable = {'Nontunable','Nontunable'};
  msys = block.DialogPrm(1).Data;             % get the model

  %  Input ports
  block.NumInputPorts  = 1;                   % only one input port 
  block.SetPreCompInpPortInfoToDynamic;       % dynamic port info
  block.InputPort(1).Dimensions = msys.m;     % dimension input vector
  block.InputPort(1).DatatypeID  = 0;         % double
  block.InputPort(1).Complexity  = 'Real';    % real
  
  % Output ports
  block.NumOutputPorts = 1;                   % 2 output ports for output and states to be used in Simulink)
  block.SetPreCompOutPortInfoToDynamic;       % dynamic port info
  block.OutputPort(1).Dimensions = msys.p;    % dimension output vector
  block.OutputPort(1).DatatypeID  = 0;        % double
  block.OutputPort(1).Complexity  = 'Real';   % real
  
  block.SampleTimes = [msys.ts 0];            % set sampling times
  if msys.ts      
     block.NumContStates = 0;                 %  discrete-time model
  else
     block.NumContStates = msys.n;            %  continuous-time model
  end

  % Define subfunctions to call
  block.SetAccelRunOnTLC(false);
  block.RegBlockMethod('InitializeConditions', @InitializeConditions);
  block.RegBlockMethod('Outputs', @Outputs);
  block.RegBlockMethod('Update', @Update);
  block.RegBlockMethod('Derivatives', @Derivatives);
  block.RegBlockMethod('PostPropagationSetup',@PostPropagationSetup);
%endfunction
    
function PostPropagationSetup(block)
  msys = block.DialogPrm(1).Data;            % get the model
  if msys.ts  
    block.NumDworks = 1;                     % discete-time state vector
    block.Dwork(1).Name            = 'x';    % state name default
    block.Dwork(1).Dimensions      = msys.n; % order
    block.Dwork(1).DatatypeID      = 0;      % double 
    block.Dwork(1).Complexity      = 'Real'; % real
    block.Dwork(1).UsedAsDiscState = true;   
  end

function InitializeConditions(block)
   x0 = block.DialogPrm(2).Data(:);  % initial conditions from 2nd arg
      if block.SampleTimes(1)       
       block.Dwork(1).Data=x0;       % initialize discrete-time state
   else
       block.ContStates.Data=x0;     % initialize continuous-time state
   end
%endfunction

function Outputs(block)
   % stack state-input column vector
   if block.SampleTimes(1)  % discrete-time model 
      xu = [block.Dwork(1).Data(:); block.InputPort(1).Data(:)];
      block.OutputPort(1).Data = [block.DialogPrm(1).Data.getOutput(xu)];
   else
      xu = [block.ContStates.Data(:); block.InputPort(1).Data(:)];
      block.OutputPort(1).Data = [block.DialogPrm(1).Data.getOutput(xu)];
   end
   %endfunction

function Update(block)
   if block.SampleTimes(1)
        % stack state-input column vector
        xu = [block.Dwork(1).Data; block.InputPort(1).Data];
        % compute the discrete-time next state (MTI Toolbox method)
        block.Dwork(1).Data = block.DialogPrm(1).Data.getNextState(xu);
   end
%endfunction

function Derivatives(block)
    if ~block.SampleTimes(1)
        % stack state-input column vector
        xu = [block.ContStates.Data(:);block.InputPort(1).Data(:)];
        % compute the continuous-time state derivative (MTI Toolbox method)
        block.Derivatives.Data = block.DialogPrm(1).Data.getNextState(xu);
    end
%endfunction