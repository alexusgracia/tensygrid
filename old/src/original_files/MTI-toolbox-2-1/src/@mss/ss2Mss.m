function msys = ss2Mss(lsys,normtype)
% <ss2Mss> - Takes a ss object and translates it to an mss object
%   
%   Input parameter    
%       - lsys:  lss-object, linear state-space model of type 'ss'
%       - normtype - type of CPN norm. Currently only the 1-norm CPN is
%           supported
%
%   Output parameter
%        -msys: mss object - multilinear time-invariant state space model
%
%   Example: mysys=mss.ss2Mss(lsys,'1')


    %check if ss object
    if ~isa(lsys,'ss')
        error("Input is not state-space model of type 'ss'. Please change to 'ss'.")
    end
    
    %check if 1-norm CPN 
    if (normtype ~='1')
        error("Normtype not yet supported. Please use '1'")
    end
    
    nm=sum(size(lsys.B)); % Number of states and input
    nq=sum(size(lsys.D)); % Number of states and feedforward terms
 
    %Convert into mss-object and pass on discretization step size
    msys=mss(CPN1(eye(nm),[lsys.A, lsys.B]),CPN1(eye(nq),[lsys.C, lsys.D]),lsys.Ts);
    
    % TODO: pass on names of variables 
end

