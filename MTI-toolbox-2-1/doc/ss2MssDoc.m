%% <ss2Mss> - Takes a ss object and translates it to an mss object
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