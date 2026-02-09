function [y,tOut,x] = msimDisc(sim, u, t, x0)
% MSIMDISC Simulate discrete-time time response of dynamic discrete-time explicit MTI model to arbitrary inputs. 

n = sim.sys.n;
p = sim.sys.p;
m = sim.sys.m;         % number of inputs

if ~isempty(sim.sys.G) && size(sim.sys.G.U,1)-n ~= m
    error('size mismatch')
end
% if(nargout>2 && p ==0)
%     warning('System does not have output tensor, the output trajecory will equal the state trajectory')
% end

% We need the state trajectory if either:
% - the user asked for state and output trajectories
% - There is no output equation, so only the state trajectory is of
%   interest (and re'labeled' as y at the end of the routine)
includeStateTrajectory = (nargout > 2) || (p == 0); 

% Initialize result matrices
if length(t)==1
    t = sim.sys.ts:sim.sys.ts:t;
end

if includeStateTrajectory
    x = zeros(length(t),n);
else
    x = zeros(1,n);
end

x(1,:)= x0(:);
y = zeros(length(t),p);

kx = 1;
kxp = 1;



% discrete time simulation
for k = 1:length(t)
    % index only needed for 2nd output arg (state trajectory)
    % in case of nargout == 1 , x is calculated in place (continually overwriting the first line of x)
    if includeStateTrajectory
        kx = k;
        kxp = k+1;
    end
    % stack state-input vector
    % TODO!! Needs a better solution
    if(m > 0)
        xu = [x(kx,:)';u(k,:)'];
    else
        xu = x(kx,:)';
    end
    % compute output


    if(p ~= 0)
        y(k,:) = sim.sys.getOutput(xu);
    end

    % compute next state
    x(kxp,:) = sim.sys.getNextState(xu);

end

if (p == 0) 
    y = x;
end

tOut = t;
end

