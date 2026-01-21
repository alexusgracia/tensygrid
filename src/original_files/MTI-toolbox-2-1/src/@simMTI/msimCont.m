function [y,tOut,x] = msimCont(sim, u, t, x0, method)
% MSIMCONT Simulate continuous-time time response of dynamic continuous-time explicit MTI model to arbitrary inputs. 

n = sim.sys.n;
p = sim.sys.p;
m = sim.sys.m;         % number of inputs
if ~isempty(sim.sys.G) && size(sim.sys.G.U,1)-n ~= m
    error('size mismatch')
end
% if(nargout>1 && p ==0)
%     warning('System does not have output tensor, the output trajecory will equal the state trajectory')
% end

[tOut, x] = sim.sys.getStateTrajectory(u, t, x0, method);

% Calculate Outputs at timestept tOut from the ode-solver,
% where also x is saved
if(p ~= 0)
    y = zeros(length(tOut),p);
    for k = 1:length(tOut)
        for j = 1:m
            u_at_tOut(j) = interp1(t, u(:,j), tOut(k), method);
        end
        xu = [x(k,:)'; u_at_tOut'];
        y(k,:) = sim.sys.getOutput(xu);
    end
end

if(nargout >= 1 && p == 0)
    y = x;
end

end

