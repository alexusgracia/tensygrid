function dxdt=PLLnonlinear(t,x,Kp,Ki,w0,timeInput,voltageInput)

if size(voltageInput,1)~=3
    error('The number of rows voltageInput must be 3, try transposing it.')
end

u=[interp1(timeInput,voltageInput(1,:),t);interp1(timeInput,voltageInput(2,:),t);interp1(timeInput,voltageInput(3,:),t)];

% PLL tuning
k_p=Kp;
k_i=Ki;
w_0=w0;

%state equations
dxdt=[x(2)+k_p*1/3*((-2*u(1)+u(2)+u(3)).*sin(x(1))+sqrt(3)*(u(2)-u(3)).*cos(x(1)))+w_0;...
     k_i*1/3*((-2*u(1)+u(2)+u(3)).*sin(x(1))+sqrt(3)*(u(2)-u(3)).*cos(x(1)))];
end

