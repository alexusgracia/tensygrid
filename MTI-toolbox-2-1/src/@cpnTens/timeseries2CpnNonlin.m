function [FU,FPhi,J]=timeseries2CpnNonlin(data,r,options)
%timeseries2CpnNonlin function for rank n parameter identification
%with normalized MTI models from iddata object with input-output data in
%time-domain
%22.02.2022
%leona.schnelle@haw-hamburg.de
%Input:                                     Output:
%   data    measurement iddata object           FU, Fphi CPN parameters                   
%   r       rank                                J   cost function result

ui=data.u;
x=data.y;
if isempty(ui)
    m=0;
else
    m = length(ui(1,:));            % number of inputs
end
if isempty(x)
    error('empty input-output data given')
else
    n =length(x(1,:)) ;             % model order 
end

ts=data.Ts;
tv = data.SamplingInstants;    % Time vector

W=nan(2*n+m,r);                % Initialize parameter matrix rank n 
    
x0 = x(1,:)';
%% Parameterid Optimization
if options.Normtype=="2-norm" 
    low=[zeros(r,n+m),options.lowerBound*ones(r,n)]';
    up=[pi/2*ones(r,n+m),100*ones(r,n)]';
elseif options.Normtype=="1-norm"
    low=[0*ones(r,n+m),options.lowerBound*ones(r,n)]';
    up=[1*ones(r,n+m),options.upperBound*ones(r,n)]';
end

% nonlinear least squares algorithm 
if options.Method=="lsqnonlin"
    lsq = 1;       % flag for lsq solver  
    loptions = optimoptions('lsqnonlin','Algorithm','levenberg-marquardt','MaxFunctionEvaluations',options.MaxFunctionNonlin,'Display',options.Display);
elseif options.Method=="fmincon"
    lsq=0;
    foptions = optimoptions('fmincon','Algorithm','interior-point','MaxFunctionEvaluations',options.MaxFunctionNonlin,'Display',options.Display);
else
    error("No valid optimization option.")
end
         
% genetic algorithm for global first guess
if options.Initialize=="ga"
    goptions = optimoptions('ga','MaxGenerations',20,'FitnessLimit',1,'Display','final');
    W0 = ga(@(W)cpnTens.OptiMTI(W,x,ui,tv,n,m,r,x0,ts,options.Focus,0,options.Normtype),length(W(:)),[],[],[],[],low(:)',up(:)',[],goptions);
elseif options.Initialize=="zero"
    W0=zeros(2*n+m,r); % initialize with zero parameter matrix
elseif options.Initialize=="random"
    W0=rand(2*n+m,r); % initialize with random parameter matrix
else
    error('No valid inizialization method given, choose "ga" for genetic algorithm, "zero" or "random"')
end
      
if lsq
    [W,J] = lsqnonlin(@(W)cpnTens.OptiMTI(W,x,ui,tv,n,m,r,x0,ts,options.Focus,0,options.Normtype),W0,low(:)',up(:)',loptions);
else
    [W,J] = fmincon(@(W)cpnTens.OptiMTI(W,x,ui,tv,n,m,r,x0,ts,options.Focus,0,options.Normtype),W0,[],[],[],[],low(:)',up(:)',[],foptions);
end
            
    
W=reshape(W,[2*n+m,r]);
FU=[W(1:n+m,:)];
FPhi=[W(n+m+1:n+m+n,:)];
disp(['identification finished:',' final cost: ',num2str(J(end))]);
end 