function [FU,FPhi,cost]=timeseries2CpnAls_nonneg(data,r,options)
%TIMESERIES2CPNALS function for rank n parameter identification
%for normalized MTI factor matrices from iddata object with input-output data in
%time-domain using alternating linear scheme algorithm
%
%   Input:                                  Output:
%                                                                                   
%   data        measurement iddata object       FU      structure parameter matrix
%   r   	    rank                            FPhi    CPN MTI Model                                                            
%   options:    AlsIterations                   J       cost function result
%               AlsTolerance
%               Display ("on", "off")
%               MaxInitialization    
%               Initialize ("random", "zero")
%               Focus ("simulation", "prediction")
%    
%   28.08.2023
%   leona.schnelle@haw-hamburg.de
%
%   See also mlgreyest
maxiter=options.AlsIterations-1;
tol=options.AlsTolerance;
disp_flag=false;

iter=0;
cost_k=inf;
cost=[nan,inf];
init=0;

us=data.u;
xs=data.y;
if isempty(us)
    m=0;
else
    m = length(us(1,:));       % number of inputs
end
if isempty(xs)
    error('empty input-output data given')
else
    n =length(xs(1,:)) ;       % model order 
end

ts=data.Ts;                    % sampling time
tv = data.SamplingInstants;    % Time vector


% check for data length and parameter number
if length(tv)-1<2*r
    error('Underdeterminded equation system. Reduce rank or provide more data.')
end
% Initialize parameters
if isstring(options.Initialize)
    if options.Initialize=="random"
    random_flag=true;
    FU=rand(n+m,r);
    FPhi=rand(n,r); 
    init=init+1;
elseif options.Initialize=="zero"
    random_flag=false;
    FU=zeros(n+m,r);
    FPhi=zeros(n,r); 
    init=init+1;
else
    error("No valid initialization method given for method ALS. Must be 'zero' or 'random' or 'diagonal'." )
    end
else 
    FU=options.Initialize.U;
    FPhi=options.Initialize.Phi;  
end
        if options.Focus=="simulation"
            xsim=msim(mss(CPN1(FU,FPhi),ts),us(tv(1:end-1),:),tv(1:end-1),xs(1,:));
        elseif options.Focus=="prediction"      % TODO: use getNextState for matrix input xu and remove loop
            xsim=nan(length(tv),n);
            xsim(1,:)=xs(1,:);
            for k=1:length(tv)-1
                 xuk=[xs(k,:) us(k,:)]';
                 xsim(k+1,:)=getNextState(mss(CPN1(FU,FPhi),ts),xuk);
            end
        else
            error("No valid focus option for method als.")
        end
        cost_k=norm(xsim(2:end,:)-xs(tv(2:end),:),"fro");
       
%% alternating linear scheme
while (cost_k>tol || isnan(cost_k)) && maxiter>=iter 
    iter=iter+1;
    for ind=1:n+m                               % index of free parameter

        indother=true(n+m,1);                   % index of fixed parameters 
        indother(ind,1)=false;

        xu=[xs(tv,:) us(tv,:)];
        % Calculation of polynomial for fixed parameters for every time step

        f_other=nan(length(tv)-1,r*n);          % initialize matrix
        for k=1:length(tv)-1
            f_other(k,:)=reshape(FPhi.*prod((1-abs(FU(indother,:)))+FU(indother,:).*xu(k,indother)'),1,n*r);
        end
        % Linear solving    

        A=reshape([f_other xu(1:end-1,ind).*f_other],n*length(f_other),2*r);
        b=xu(2:end,1:n);
        G=reshape(lsqnonneg(A,b(:)),[r 2])';
        %G=reshape(A\b(:),[r 2])';               % [g1r1 g1r2 g1r2 ... ;g2r1 g2r2 g2r3 ...]
        norm1=vecnorm(G,1); 
        U_temp=FU;
        b=G(1,:)<0;                             % TODO: use general function cp2cpn()  
        U_temp(ind,:)=G(2,:)./norm1.*(~b-b);    % normalize paramater in 1-norm
        FPhi_temp=norm1.*FPhi.*(~b-b);
        if options.Focus=="simulation"
            xsim=msim(mss(CPN1(U_temp,FPhi_temp),ts),us(tv(1:end-1),:),tv(1:end-1),xs(1,:));
        elseif options.Focus=="prediction"      % TODO: use getNextState for matrix input xu and remove loop
            xsim=nan(length(tv),n);
            xsim(1,:)=xs(1,:);
            for k=1:length(tv)-1
                 xuk=[xs(k,:) us(k,:)]';
                 xsim(k+1,:)=getNextState(mss(CPN1(U_temp,FPhi_temp),ts),xuk);
            end
        else
            error("No valid focus option for method als.")
        end
        cost_k=norm(xsim(2:end,:)-xs(tv(2:end),:),"fro");
        % Update FU and FPhi, if cost better than the best fit
        if cost(end)>=cost_k && ~isnan(cost_k)
        FU=U_temp;
        FU(FU<0)=0;
        FPhi=FPhi_temp;
        cost=[cost(end);cost_k];
        end

        if options.Display=="on"
            disp_flag=true;
            disp(['iterations: ',num2str(iter),' Index: ',num2str(ind),' Cost: ',num2str(cost_k)]);
        elseif options.Display~="off"
            error("No valid display option value given. Display option value for method als are 'on' or 'off'.")

        end
    end

    

% Solve for phi parameter
A_phi=nan(length(tv)-1,r);                  % initialize A-matrix for linear problem Ax=b
for k=1:length(tv)-1
    A_phi(k,:)=prod((1-abs(FU))+FU.*xu(k,:)');
end
b=xu(2:end,1:n);                            % build b-matrix fro linear problem Ax=b
FPhi_temp=(A_phi\b)';                       % solve linear equation system
        if options.Focus=="simulation"
            xsim=msim(mss(CPN1(FU,FPhi_temp),ts),us(tv(1:end-1),:),tv(1:end-1),xs(1,:));
        elseif options.Focus=="prediction"  % TODO: use getNextState for matrix input xu and remove loop
            xsim=nan(length(tv),n);
            xsim(1,:)=xs(1,:);
            for k=1:length(tv)-1
                 xuk=[xs(k,:) us(k,:)]';
                 xsim(k+1,:)=getNextState(mss(CPN1(FU,FPhi_temp),ts),xuk);
            end
        end

cost_k=norm(xsim(2:end,:)-xs(tv(2:end),:),"fro");

% Update FPhi, if cost better than the best fit
if cost(end)>=cost_k && ~isnan(cost_k)
FPhi=FPhi_temp;                  
cost=[cost(end);cost_k];
end 

if disp_flag
    disp(['iterations: ',num2str(iter),' Index: ','Phi',' Cost: ',num2str(cost_k)]);
end

% breaking loop criteria and new initialization in case of bad results
    if abs(cost(2)-cost(1))<=0.0001 && iter<maxiter
        if cost(2)<=options.ToleranceInitialization || init>=options.MaxInitialization
            fprintf('Change from one iteration to the next is lower than the limit. \n')
            drawnow
            break;
        else
        % New initialization parameters
        if random_flag
            FU=rand(n+m,r);
            FPhi=rand(n,r);
        else
            fprintf('Change from one iteration to the next is lower than the limit. \n')
            drawnow
            break;
        end
        init=init+1;
        iter=0;
        cost(2)=inf;
        if disp_flag
            fprintf('New initialization. \n')
        end
        end
    end

    if isinf(cost(2)) && iter > 3 && iter<maxiter
        % New initialization parameters
        if random_flag
            FU=rand(n+m,r);
            FPhi=rand(n,r);
        else
            FU=zeros(n+m,r);
            FPhi=zeros(n,r);  
        end
        init=init+1;
        iter=0;
        cost(2)=inf;
        fprintf('New initialization. \n')
    end
end
cost=cost(2);
disp(['Identification finished: ' 'number of iterations: ',num2str(iter),', final cost: ',num2str(cost(end))]);
end