function cost=OptiMTI(W,x,u,tv,n,m,r,x0,ts,os,lsq,ntype)
    %OptiMTI cost function to calculate the squared error sum between measured state and
    %simulated state with discrete time rank n normalized MTI model
    %22.02.2022
    %leona.schnelle@haw-hamburg.de
    %Inputs:                            Outputs:
    %   W   MTI model parameters            cost    squared error sum
    %   x   state measurement
    %   u   input measurement
    %   tv  time vector
    %   n   order of the system (number of states)
    %   m   number of inputs  
    %   r   CP rank of model
    %   x0  initial state
    %   lsq flag for lsq
    %   os =1 one step, =0 full simulation
    %   ntype  normalization p-norm ("1-norm","2-norm")

    % create normalized MTI model from parameters
    W=reshape(W,[2*n+m,r]);          % reshape parameter vector to matrix [F.U;F.phi]
 
    F.U=[W(1:n+m,:)];
    F.phi=[W(n+m+1:n+m+n,:)];         
    if(ntype ~= "1-norm")
        error('Currently unsupported normtype');
    end
    si = simMTI(mss(CPN1(F.U,F.phi),ts));

    if os=="simulation"                                          % full simulation 
       xsim=simulate(si,u,tv(end-1),x0);                         % simulation
    
        if lsq
            cost = xsim-x;
        else
            cost = norm(xsim-x,"fro");                                 % cost TODO clear 
        end
    elseif os=="prediction"                                      % one step simulation 
        %error('one step simulation not implemented yet')
            xsim=x;
            for k=1:length(tv)
                 xuk=[x(k,:) u(k,:)]';
                 %sys = mss(CPN1(U_temp,FPhi),ts);
                 xsim(k+1,:)=getNextState(mss(CPN1(F.U,F.phi),ts),xuk);
            end

%         % plot([xsim,x]), drawnow, disp(W)
         if lsq
             cost = xsim(1:end-1,:)-x;
         else
             cost = norm(xsim(1:end-1,:)-x,"fro");
         end
    end
end

