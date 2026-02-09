function [sys,J] = mlgreyest(data,r,structuralRestrictions,options)
%MLGREYEST function for rank n parameter identification
%with normalized MTI models or structurally restricted MTI models
%from iddata object with input-output data in time-domain
%
%   Input:                                      Output:
%   data    measurement iddata object           sys     mss object                                                            
%   r   	rank                                J       cost function result
%   structuralRestrictions otvlTens object [optional]
%   options
%
%   [sys,J] = mlgreyest(data,r) identifies a multilinear Model in 1-norm
%   CPN representation with default options
%
%   [sys,J] = mlgreyest(data,1,structuralRestrictions) identifies a
%   multilinear model with a restricted multilinear structure, which is
%   given as an otvlTens object
%
%   [sys,J] = mlgreyest(data,r, options) identifies a multilinear model
%   with custom options
%   The supported options are:
%  
%   Initialize  Handling of initial parameters (F0) during estimation.
%               Specify as one of the following:
%               'random': F0 is initialized random. Default.
%               'estimate': Initial parameters are estimated with global
%                optimization algorithm. Not supported for alternating
%                linear scheme algorithm
%
%   Display       View estimation progress ('on'/'off'). Default: 'off'.
%                
%
%   Normtype      spezifies the normalization typ of the mti cpn object. 
%                 Default: "1-norm"
%        
%   Method        The method for the optimization in set here. 
%                 'als': for altenating linear scheme optimization
%                 'fmincon': for general nonlinear optimization using the interior point algorithm
%                 'lsqnonlin': for nonlinear least squares optimization.
%                 Default: 'als'
%   Nonnegative   set to 'true' to perform nonnegative parameter estimation
%                 with ALS algorithm. 
%                 Default: 'false'
%
%   AlsTolerance  criteria to terminate the alternating linear scheme
%                 optimization if the value in AlsTolerance is reached.
%                 Default: 1
%
%   AlsIterations criteria to terminate the alternating linear scheme
%                 optimization if themaximum number of iterations is reached.
%                 Default: 200
%
%   MaxFunctionNonlin maximal number of function evaluation in nonlinear
%                  optimization
%                  Default: 1000000
% 
%
%   [sys,J] = mlgreyest(data,1,structuralRestrictions, options) identifies a
%   multilinear model with a multilinear structure given by an otvlTens-object. 
%   The supported options are:
%   
%   Initialize  Handling of initial parameters (FPhi, Fa, Fc) during estimation.
%               Specify as one of the following:
%               'random': parameters are initialized random. Default.
%               'none': parameters are initialized by ones and zeros
%               preserving the structure contained in the otvlTens object
%               'user': parameters are read from the otvlTens object
%
%
%   Display       View estimation progress ('on'/'off'). Default: 'off'.
%
%   AlsTolerance  criteria to terminate the alternating linear scheme
%                 optimization if the value in AlsTolerance is reached.
%                 Default: 1
%
%   AlsIterations criteria to terminate the alternating linear scheme
%                 optimization if themaximum number of iterations is reached.
%                 Default: 200
%
%   
% For detailed documentation see <a href="matlab:open((which('mlgreyestDoc.html')))">here</a>
% 
%
% See also timeseries2CpnAls timeseries2CpnNonlin otvlTens

%  Marah Engels, Leona Schnelle, Enrico Uhlenberg - 12/06/2024

arguments
         data 
         r
         structuralRestrictions (1,:) otvlTens=otvlTens().empty
         options.lowerBound  double = -10
         options.upperBound  double = 10
         options.Normtype (1,1) string = "1-norm"
         options.Focus (1,1) string = "simulation"      % "prediction" for one step optimization
         options.Method (1,1) string = "als"            %"fmincon" %"lsqnonlin" %"als"
         options.Initialize ="random"                   %"estimate" "zero" %TODO: possibilty of different default parameter for als and fmincon 
         options.Display (1,1) string ="off"            %"on" for fmincon: "iter"
         options.AlsTolerance (1,1) double=0.1
         options.AlsIterations (1,1) double=200
         options.MaxInitialization (1,1) int16=5 
         options.ToleranceInitialization (1,1) double=0.5 %Tolerance until new initialization is done
         options.Nonnegative logical=false %non negative least squares
         options.MaxFunctionNonlin  (1,1) double=100000
         options.Identification = "states"
end

if ~isempty(structuralRestrictions) % if user requests otvlTens ALS 
    warning('Rank limitation is currently only implemented through the strucutral restrictions');
    if length(data)>1
        error('tvl Identification with output Tensor not yet implemented')
    end
    boolmssobj = structuralRestrictions;
    switch lower(options.Identification)
        case 'states'
            if isa(boolmssobj,'otvlTens')
                [boolmssobj.FPhi, boolmssobj.Fa, boolmssobj.Fc, J] = ...
                otvl2AlsA(boolmssobj,data, options);
     
                sys = boolmssobj;
            else
                error(['Wrong data type in first position for chosen...' ...
                    ' identification type. Argument has to be of type ' ...
                    'otvlmss for identification of type states.'])
            end

        case 'outputs'
            if isempty(data.u)
                [boolmssobj.Fp, boolmssobj.Fc,  boolmssobj.Fa, ...
                    boolmssobj.Fphi, J] = ...
                bvl2Als(boolmssobj, stateEstimate, data, options);
            else
                [boolmssobj.Fp, boolmssobj.Fc,  boolmssobj.Fa, boolmssobj.Fd, ...
                    boolmssobj.Fphi, J] = ...
                bvl2AlsInput(boolmssobj, stateEstimate, data, options);
            end

            sys = boolmssobj;
            
    end
else % User requested CPN1 ALS from data - no structural restritction via tvl

     if options.Method=="als"
         if options.Nonnegative==true
         [FU,FPhi,J] =cpnTens.timeseries2CpnAls_nonneg(data,r,options);
         else
         [FU,FPhi,J] =cpnTens.timeseries2CpnAls(data,r,options);
         end
     else
         [FU,FPhi,J] =cpnTens.timeseries2CpnNonlin(data,r,options);
     end
    

    sys=mss(CPN1(FU,FPhi),data.Ts);
end
end

