function [Fphi, FaInter, Fc, cost] = otvl2AlsA(otvlTensObj, data, options)
%OTVL2ALSA Function performing parameter identification with an ALS
%algorithm for a given otvlTens object.
% inputs: 
%   otvlTensObj: otvlTens object containing structure of multilinear state
%   equations
%   data: measurement iddata object
%   options
%
% outputs:
%   Fphi: estimated values for parameter FPhi
%   FaInter: estimated values for internal parameter of Fa
%   Fc: estimated valued for parameter Fc
%
% See also otvl otvlTens mlgreyest

% Marah Engels - 31/05/2024
    
    if isempty(otvlTensObj.structure)
        error('No structure given for otvlmssObj.')
    end

    %read information from otvlmss object
    N_EQUATIONS = otvlTensObj.NEquations; %number of state equations given in TVLs
    nRowsTotal = otvlTensObj.NRowsTotal; % sum of rows in all TVLs of the system
    dcTvl = otvlTensObj.DontCareArray; % merged Dont-Care-arrays of all TVLs
    bTvl = otvlTensObj.BooleanArray; % merged Boolean-arrays of all TVLs
    index = otvlTensObj.IndexVector; % index array to allocate TVL rows to equations
    FPhiStruct = otvlTensObj.FPhi; % FPhi for model structure scaling
    %FlambdaStruct = otvlmssObj.Falpha; %Falpha for model structure
    FaInterStruct = otvlTensObj.FaInter; % internal Fa for model structure adjustment
    FcStruct = otvlTensObj.Fc; % Fc for model structure offset
    %[lambdaRow, lambdaCol] = find(FlambdaStruct); %get indices of 0 in
    %model structure
    %lambdaIndices = [lambdaRow lambdaCol]';
    

    %read options for ALS algorithm
    MAX_ITERATIONS = options.AlsIterations-1;
    TOLERANCE = options.AlsTolerance;
    MAX_INITIALIZATIONS = options.MaxInitialization;
    

    %initialize parameters
    iteration = 0;
    cost=[nan; inf];
    costPhi = [nan; inf];
    costA = [nan; inf];
    cost_k = inf;
    initialization = 0;

    %create test-arrays to store cost of identification
    costArr = nan(MAX_ITERATIONS,1);
    costPhiArr = nan(MAX_ITERATIONS,1);
    costaArr = nan(MAX_ITERATIONS,1);


    xSys = data.y; %state trajectories
    timeVector = data.SamplingInstants;
    timeSteps = length(timeVector);

    if isempty(xSys)
        error('empty input-output data given')
    else
        xInit = xSys(1,:);
    end

    if ~(size(xSys,2) == N_EQUATIONS)
        error(['size of state trajectory is not consistent with number of ' ...
            'equations given in otvlmss'])
    end
    
    % initialize parameters
    if options.Initialize == "random"
        flaginit = 0;
        Fphi = rand(size(FPhiStruct));
        %Flambda = sprand(FlambdaStruct);
        FaInter = rand(size(FaInterStruct));
        Fc = rand(size(FcStruct));
        initialization = initialization+1;
    elseif options.Initialize == "none"
        flaginit = 1;
        Fphi = ones(size(FPhiStruct));
        %Flambda = spones(FlambdaStruct);
        FaInter = zeros(size(FaInterStruct));
        Fc = zeros(size(FcStruct));
        initialization = initialization+1;
    elseif options.Initialize == "user"
        flaginit = 1;
        Fphi = FPhiStruct;
        %Flambda = FlambdaStruct;
        FaInter = FaInterStruct;
        Fc = FcStruct;
        initialization = initialization+1;
    else
        error("No valid initialzation method given for method ALS. " + ...
            "Must be 'zero' or 'random'.\n" )
    end
    
    %update Parameters of otvlss object with initialized parameters
    otvlTensObj.FPhi = Fphi;
    %otvlmssObj.Falpha = Flambda;
    otvlTensObj.FaInter = FaInter;
    otvlTensObj.Fc = Fc;

    %Initialize best found parameters and cost
    FphiBest = Fphi;
    FaInterBest = FaInter;
    FcBest = Fc;
    bestSim = otvlsimA(otvlTensObj,[],timeVector, xInit);
    bestCost = norm(bestSim(:) - xSys(:), 'fro');

    while (cost_k>TOLERANCE || isnan(cost_k)) && MAX_ITERATIONS > iteration %perform ALS algorithm
         
         iteration = iteration + 1;
         
         for indPhi = 1:nRowsTotal
            %estimate Fphi

            %create constant equation parts. const1 = part of equation 
            % adjacent to Fphi(ind), const2 = rest of equation
            
            %const1
            const1Phi = prod(dcTvl(indPhi,:) ...
                + (~dcTvl(indPhi,:))...
                .*(bTvl(indPhi,:).*(xSys(1:end-1,:) + FaInter(1,:))...
                + (~bTvl(indPhi,:)).*(FaInter(2,:) - xSys(1:end-1,:)))...
                ,2); %product of TVL-elements in row given by indPhi 

            %const2
            nonrelInd = find(~(index == index(indPhi,1)))';
            const2Phi = reshape(Fphi(setdiff(1:end,[indPhi nonrelInd]),:).*prod( ...
                dcTvl(setdiff(1:end,[indPhi nonrelInd]),:) ...
                + (~dcTvl(setdiff(1:end,[indPhi nonrelInd]),:)) ...
                .*(bTvl(setdiff(1:end,[indPhi nonrelInd]),:) ...
                .*permute(xSys(1:end-1,:) + FaInter(1,:), [3,2,1]) ...
                + (~bTvl(setdiff(1:end,[indPhi nonrelInd]),:)) ...
                .*(permute(FaInter(2,:) - xSys(1:end-1,:), [3,2,1]))) ...  
                ,2),[],timeSteps-1,1); %product of TVL-elements in all other rows

            const2Phi = sum(const2Phi,1)';
        
            % for deterministic systems, parts of the A matrix equal 0. Might
            % have to be adjusted for non-deterministic systems
            A = reshape(const1Phi,[],1);
            b = reshape(xSys(2:end,index(indPhi,:)) - const2Phi - Fc(index(indPhi,:),:),[],1);
    
            paramPhi = A\b;
            
            %update Fphi
            Fphi(indPhi,:) = paramPhi;
            otvlTensObj.FPhi = Fphi;
    
         end

         %perform simulation with updated parameters and get cost
        if options.Focus=="simulation"
            xOtvlSim = otvlsimA(otvlTensObj,[],timeVector, xInit);
        elseif options.Focus=="prediction"
            xPre = nan(timeSteps,nRowsTotal); 
            xOtvlSim = nan(timeSteps, N_EQUATIONS);
            xOtvlSim(1,:) = xInit;
            for k = 1:timeSteps-1
                xPre(k+1,:) = (otvlTensObj.FPhi.*prod(dcTvl + (~dcTvl) ...
                    .*(bTvl.*(xSys(k,:) + otvlTensObj.FaInter(1,:)) ...
                    +(~bTvl).*(otvlTensObj.FaInter(2,:) - xSys(k,:))),2))';
                xOtvlSim(k+1,:) = accumarray(index, xPre(k+1,:)')' + otvlTensObj.Fc';
            end
        else
            error("No valid focus option for method als.")
        end
         cost_k = norm(xOtvlSim - xSys, 'fro');
         costPhi = [costPhi(end);cost_k];

         costPhiArr(iteration,1) = costPhi(2);

         for ind = 1:N_EQUATIONS
            %estimate Fa
            %case 1: x+a1

            %const1
            constVals1a1 = reshape(Fphi(bTvl(:,ind) == 1,:) ...
                 .*prod(dcTvl(bTvl(:,ind) == 1,setdiff(1:end,ind)) ...
                + (~dcTvl(bTvl(:,ind) == 1,setdiff(1:end,ind))) ...
                .*(bTvl(bTvl(:,ind) == 1,setdiff(1:end,ind)) ...
                .*permute(xSys(1:end-1,setdiff(1:end,ind)) ...
                + FaInter(1,setdiff(1:end,ind)), [3,2,1]) ...
                + (~bTvl(bTvl(:,ind) == 1,setdiff(1:end,ind))) ...
                .*(FaInter(2,setdiff(1:end,ind)) ...
                - permute(xSys(1:end-1,setdiff(1:end,ind)), [3,2,1]))) ...  
                ,2),[],timeSteps-1,1)'; %products of all rows, where TVL(:,ind)=1 excluding column ind
            %create indexes for accumarray const1
            subs1Timea1 = repmat((1:timeSteps-1)',sum(bTvl(:,ind) == 1),1); %timesteps repeated as many times, as rows where TVL(:,ind) = 1 exist
            subs1Indexa1 = repelem(index(bTvl(:,ind) == 1,:), timeSteps-1); %all indexes where TVL(:,ind) = 1 repeated 

             const1a1 = ...
                accumarray([subs1Timea1(:), subs1Indexa1(:)], ...
                constVals1a1(:), [], @(x) sum(x));%sum based on Index vector
             const1a1 = const1a1(:,unique(subs1Indexa1));

            %const2
            bTvlRelevant = bTvl(ismember(index,unique(index(bTvl(:,ind) == 1,:))),:);
            dcTvlRelevant = dcTvl(ismember(index,unique(index(bTvl(:,ind) == 1,:))),:);
            FphiRelevant = Fphi(ismember(index,unique(index(bTvl(:,ind) == 1,:))),:);
            FcRelevant = Fc(unique(index(bTvl(:,ind) == 1,:)),:);

            constVals2a1 = reshape(FphiRelevant(~(bTvlRelevant(:,ind) == 1),:) ...
                 .*prod(dcTvlRelevant(~(bTvlRelevant(:,ind) == 1),:) ...
                + (~dcTvlRelevant(~(bTvlRelevant(:,ind) == 1),:)) ...
                .*(bTvlRelevant(~(bTvlRelevant(:,ind) == 1),:) ...
                .*permute(xSys(1:end-1,:) ...
                + FaInter(1,:), [3,2,1]) ...
                + (~bTvlRelevant(~(bTvlRelevant(:,ind) == 1),:)) ...
                .*(FaInter(2,:) - permute(xSys(1:end-1,:), [3,2,1]))) ...  
                ,2),[],timeSteps-1,1)'; %products of all rows, where TVL(:,ind)=0
            %create indexes for accumarray const2
            subs2Timea1 = repmat((1:timeSteps-1)',sum(~(bTvlRelevant(:,ind) == 1)),1); %timesteps repeated as many times, as rows where TVL(:,ind) = 0 exist
            subs2Indexa1 = repelem(index(~(bTvlRelevant(:,ind) == 1),:), timeSteps-1); %all indexes where TVL(:,ind) = 0 repeated 

            const2a1 = ...
                accumarray([subs2Timea1(:), subs2Indexa1(:)], ...
                constVals2a1(:), [], @(x) sum(x)); %sum based on Index vector
            const2a1 = const2a1(:,unique(subs2Indexa1));

            %estimate Fa1
            A = reshape(const1a1,[],1);
            b = reshape(xSys(2:end,unique(index(bTvl(:,ind) == 1,:))) - const2a1 - FcRelevant' - ...
                repmat(xSys(1:end-1,ind),1,length(unique(index(bTvl(:,ind) == 1,:)))).*const1a1,[],1);
            parama1 = lsqnonneg(A,b);

            % update Fa1
            FaInter(1, ind) = parama1;
            otvlTensObj.FaInter = FaInter;

            %case 2: a2-x
            %const1
            constVals1a2 = reshape(Fphi(or(dcTvl(:,ind),bTvl(:,ind)) == 0,:) ...
                 .*prod(dcTvl(or(dcTvl(:,ind),bTvl(:,ind)) == 0,setdiff(1:end,ind)) ...
                + (~dcTvl(or(dcTvl(:,ind),bTvl(:,ind)) == 0,setdiff(1:end,ind))) ...
                .*(bTvl(or(dcTvl(:,ind),bTvl(:,ind)) == 0,setdiff(1:end,ind)) ...
                .*permute(xSys(1:end-1,setdiff(1:end,ind)) ...
                + FaInter(1,setdiff(1:end,ind)), [3,2,1]) ...
                + (~bTvl(or(dcTvl(:,ind),bTvl(:,ind)) == 0,setdiff(1:end,ind))) ...
                .*(FaInter(2,setdiff(1:end,ind)) ...
                - permute(xSys(1:end-1,setdiff(1:end,ind)), [3,2,1]))) ...  
                ,2),[],timeSteps-1,1)'; %products of all rows, where TVL(:,ind)=0 excluding column ind
            %create indexes for accumarray const1
            subs1Timea2 = repmat((1:timeSteps-1)',sum(or(dcTvl(:,ind),bTvl(:,ind)) == 0),1); %timesteps repeated as many times, as rows where TVL(:,ind) = 0 exist
            subs1Indexa2 = repelem(index(or(dcTvl(:,ind),bTvl(:,ind)) == 0,:), timeSteps-1); %all indexes where TVL(:,ind) = 0 repeated 
            const1a2 = zeros(timeSteps-1, size(FcRelevant,1));

            const1a2_pre = ...
                accumarray([subs1Timea2(:), subs1Indexa2(:)], ...
                constVals1a2(:), [], @(x) sum(x));%sum based on Index vector
            if size(const1a2_pre,2) > size(unique(subs1Indexa2))
                const1a2_pre(:,all(const1a2_pre == 0)) = [];
                const1a2 = const1a2_pre;
            else
                const1a2(:,unique(subs1Indexa2)) = const1a2_pre;
            end

            %const2
            bTvlRelevant = bTvl(ismember(index,unique(index(or(dcTvl(:,ind),bTvl(:,ind)) == 0,:))),:);
            dcTvlRelevant = dcTvl(ismember(index,unique(index(or(dcTvl(:,ind),bTvl(:,ind)) == 0,:))),:);
            FphiRelevant = Fphi(ismember(index,unique(index(or(dcTvl(:,ind),bTvl(:,ind)) == 0,:))),:);
            FcRelevant = Fc(unique(index(or(dcTvl(:,ind),bTvl(:,ind)) == 0,:)),:);

            constVals2a2 = reshape(FphiRelevant(or(dcTvlRelevant(:,ind),bTvlRelevant(:,ind)) == 1,:) ...
                 .*prod(dcTvlRelevant(or(dcTvlRelevant(:,ind),bTvlRelevant(:,ind)) == 1,:) ...
                + (~dcTvlRelevant(or(dcTvlRelevant(:,ind),bTvlRelevant(:,ind)) == 1,:)) ...
                .*(bTvlRelevant(or(dcTvlRelevant(:,ind),bTvlRelevant(:,ind)) == 1,:) ...
                .*permute(xSys(1:end-1,:) ...
                + FaInter(1,:), [3,2,1]) ...
                + (~bTvlRelevant(or(dcTvlRelevant(:,ind),bTvlRelevant(:,ind)) == 1,:)) ...
                .*(FaInter(2,:) ...
                - permute(xSys(1:end-1,:), [3,2,1]))) ...  
                ,2),[],timeSteps-1,1)';%products of all rows, where TVL(:,ind)=1
             %create indexes for accumarray
            subs2Timea2 = repmat((1:timeSteps-1)',sum(or(dcTvlRelevant(:,ind),bTvlRelevant(:,ind)) == 1),1); %timesteps repeated as many times, as rows where TVL(:,ind) = 1 exist
            subs2Indexa2 = repelem(index(or(dcTvlRelevant(:,ind),bTvlRelevant(:,ind)) == 1,:), timeSteps-1); %all indexes where TVL(:,ind) = 1 repeated 
            const2a2 = zeros(timeSteps-1, size(FcRelevant,1));

            const2a2_pre = ...
                accumarray([subs2Timea2(:), subs2Indexa2(:)], ...
                constVals2a2(:), [], @(x) sum(x));%sum based on Index vector;
            if size(const2a2_pre,2) > size(unique(subs2Indexa2))
                const2a2_pre(:,all(const2a2_pre == 0)) = [];
                const2a2 = const2a2_pre;
            else
                const2a2(:,unique(subs2Indexa2)) = const2a2_pre;
            end

            %estimate Fa2
            A = reshape(const1a2,[],1);
            b = reshape(xSys(2:end,unique(index(or(dcTvl(:,ind),bTvl(:,ind)) == 0,:))) - const2a2 - FcRelevant' + ...
                repmat(xSys(1:end-1,ind),1,length(unique(index(or(dcTvl(:,ind),bTvl(:,ind)) == 0,:)))).*const1a2,[],1);
            parama2 = lsqnonneg(A,b);

            % update Fa2
            FaInter(2, ind) = parama2;
            otvlTensObj.FaInter = FaInter;

         end

         %perform simulation with updated parameters and get cost
        if options.Focus=="simulation"
            xOtvlSim = otvlsimA(otvlTensObj,[],timeVector, xInit);
        elseif options.Focus=="prediction"
            xPre = nan(timeSteps,nRowsTotal); 
            xOtvlSim = nan(timeSteps, N_EQUATIONS);
            xOtvlSim(1,:) = xInit;
            for k = 1:timeSteps-1
                xPre(k+1,:) = (otvlTensObj.FPhi.*prod(dcTvl + (~dcTvl) ...
                    .*(bTvl.*(xSys(k,:) + otvlTensObj.FaInter(1,:)) ...
                    +(~bTvl).*(otvlTensObj.FaInter(2,:) - xSys(k,:))),2))';
                xOtvlSim(k+1,:) = accumarray(index, xPre(k+1,:)')' + otvlTensObj.Fc';
            end
        else
            error("No valid focus option for method als.")
        end
         cost_k = norm(xOtvlSim - xSys, 'fro');
         costA = [costA(end);cost_k];

         costaArr(iteration,1) = costA(2);

         %estimate Fc
        constAllVals = reshape(Fphi(:,:) ...
             .*prod(dcTvl(:,:) ...
            + (~dcTvl(:,:)) ...
            .*(bTvl(:,:) ...
            .*permute(xSys(1:end-1,:) ...
            + FaInter(1,:), [3,2,1]) ...
            + (~bTvl(:,:)) ...
            .*(FaInter(2,:) ...
            - permute(xSys(1:end-1,:), [3,2,1]))) ...  
            ,2),[],timeSteps-1,1)'; %sum over products of all TVL rows
        subs1all = repmat((1:timeSteps-1)',nRowsTotal,1); %timesteps repeated as many times, as there are TVL-rows
        subs2all = repelem(index, timeSteps-1); %all Indexes repeated
        constAll= accumarray([subs1all(:), subs2all(:)], ...
            constAllVals(:), [], @(x) sum(x)); %sum based on Index vector

        A = ones(timeSteps-1,1);
        b = xSys(2:end,:) - constAll;

        %estimate and update Fc
        Fc = (A\b)';
        otvlTensObj.Fc = Fc;

        %perform simulation with updated parameters and get cost
        if options.Focus=="simulation"
            xOtvlSim = otvlsimA(otvlTensObj,[],timeVector, xInit);
        elseif options.Focus=="prediction"
            xPre = nan(timeSteps,nRowsTotal); 
            xOtvlSim = nan(timeSteps, N_EQUATIONS);
            xOtvlSim(1,:) = xInit;
            for k = 1:timeSteps-1
                xPre(k+1,:) = (otvlTensObj.FPhi.*prod(dcTvl + (~dcTvl) ...
                    .*(bTvl.*(xSys(k,:) + otvlTensObj.FaInter(1,:)) ...
                    +(~bTvl).*(otvlTensObj.FaInter(2,:) - xSys(k,:))),2))';
                xOtvlSim(k+1,:) = accumarray(index, xPre(k+1,:)')' + otvlTensObj.Fc';
            end
        else
            error("No valid focus option for method als.")
        end

         cost_k = norm(xOtvlSim - xSys, 'fro');
         cost = [cost(end);cost_k];

         costArr(iteration,1) = cost(2);
        % if cost(2) - cost(1) > 0.0001
        %     warning('Cost Tvl Als increased\n')
        %     fprintf('Cost increased with %f at otvlals at iteration %i, initialization %i. \n',cost(2) - cost(1), iteration, initialization)
        % end

         %break loop if cost does converge at high value
            if ((abs(cost(2)-cost(1))<=0.0001)||(isinf(cost(2)) && iteration>3)) && iteration<MAX_ITERATIONS
                if cost(2)<0.1 || initialization>=MAX_INITIALIZATIONS || flaginit
                %fprintf('Change from one iteration to the next is lower than the limit. \n')
                break
                else
                %Compare found parameters with best found parameters
                if cost(2) < bestCost
                    %update best parameters if cost has decreased
                    FphiBest = Fphi;
                    FaInterBest = FaInter;
                    FcBest = Fc;
                    bestCost = cost(2);
                    if options.Display == "on"
                        fprintf('Best parameters updated at iteration %i, initialization %i. \n', iteration, initialization');
                    end
                end

                % New initialization parameters
                Fphi = rand(size(FPhiStruct));
                FaInter = rand(size(FaInterStruct));
                Fc = rand(size(FcStruct));
                %Flambda = sprand(FlambdaStruct);
                otvlTensObj.FPhi = Fphi;
                otvlTensObj.FaInter = FaInter;
                otvlTensObj.Fc = Fc;
                initialization = initialization+1;
                if options.Display == "on"
                    fprintf('New initialization otvl at iteration %i, initialization %i. \n', iteration, initialization);
                end
                    iteration = 0;
                costArr = nan(MAX_ITERATIONS,1);
                costPhiArr = nan(MAX_ITERATIONS,1);
                costaArr = nan(MAX_ITERATIONS,1);
                cost(2) = inf;
                costA(2) = inf;
                costPhi(2) = inf;
                end
            end
         
    end

    %Get best identified parameters in case cost is lower
    if bestCost < cost(2)
        otvlTensObj.FPhi = FphiBest;
        otvlTensObj.FaInter = FaInterBest;
        otvlTensObj.Fc = FcBest;
        cost = bestCost;
    else
        cost = cost(2);
    end
    
    otvlTensObj.Fa(2,:) = 1-otvlTensObj.FaInter(2,:); %update external Fa

    
    disp(['Identification finished: ' 'number of iterations: ', ...
        num2str(iteration),', number of initializations: ', ...
        num2str(initialization),', final cost: ',num2str(cost(end))]);
    if options.Display == "on"
        %create test plots
        figure('Name','Cost ALS TVL')
        plot(costArr, 'DisplayName','Overall cost')
        title('Convergence behavior of ALS algorithm')
        ylabel('cost')
        xlabel('iteration')
        ylim auto
        legend show
    end

end