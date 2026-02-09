% Extracting operating point from the Simulink model
% KauChr, 06.06.2025

tSimIdx=find(out.tout>=timeLinearization,1);
%%
xEP=zeros(MTImodel.n,1);
xEP_transfer=zeros(MTImodel.n,1);
% strNames=string(out.xout.getElementNames);

matchingElements.SL.str=out.xout.getElementNames;
matchingElements.SL.strIdx=[1:out.xout.numElements]; % get the index for the string

matchingElements.MTI.str=MTImodel.stateName;
matchingElements.MTI.strIdx=[1:MTImodel.n];

% 1. delete <missing>/"" entries
matchingElements.SL.strIdx=matchingElements.SL.strIdx(~ismissing(matchingElements.SL.str));
matchingElements.SL.str=matchingElements.SL.str(~ismissing(matchingElements.SL.str));


% 2. find matching elements

matchingElements.m1.names=matchingElements.SL.str(contains(matchingElements.SL.str,matchingElements.MTI.str));
matchingElements.m1.idx=matchingElements.SL.strIdx(contains(matchingElements.SL.str,matchingElements.MTI.str));

  % Create a logical array where each element indicates if it's found
    isContained = cellfun(@(x) any(contains(matchingElements.SL.str, x)), matchingElements.MTI.str);
    
    % Get the indices of matchingElements.MTI.str where isContained is true
     matchingElements.m1.MTI.idx= find(isContained);

% 3. assign values to MTI model and check which elements are remaining
  



% XX.  find general component names

[matchingElements.components.names,matchingElements.components.idx,t3]=unique(extractBefore(matchingElements.SL.str,"_"));

% %%
% for k=1:length(matchingElements.m1.names)
% 
%     if length(out.xout{matchingElements.m1.idx(k)}.Values.Data(tSimIdx,:))==1
%         xEP(matchingElements.MTI.strIdx(matches(matchingElements.MTI.str,matchingElements.m1.names(k))))=out.xout{matchingElements.m1.idx(k)}.Values.Data(tSimIdx,:);
%     else
%         for l=1:length(matchingElements.m1.MTI.idx)
%             if contains(matchingElements.m1.names(k),matchingElements.MTI.str(matchingElements.m1.MTI.idx(l)))
%                    xEP(matchingElements.m1.MTI.idx(l):matchingElements.m1.MTI.idx(l)+1)=out.xout{matchingElements.m1.idx(k)}.Values.Data(tSimIdx,:);
%             end
%         end 
%     end
% 
% end
%%

%% get indices
% 1. all exactly matching elements, refering to the elements coming from
% the simulink file
matchingElements.SL.exactMatchingIdx=matchingElements.m1.idx(matches(matchingElements.m1.names,matchingElements.MTI.str(matchingElements.m1.MTI.idx)));
matchingElements.SL.exactMatchingNames=matchingElements.m1.names(matches(matchingElements.m1.names,matchingElements.MTI.str(matchingElements.m1.MTI.idx)));

% 2. all remaining elememts, where an exact matching could not be found,
% for example, grid_IDQ 
matchingElements.SL.notExactMatchingIdx=matchingElements.m1.idx(~matches(matchingElements.m1.names,matchingElements.MTI.str(matchingElements.m1.MTI.idx)));
matchingElements.SL.notExactMatchingNames=matchingElements.m1.names(~matches(matchingElements.m1.names,matchingElements.MTI.str(matchingElements.m1.MTI.idx)));

%% transfer data
% 3. storing the elements for the exact matching
for k=1:length(matchingElements.SL.exactMatchingIdx)
    xEP(matchingElements.MTI.strIdx(matches(matchingElements.MTI.str,matchingElements.SL.exactMatchingNames(k))))=out.xout{matchingElements.SL.exactMatchingIdx(k)}.Values.Data(tSimIdx,:);
    xEP_transfer(matchingElements.MTI.strIdx(matches(matchingElements.MTI.str,matchingElements.SL.exactMatchingNames(k))))=1;
end

% 4. all remaining elements like IDQ

remainFillxEPidx=matchingElements.m1.MTI.idx(find(xEP(matchingElements.m1.MTI.idx)==0));

for k=1:length(remainFillxEPidx)       
    xEP(remainFillxEPidx(k):remainFillxEPidx(k)+1)=out.xout{matchingElements.SL.strIdx(contains(matchingElements.SL.str,matchingElements.MTI.str(remainFillxEPidx(k))))}.Values.Data(tSimIdx,:);
    xEP_transfer(remainFillxEPidx(k):remainFillxEPidx(k)+1)=[1;1];
end



%% cos and sin-values

% find the indices of the states of cos/sin of the MTI model
indexCosine=find(contains(matchingElements.MTI.str,'cos'));
indexSine=find(contains(matchingElements.MTI.str,'sin'));

% find the indices of the states of the angles in the nonlinear model
indexTheta=find(contains(matchingElements.SL.str,'theta'));    


for k=1:length(indexCosine)
    for l=1:length(indexTheta)
        if contains(matchingElements.MTI.str(indexCosine(k)),extractAfter(matchingElements.SL.str(indexTheta(l)), "_"))
              % insert value
              xEP(indexCosine(k):indexSine(k))=[cos(out.xout{matchingElements.SL.strIdx(indexTheta(l))}.Values.Data(tSimIdx));sin(out.xout{matchingElements.SL.strIdx(indexTheta(l))}.Values.Data(tSimIdx))];
              %  
              xEP_transfer(indexCosine(k):indexSine(k))=[1;1]; % tick off 
                


        end
    end 
end

%% Check point if all initial values could obtained 

remainingElementsName=matchingElements.MTI.str(find(xEP_transfer==0));
remainingElementsIdx=matchingElements.MTI.strIdx(find(xEP_transfer==0));

if isempty(remainingElementsName)
    display('All initial state values are obtained.')
else 
    error('Not all initial state values could be obtained.')
end


% %%
% 
% %% Com
% xEPcomp=[out.xout{49}.Values.Data(tSimIdx,:),out.xout{33}.Values.Data(tSimIdx,:),cos(out.xout{1}.Values.Data(tSimIdx,:)),sin(out.xout{1}.Values.Data(tSimIdx,:)),out.xout{1}.Values.Data(tSimIdx,:),out.xout{41}.Values.Data(tSimIdx,:),cos(out.xout{42}.Values.Data(tSimIdx,:)),sin(out.xout{42}.Values.Data(tSimIdx,:)),out.xout{42}.Values.Data(tSimIdx,:),out.xout{44}.Values.Data(tSimIdx,:),out.xout{45}.Values.Data(tSimIdx,:),out.xout{37}.Values.Data(tSimIdx,:),out.xout{38}.Values.Data(tSimIdx,:),...
%        out.xout{31}.Values.Data(tSimIdx,:),out.xout{48}.Values.Data(tSimIdx,:),cos(out.xout{35}.Values.Data(tSimIdx,:)),sin(out.xout{35}.Values.Data(tSimIdx,:)),out.xout{34}.Values.Data(tSimIdx,:),out.xout{39}.Values.Data(tSimIdx,:),out.xout{40}.Values.Data(tSimIdx,:)].';
% converter1=squeeze(out.converter1.signals.values);
% converter2=squeeze(out.converter2.signals.values);
% 
%  yEPcomp=[nonlinearSimulationData.VloadDQ(:,tSimIdx).',converter1(:,tSimIdx).',xEPcomp(19:20).',converter2(:,tSimIdx).'];
% 
% 
% 
% initialValuesComparison=table(xEP,xEPcomp,'RowNames',MTImodel.stateName,'VariableNames',{'MTI','Nonlinear'})
% 
% initialValuesComparisonAlgebraic=table(multilinearSimulationData.y(1,:).',yEPcomp.','RowNames',MTImodel.algebraicName,'VariableNames',{'MTI','Nonlinear'})
% 
% 

 
 