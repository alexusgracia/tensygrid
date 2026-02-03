% Computing the chain rule 
% KauChr, 29.04.2025
%
%%

% Obtain symbolic Jacobian from 
JMTImodel=jacobian(MTImodel);
%%
% Extracting state-space matrices
E=[JMTImodel.equality.stateDerivative,zeros(MTImodel.n+MTImodel.p,MTImodel.p)];
A=-[JMTImodel.equality.state,JMTImodel.equality.algebraic];
B=-JMTImodel.equality.input;
C=eye(size(-[JMTImodel.equality.state,JMTImodel.equality.algebraic],1));
D=zeros(size(B));

% Getting symbolic variables
syms x [MTImodel.n 1]
syms xp [MTImodel.n 1]
syms yp [MTImodel.p 1]
syms y [MTImodel.p 1]
syms u [MTImodel.m 1]
syms t


% Getting position of recasted variables
indexCosine=find(contains(MTImodel.stateName,'cos'));
indexSine=find(contains(MTImodel.stateName,'sin'));
indexAll=sort([indexCosine, indexSine]);


indexTheta=find(contains(MTImodel.stateName,'theta')); % get all the indices where there is theta
indexTheta=setdiff(indexTheta,indexAll); % take only the indices that are not (!) in indexAll



%% Left Jacobian 
JacobianLeft=A(:,indexAll);

%% Right Jacobian

syms theta [length(indexTheta) 1]

functionVector=[];

for k=1:length(indexTheta)
    functionVector=[functionVector, cos(theta(k)), sin(theta(k))];            
end

JacobianRight=[];
for k=1:length(indexTheta)
    JacobianRight=[JacobianRight,diff(functionVector.',theta(k))];
end


 JacobianRight=subs(JacobianRight,functionVector.',x(indexAll));

%% Compute Jacobian Chain Rule 
JacobianChainRule=JacobianLeft*JacobianRight;


%% Correcting A matrix
A(:,indexAll)=zeros(size(A(:,indexAll)));

for k=1:length(indexTheta)
    A(:,indexTheta(k))=JacobianChainRule(:,k);           
end

%% Get filepath to store the matlabFunction

 %functionPath = filepath(pathMTImodels);%[pathMTImodels,'\demo\Analysis\Applications\PowerSystems\SSA2bus\src\model\MTImodels']; 
 [functionPath,~,~] = fileparts(pathMTImodels);


%% Generate matlabFunction
% in4/u needs to be replaced with 2*pi*50 to reduce the u as an argument of the
% function in case these Jacobians are used for siumlation where function
% handle should have only the arguments (t,x,xp)

outFile = fullfile(functionPath, 'JacobianFunctionA.m'); 
JacobianFunctionA=matlabFunction(A,'Vars',{t,[x; y;],[xp; yp;],u},'File',outFile);

outFile = fullfile(functionPath, 'JacobianFunctionE.m'); 
JacobianFunctionE=matlabFunction(E,'Vars',{t,[x; y;],[xp; yp;],u},'File',outFile);

outFile = fullfile(functionPath, 'JacobianFunctionB.m'); 
JacobianFunctionB=matlabFunction(B,'Vars',{t,[x; y;],[xp; yp;],u},'File',outFile);



