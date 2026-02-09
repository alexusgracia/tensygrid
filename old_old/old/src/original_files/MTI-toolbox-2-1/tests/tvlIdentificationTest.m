clear
%create TVL structures for state equations
myTvlStruct1(:,:,1) = [0 0 0; 0 1 0; 0 0 0]; %position of dont cares
myTvlStruct1(:,:,2) = [1 1 0; 1 0 1; 0 0 0]; %Boolean values
myTvlStruct2(:,:,1) = [1 0 0; 1 0 0]; %position of dont cares
myTvlStruct2(:,:,2) = [0 1 1;0 0 0]; %Boolean values
myTvlStruct3(:,:,1) = [1 1 0]; %position of dont cares
myTvlStruct3(:,:,2) = [0 0 0]; %Boolean values

%create  OTVLs
myOTvl1 = otvl(myTvlStruct1);
myOTvl2 = otvl(myTvlStruct2);
myOTvl3 = otvl(myTvlStruct3);

%create otvlmss object that contains structural model information formed by
%TVLs
myC = zeros(3,1); %dimension: number of variables x1
myPhi = ones(6,1); %dimension: total TVL rows x 1
myA = zeros(2,3);
timestep = 1;
oTens = otvlTens(myPhi, myA, myC, myOTvl1, myOTvl2, myOTvl3);
omss = mss(oTens,timestep);
%% exampe otvlsim

%initial states
x0 = [1 1 1];

tEnd = 150;
%simulate state trajectory
myNewX = msim(omss,[],1:tEnd-1,x0);
xcv = 1;
%% example ALS
load("xdata.mat")

%identification with ALS algorithm
[otvlSysID, costID] = mlgreyest(xDataEst,1,oTens,AlsTolerance=5e-3);