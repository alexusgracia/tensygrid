function msys = c2d(obj,Tsample,options)
%function msys = c2d(obj,Tsample)
% Transformation of the continuous time MTI system to its
% discrete time version by first order euler approximation dx =
% (x(k+1) - x(k))/Tsample, with sampling time Tsample.  
% Input parameters: 
%  - obj: mss object
%  - Tsample: discrete sampling time
%
% Output parameter: 
%  - msys: mss object 
%
%
% Example
% msys = mss.c2d(obj,Tsample)
%
%   31.03.2024
%   leona.schnelle@haw-hamburg.de
% For more detailed documentation see: <a href="../../documentation/mss/html/c2dDoc.html">click here</a>
arguments 
    obj
    Tsample
    options.fullrank (1,1) int16 = 1
end

if  options.fullrank
    Fx=diag(ones(obj.n,1));
    Fu=zeros(obj.m,obj.n);
    F.U=[Tsample*obj.F.U [Fx;Fu]];
    F.phi=[obj.F.phi Fx];
else
% low rank TODO: efficiency with avoiding loops 
% if the term x(k) is already in the equation, no rank is added
    F.phi=obj.F.phi;
    F.U=obj.F.U;
    %ozFU=zeros(1,size(obj.F.U,2));
    statenum=zeros(1,size(obj.F.U,2));
    for k=1:size(obj.F.phi,2)
        indnz=find(obj.F.U(:,k));
        if size(indnz,1)==1&& ~any(statenum==indnz)&& indnz<=obj.n
        %ozFU(k)=indnz;
            F.phi(indnz,k)=F.phi(indnz,k)+1;
            statenum(k)=indnz;
        end
    end

%ozFU(:,~all([any(ozFU,1); ozFU<=obj.n]))= 0;
% F.phi(ozFU(any(ozFU,1)),find(ozFU))=F.phi(ozFU(any(ozFU,1)),(find(ozFU)))+1;
%F.phi([ozFU(find(ozFU)) find(ozFU)])=obj.F.phi(ozFU(find(ozFU)),(find(ozFU)))+1;

    missing=find(~ismember([1:obj.n],statenum));
    if missing
        Fx=zeros(obj.n,length(missing));
        Fu=zeros(obj.m,length(missing));
        %Fx(missing,[1:length(missing)])=1; %bug sets too many ones
        for k=1:length(missing)
            Fx(missing(1),k)=1;
        end
        F.U=[Tsample*obj.F.U [Fx;Fu]];
        F.phi=[obj.F.phi Fx];
    end
    
end
msys=mss(CPN1(F.U,F.phi),CPN1(obj.G.U,obj.G.phi),Tsample);
end  % of c2d
