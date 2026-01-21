function [f0, f1, f2, f3] = ktensor2Cpn(T)
%function [f0, f1, f2, f3] = ktensor2Cpn(T)
% <ktensor2Cpn> - Prepares the factors necessary for the
% transformation of ktensor T (@tensortoolbox) into a CPN tensor. 
%   
%   input parameter:
%   - F: [input tensor of class 'ktensor'@tensortoolbox].
%
%   output parameters:
%   - f0: [factors of lambda]
%   - f1: [factors of "1"]
%   - f2: [factors of signals]
%   - f3: [factors of parameters]

%   Example:
% 
%   T=cp_als(tensor(rand(2,2,2,2,3)),4);
%   [f0, f1, f2, f3] = ktensor2Cpn(T);


    MF = cell2mat(T.U);     % Matrix of factors
    N  = size(T.U,1)-1;     % Number of signals
    f0 = T.lambda;          % Lambda vector
    f1 = MF(1:2:2*N,:);     % Factors of "1"
    f2 = MF(2:2:2*N,:);     % Factors of signals
    f3 = MF(2*N+1:end,:);   % Factors of parameters
    zf = all(f1|f2);
    if any(~zf)
        f0=f0(zf,:);
        f1=f1(:,zf);
        f2=f2(:,zf);
        f3=f3(:,zf);
        warning('Zero factors detected, tensor reduced.')
    end
end