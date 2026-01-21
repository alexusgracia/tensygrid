function [U,phi] = ktensor2CPN1(T)
% ktensor2CPN1 calculates the 1-norm normalized cpn factors from a
% ktensor T (@tensortoolbox)
%   input parameter:
%   - T: [input tensor of class 'ktensor'@tensortoolbox].
%
%   output parameters:
%   - U: [factor matrix of cpn tensor]
%   - phi: [phi matrix of cpn tensor]
%   Example:
% 
%   T=cp_als(tensor(rand(2,2,2,2,3)),4);
%   [U,phi] = CPN1.ktensor2CPN1(T);

[f0, f1, f2, f3] = cpnTens.ktensor2Cpn(T);

fn = abs(f1)+abs(f2); 
b=f1<0;
U = f2./fn.*(~b-b);        % switch sign of U if f1 was negative
% normalized Phi matrix for T

phi =[f3.*(f0'.*prod(fn))]; 
odd=mod(sum(b),2)>0;        % TODO more efficient?
phi(:,odd)=-phi(:,odd);     % switch sign of phi if number of negatives in a column of f1 is odd
end