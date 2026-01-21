function Jsparse = cpn2LinSparse(obj,x)
% <cpn2LinSparse> - Computes Jacobian around operating point
%
% Input parameter: 
%  - obj: CPN1 object 
%  - x  : operating point 
% 
% Output parameter:
%  - Jsparse: sparse Jacobian matrix
%
% Example:
% Jsparse = cpn2LinSparse(obj,x_op)
%
% For detailed documentation see <a href="matlab:open((which('jacobianDoc.html')))">here</a>

% Enrico Uhlenberg, Marah Engels, Torben Warnecke, Leandro Samaniego,
% Carlos Cateriano Yáñez, Christoph Kaufmann, Leona Schnelle, Georg
% Pangalos, and Gerwald Lichtenberg - 12/06/2024

%% Checking dimensions
% Getting dimensions 
[n,rf]=size(obj.U); % n : number of states
                   % rf: rank

% Checking with operating point
if n ~= length(x)
    error("Dimensions of the operating point argument x do not match with the CPN1 object dimensions.")
end


%% Getting indices of sparse system

    [rfi,cfi,vfi] = find(obj.U); % rfi: row index of the factor 
                                 % cfi: column index of the factor 
                                 % vfi: value of the factor  


%% Logic 
% state
    vf = vfi.*(x(rfi)-sign(vfi))+1;
% end

zf = accumarray(cfi,~vf,[rf 1],@sum);     % Number of zeros in elements
pf = accumarray(cfi,vf,[rf 1],@prod,1);   % Product of nonzero factors

%% Special cases

vf(~vf) = Inf;       % to avoid div/0 
FD = pf.*sparse(cfi,rfi,1./vf,rf,n);   % "Product of others" state matrix

% Special cases of exact one zero factor in indexed elements
rowsf = find(zf==1);
for k = 1:length(rowsf)                       
    ri = rowsf(k);                     % row number where this happens
    ci = rfi((cfi==ri)&(vf==Inf));     % column number of the only nonzero element
    FD(ri,ci) = prod(vf((cfi==ri)&(vf~=Inf)));  % Set it to the product of others 
end

Jsparse=obj.phi * (FD.*obj.U'); %% TODO check dimensions
end % of cpn2LinSparse


