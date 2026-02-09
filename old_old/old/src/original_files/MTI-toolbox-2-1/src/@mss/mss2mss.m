function [mtiSystemTransformed] = ...
    mss2mss(mtiSystem, transformation, offset)
% 
% <mss2mss> - Linear state transformation of an explicit MTI model in CPN 
%             format
%   
%
% Input parameters: 
%
%  - mtiSystem: mss object 
%               (with number of states n and number of inputs m). 
%
%  - transformation: transformation matrix 
%                    This must be square and can be n x n, in which case 
%                    only the states are transformed; or (n + m) x (n + m), 
%                    in which case both states and inputs are transformed. 
%                    Currently, only scaling and parameter permutations are
%                    allowed, i.e. there can be only 1 non-zero entry per 
%                    row and column. A diagonal matrix results in scaling. 
%
%  - offset: vector of offsets 
%            (optional, defaults to zero). This must be of length n or 
%            (n + m), matching the dimensions of transformation. 
%
%
% Output parameter: 
%
%  - mtiSystemTransformed: mss object (transformed). 
%
%
% Example
% sys_t = mss2mss(sys, T, c)
%
%
% For detailed documentation see <a href="matlab:open((which(...
% 'mss2mssDoc.html')))">here</a>

% Leona Schnelle - 31/03/2024
% Adapted by Sarah Casura - February 2025


% Define arguments and their defaults
arguments
         mtiSystem
         transformation
         offset = zeros(size(transformation, 1), 1) % default offset to 0
end


% Basic checks for function to work ---------------------------------------

% Dimension consistency checks
% check size of offset and transformation makes sense relative to what is 
% in mtiSystem. 
% TODO check if this needs update after CPN is redefined.

if size(transformation, 1) ~= size(transformation, 2) 
    error('transformation matrix must be square.')
end

if size(transformation, 1) ~= mtiSystem.n && ...
        size(transformation, 1) ~= (mtiSystem.n + mtiSystem.m)
    error(['transformation matrix must be of dimension n x n ', ...
        'or (n + m) x (n + m) \n (where n is the number of states ', ...
        'and m the number of inputs of the mss object)'])
end

% for the offset, we expect a n x 1 vector [or (n + m) x 1]
% but also allow 1 x n vector (just flip if correct length)
if size(offset, 1) ~= size(transformation, 1) 
    if size(offset, 2) == size(transformation, 1) & size(offset, 1) == 1
        offset = offset'; 
    else
        error(['length of offset vector must match size of ', ...
            'transformation matrix'])
    end
end

% Check that transformation matrix has only one entry per row and column
% (we do this by ensuring that no row or column contains zeros only and 
% there are the right number of elements in total). 
isAnyColumnZero = any(sum(transformation, 1) == 0);
isAnyRowZero = any(sum(transformation, 2) == 0);
isCorrectNumberElements = (nnz(transformation) == size(transformation, 1));
if isAnyColumnZero | isAnyRowZero | ~isCorrectNumberElements
    error(['transformation matrix can have only 1 non-zero element ', ...
        'per row and column (allowing only parameter permutation and ', ...
        'scaling). General linear transformations not implemented yet.'])
end

% Check for type
if ~isa(mtiSystem.F, 'CPN1') % TODO update this once CPN is redefined
    error(['transformation in given mti object type not implemented ', ...
        'yet. F.U (and optional G.U) must be CPN1!'])
end



% set G in case of no output equation -------------------------------------
% TODO update this too after new CPN
% I want to avoid modifying the original object (mtiSystem) in-place. 
% Taking a copy of it doesn't do the trick, because it's a handle, so it 
% points to the same object and still modifies the original even if I just 
% work on the copy. Hence I need to create a new temporary object entirely. 
% TODO solve this more elegantly (e.g. by initializing output earlier)
temporaryG = mtiSystem.G; % here I only copy G, not the entire handle... 
if size(mtiSystem.G, 1) == 0
    G.U = [diag(ones(mtiSystem.n, 1)); zeros(mtiSystem.m, mtiSystem.n)];
    G.phi = diag(ones(mtiSystem.n, 1));
    temporaryG = CPN1(G.U, G.phi);
end




% transformation of mss state equation ------------------------------------
% note relative to Leonas original version: a' = diag(transformation) and
% b' = offset

sizeT = size(transformation, 1); % this is the size of the matrix, which
% tells me whether or not we transform the input so I don't need an if-else 
% clause. 
% note this assumes that F.U, G.U always has states first, then inputs!!!
% (which is what Leona assumed. TODO check that assumption with new CPN). 

Fsc.U = mtiSystem.F.U;
% first values of Fx and Fu: 
firstrow = 1 - abs(mtiSystem.F.U(1:sizeT, : )) + ...
    mtiSystem.F.U(1:sizeT, : ) .* offset;
% second values of Fx and Fu: 
secondrow = mtiSystem.F.U(1:sizeT, : ) .* diag(transformation); 
norm1 = abs(firstrow) + abs(secondrow);    % norm 1 of both values
sign1 = firstrow < 0;                      % check for negative values 

% build norm 1 representation of F
Fsc.U(1:sizeT, : ) = secondrow ./ norm1 .* (~sign1 - sign1);
% add column for absolut c substraction from T^-1(x-c) 
% TODO (Leona): check, if zero row for absolut value already exists in 
% F.U and add
Fsc.U(1:sizeT, end+1) = 0; 

Fsc.phi = transformation(1:mtiSystem.n, 1:mtiSystem.n) \ ...
    mtiSystem.F.phi .* prod(~sign1 - sign1) .* prod(norm1);
Fsc.phi( : , end+1) = transformation(1:mtiSystem.n, 1:mtiSystem.n) \ ...
    -offset(1:mtiSystem.n);   
% note T \ X is the same as inv(T) * X which works here because T is 
% diagonal. Not sure whether this will work in general... 
% note also that in this case, we *do* need the cutting of T to 1:n for the 
% case when the input is transformed (it will have no effect when the input
% is not transformed, since T is already n x n then).

% transformation of output equation
% note this repeats quite a few lines from above. 
% TODO once implementing new CPN check how to condense this further. 
Gsc.U = temporaryG.U;
firstrow = 1 - abs(temporaryG.U(1:sizeT, : )) + ...
    temporaryG.U(1:sizeT, : ) .* offset;
secondrow = temporaryG.U(1:sizeT, : ) .* diag(transformation);
norm1 = abs(firstrow) + abs(secondrow);
sign1 = firstrow<0;
Gsc.U(1:sizeT, : ) = secondrow./ norm1 .* (~sign1 - sign1);
Gsc.phi = temporaryG.phi .* prod(~sign1-sign1) .* prod(norm1);


% Build scaled mss object
 mtiSystemTransformed = mss(CPN1(Fsc.U, Fsc.phi), CPN1(Gsc.U, Gsc.phi), ...
     mtiSystem.ts);

end  % of mss2mss
