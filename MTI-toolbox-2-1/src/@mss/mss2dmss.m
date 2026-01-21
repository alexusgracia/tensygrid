function sysOut = mss2dmss(sysIn)
%MSS2DMSS Conversion of mss-class to hyDmss-class
%   converts an explicit multilinear timeinvariant model into an implicit
%   (hybrid) multilinear timeinvariant model

sysOut = dmss(sysIn);
end

