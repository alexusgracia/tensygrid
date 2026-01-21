function dmss = mss2hyDmss(sys)
%MSS2HYDMSS Conversion of mss-class to hyDmss-class
%   converts an explicit multilinear timeinvariant model into an implicit
%   (hybrid) multilinear timeinvariant model

dmss = hyDmss(sys);
end

