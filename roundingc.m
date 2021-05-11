function [f, w] = roundingc(x)

%this function performs correct rounding and rounding "in the wrong way"
%(Conway, Sloane, IT-28, No 2, March 1982). Ties are broken so as to give
%preference to points of smaller norm.



