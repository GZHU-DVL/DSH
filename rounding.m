function [f, w] = rounding(x)

%this function performs correct rounding and rounding "in the wrong way"
%(Conway, Sloane, IT-28, No 2, March 1982). Ties are broken so as to give
%preference to points of smaller norm.

if x==0
    f = 0;
    w = 1;
elseif x>0
    if (x - floor(x)) <= 1/2
        f = floor(x);
        w = f + 1;
    else
        f = ceil(x);
        w = f - 1;
    end
else
    if (ceil(x)-x) <= 1/2
        f = ceil(x);
        w = f - 1;
    else
        f = floor(x);
        w = f + 1;
    end
end


