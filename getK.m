function [ K ] = getK( T )
TmeltStart=330;
TmeltEnd=335;

if T>TmeltEnd
    K=0.13;
elseif T>TmeltStart
    K=-0.018*T+6.16;
else
    K=0.22;
end

end

