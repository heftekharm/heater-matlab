function [ Cp ] = getCp( T )
TmeltStart=330;
TmeltEnd=335;

if T>TmeltEnd
    Cp=2000;
elseif T>(TmeltStart+TmeltEnd)/2
    Cp=-28800*T+9650000;
elseif T>TmeltStart
    Cp=28800*T-9502000;
else 
    Cp=2000;
end
end

