function [ R ] = ResistancesSum( type,resisdances )
%RESISTANCESSUM Summary of this function goes here
%   Detailed explanation goes here
R=0;
if strcmp(type,'series')
    
    for i=1:length(resisdances)
    R=R+resisdances(i);
    end    
    
elseif strcmp(type,'parallel')
    
for i=1:length(resisdances)
R=R+1/resisdances(i);
end

R=R.^-1;
end



end

