function [ R ] = getResistance(type,args)%if type=="cy" args:[ri,ro,k,l]     if type=="ca" args:[dx,k,r]
if strcmp(type,'cy')
R=log(args(2)/args(1))/(2*pi*args(3)*args(4));
R;
elseif strcmp(type,'ca')
R=args(1)/(args(2)*pi*args(3).^2);
end
end

