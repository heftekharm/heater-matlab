clear
clc
format shortG
LoopsNumber=input('How many times?');%How mnay times would you like to run the heater script?
t=0;%initial t
figure
for i=1:LoopsNumber
    [time,Temps,time2,Temps2,t]=heater3withPCM(t);%the function returns 4 outputs
    plot(time,Temps);
    hold on
    plot(time2,Temps2);
    hold on
end