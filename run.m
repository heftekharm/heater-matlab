clear
clc
format shortG
LoopsNumber=input('How many times?');%How mnay times would you like to run the heater script?
t=0;%initial t
T_OFF=65;
T_ON=60;
figure
heater3withPCM(LoopsNumber,T_OFF,T_ON);