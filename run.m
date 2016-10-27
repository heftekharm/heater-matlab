clear
clc
format shortG
CoolHeatCyclesNum=4%input('How many times?');%How mnay times would you like to run the heater script?
t=0;%initial t
T_OFF=65+273;
T_ON=55+273;
figure
heater3withPCM(CoolHeatCyclesNum,T_OFF,T_ON);