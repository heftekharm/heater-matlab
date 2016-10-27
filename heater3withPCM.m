%{
In this script a heater covered with a layer of PCM round cylinderical water
tank has been simulated.The script has two main part :heating and
cooling. Each one has divided into two section.330 degree determines the
boundary between sections.Also PCM layer is divided into some sublayers.
In heating part.The water Tempratue increases until
the temprature equals to 330 degree, after that (second ection) if a
sublayer temprature be more than the degree ,that sublayer will be added to
the melting sublayer.That means mass*cp of that sublayer  is added to total
mass*cp (water+pcm).on the other hand this process results in decrasing the
velocity of elevating temprature.
At second part(cooling part) this process occures but in opposite direction.
%}

%Tempratures are in Kelvin
%Lengths,Heights and Radiuses are in meter

function [] = heater3withPCM(CoolHeatCyclesNum,T_OFF,T_ON )

T0=298;     %The air and initial Temprature in Kelvin
L1=0.6;     %The height of the water tank in meter

drPCM=0.02;     %PCM thikness in meter
r1=0.1528;      %water tank radius in meter
rp=r1+drPCM;    %total radius in meter
r3=0.1862;      %outer surface radius in meter

%PCM
Tmelt=333;      %(55C) The melting point of PCM in kelvin
LHeatperMass=180000;        %Latent heat per mass in J/Kg
Vol=((pi*(rp^2))-(pi*(r1^2)))*L1;       %PCM volume in m3
dens=900;       %Density in kg/m3
Mass=Vol*dens;      %PCM mass in kg
LHeat=LHeatperMass*Mass;        %Total latent heat in J
PCM=struct('TmeltStart',328,'TmeltEnd',333,'density',dens,'LHeatperMass',180000,'TotalLHeat',LHeat,'Volume',Vol,'Mass',Mass,'Ksolid',0.22,...
    'Kliquid',0.13,'Cp',2000);   %     K    , kg/m3 ,    J/kg     ,    J    ,   m3    ,  m3  ,  kg   , W.m/K ,  W.m/K  ,  J/kg.K

%Air
kcyAir=0.5;     %effective k cylender
kcaAir=0.5;     %effective k cartesian
dxAir=0.05;     %the distance between inner tank and outer surface at top and below of tank



R1pcy=getResistance('cy',[r1,r1+drPCM,PCM.Ksolid,L1]);%The resistance of PCM in cylinderical direction
Rp3cy=getResistance('cy',[r1+drPCM,r3,kcyAir,L1]);%The resistance of the air trapped in the heater in cylinderical direction

R13ca=getResistance('ca',[dxAir,kcaAir,r1]);%The resistance of the air trapped in the heater in cartesian direction

R13cy=ResistancesSum('series',[R1pcy,Rp3cy]);%The total resistance between the tank and the outer surface in  cylinderical direction
R13=ResistancesSum('parallel',[R13cy,R13ca,R13ca]);%The total resistance between the tank and the outer surface


%R infinity
hinf=1.947;
L3=0.7;
Rinf=1/(hinf*(2*pi*r3*L3+2*pi*rp^2));
Rt=R13+Rinf;

%water
CpWater=4188;%j/kg.k
vWater=pi*(r1.^2)*L1;%m3;
densityWater=980;%kg/m3
mCpWater=densityWater*vWater*CpWater;
alpha=1/(mCpWater*Rt);
Qin=2000;
c=(Qin/mCpWater)+(T0/(mCpWater*Rt));

dt=10;%in second;
t=0;
i=1;

T1i1=T0;T1i=T0;Tpi=T0;T3i=T0; %T in all parts at t=0 equals to T0

PCMisSolid=1;

%{
the following code calculates the tempratures of the water ,the outer surface of
pcm and the outer surface of the heater aare in steady state conditions and when
the next temprature of the water(T1i1) goes to above the pcm melting point
the loop stops.
%}

while T1i1<Tmelt
    T1i=T1i1;
    Q_PCM=(T1i-T0)/Rt;
    T3i=T1i-Q_PCM*R13;
    q1p=(T1i-T3i)/R13cy;
    Tpi=T1i-q1p*R1pcy;
    Temps(i,:)=[T1i,Tpi,T3i];
    tmin=t/60;
    time(i)=tmin;
    LastI=i;
    
    T1i1=T1i*(1-alpha*dt)+c*dt;
    i=i+1;
    t=t+dt;
end

PCM_Layers=10;%the number of pcm sublayers
%meltingArray=zeros(1,PCM_Layers);
%meltedArray=zeros(1,PCM_Layers);
for i=1 : PCM_Layers %temprature and radius for each sublayer is calculated in this loop
    T_PCMLayers(i) = T1i-i*(T1i-Tpi)/PCM_Layers;
    r_PCMLayers(i) = (drPCM/PCM_Layers)*i+r1;
end



dr_Layers=r_PCMLayers(2)-r_PCMLayers(1);
%meltedLayers=0;
meltingLayers=0;
Counter=1;%the next while loop counter
time2(Counter)=time(LastI);
Temps2(Counter,:)=[T1i,Tpi];
Rinf_ca=1/(hinf*(pi*rp^2));
Rinf_cy=1/(hinf*(2*pi*r3*L3));
Counter=Counter+1;
for cn=1:CoolHeatCyclesNum
    while 1
        
        %Water
        Qin=2000;
        R_1_infinity_ca=R13ca+Rinf_ca; %the resistance between the water tank and the surronding air in cartesian direction
        
        R_PCMLayers_total=0;
        R_PCMLayers=[];
        for i=1:PCM_Layers
            R_PCMLayers(i)=getResistance('cy',[r_PCMLayers(i)-dr_Layers,r_PCMLayers(i),getK(T_PCMLayers(i)),L1]);%the resistance of each sublayer is calculated except the melted sublayers
            R_PCMLayers_total=ResistancesSum('series',[R_PCMLayers_total,R_PCMLayers(i)]);%the total resistance of the unmelted sublayers
        end
        R_1_infinity_cy=ResistancesSum('series',[R_PCMLayers_total,Rp3cy,Rinf_cy]);
        
        R_1_infinity=ResistancesSum('parallel',[ R_1_infinity_cy, R_1_infinity_ca,R_1_infinity_ca]);
        
        QoutTotal=(T1i-T0)/R_1_infinity;
        
        meltingl=meltingLayers;
        %Tpi=T1i;
        for i=1:PCM_Layers %the tempratures of all sublayers is calculated
            
            Q_PCM=(T1i-T0)/R_1_infinity_cy;
            
            if i==1
                T_PCMLayers(i)=T1i-Q_PCM*R_PCMLayers(i);
            else
                T_PCMLayers(i)=T_PCMLayers(i-1)-Q_PCM*R_PCMLayers(i);
            end
            if i==PCM_Layers
                Tpi =T_PCMLayers(i);
            end
            
            if T_PCMLayers(i)> Tmelt%if the temprature of a layer is greater than the pcm melting point ,that sublayers is added to the melting sublayers
                meltingl=i;
                meltingArray(i)=1;
            end
            
        end
        
        meltingLayers=meltingl;
        mCpPCM=0;
        for i=1:meltingLayers %this block calculates  the mcp of each sublayer then adds to the PCM mcp total.
            if i==1
                vPCM_Layer=pi*(r_PCMLayers(i).^2)*L1 - pi*(r1.^2)*L1;%m3;
            else
                vPCM_Layer=pi*(r_PCMLayers(i).^2)*L1 - pi*(r_PCMLayers(i-1).^2)*L1;%m3;
            end
            mCpPCM_Layer=PCM.density*vPCM_Layer*getCp(T1i);
            mCpPCM=mCpPCM + mCpPCM_Layer;%the PCM mcp total
        end
        
        mCptotal=mCpWater+mCpPCM;
        T1i1=T1i+(dt/mCptotal)*(Qin-QoutTotal);%the next water temprature
        T1i=T1i1;
        if T1i>T_OFF
            break;
        end
        
        Temps2(Counter,:)=[T1i,Tpi];
        time2(Counter)=t/60;
        Counter=Counter+1;
        t=t+dt;
        
    end
    
    %Cooling section
    while 1
        
        %Water
        Qin=0;%the heater is turned off
        R_1_infinity_ca=R13ca+Rinf_ca; %the resistance between the water tank and the surronding air in cartesian direction
        
        R_PCMLayers_total=0;
        for i=1:PCM_Layers%calculating the resistances of all sublayers
            R_PCMLayers(i)=getResistance('cy',[r_PCMLayers(i)-dr_Layers,r_PCMLayers(i),getK(T_PCMLayers(i)),L1]);%the resistance of each sublayer is calculated
            R_PCMLayers_total=ResistancesSum('series',[R_PCMLayers_total,R_PCMLayers(i)]);%the total resistance of the sublayers
        end
        R_1_infinity_cy=ResistancesSum('series',[R_PCMLayers_total,Rp3cy,Rinf_cy]);%the resistance between the water tank and the surronding air in cylindrical dir
        
        R_1_infinity=ResistancesSum('parallel',[ R_1_infinity_cy, R_1_infinity_ca,R_1_infinity_ca]);%the total resistance between the water tank and surronding air
        
        QoutTotal=(T1i-T0)/R_1_infinity;%the heat that leaves the heater
        
        meltingl=meltingLayers;
        
        for i=1:PCM_Layers %the tempratures of all sublayers is calculated
            
            Q_PCM=(T1i-T0)/R_1_infinity_cy;
            
            if i==1
                T_PCMLayers(i)=T1i-Q_PCM*R_PCMLayers(i);
            else
                T_PCMLayers(i)=T_PCMLayers(i-1)-Q_PCM*R_PCMLayers(i);
            end
            if i==PCM_Layers
                Tpi =T_PCMLayers(i);%the outer surface of pcm
            end
            
            if T_PCMLayers(i)< Tmelt%if the temprature of a sublayer is less than the pcm melting point ,that sublayer will leave he melted sublayers
                meltingl=i-1;
            end
            
            
        end
        
        meltingLayers=meltingl;
        mCpPCM=0;
        for i=1:meltingLayers %this block calculates  the mcp of each sublayer then adds to the PCM mcp total.
            if i==1
                vPCM_Layer=pi*(r_PCMLayers(i).^2)*L1 - pi*(r1.^2)*L1;%m3;
            else
                vPCM_Layer=pi*(r_PCMLayers(i).^2)*L1 - pi*(r_PCMLayers(i-1).^2)*L1;%m3;
            end
            mCpPCM_Layer=PCM.density*vPCM_Layer*getCp(T1i);
            mCpPCM=mCpPCM + mCpPCM_Layer;%the PCM mcp total
        end
        
        mCptotal=mCpWater+mCpPCM;
        T1i1=T1i+(dt/mCptotal)*(Qin-QoutTotal);%the next water temprature
        T1i=T1i1;
        if (T1i<Tmelt) 
            break;
        end
        
        Temps2(Counter,:)=[T1i,Tpi];
        time2(Counter)=t/60;
        
        t=t+dt;
        Counter=Counter+1;
    end
    
    
    Qin=0;
    c=(Qin/mCpWater)+(T0/(mCpWater*Rt));
    
    %this block reduces the tempratures
    while (T1i1>T0+0.1 && cn == CoolHeatCyclesNum) || (T1i1>T_ON )
        T1i=T1i1;
        Q_PCM=(T1i-T0)/Rt;
        T3i=T1i-Q_PCM*R13;
        q1p=(T1i-T3i)/R13cy;
        Tpi=T1i-q1p*R1pcy;
        Temps2(Counter,:)=[T1i,Tpi];
        time2(Counter)=t/60;
        t=t+dt;
        Counter=Counter+1;
        T1i1=T1i*(1-alpha*dt)+c*dt;
    end
end
plot(time,Temps);
hold on;
plot(time2,Temps2);
legend('Water','PCM_{out}','Shell','Location','SouthEast')
grid on
end