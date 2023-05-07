%Second exercise of Thermal equipment
%written by Srijan Dasgupta and Soroush Rostami
%This code has been wriiten for question 2
clear all;
close all;
clc;


pyversion
[v,e] = pyversion; system([e,' -m pip uninstall -y CoolProp']);
[v,e] = pyversion; system([e,' -m pip install --user -U CoolProp']);

import py.CoolProp.CoolProp.*
fluid='R141b';

%% Input Data
Dt=2.64e-3; %section diameter at throat
D1=4.50e-3; %section diameter for the section 1
D2=6.70e-3; %section diameter for section 2
np=0.95; %efficiency
npy=0.88; %efiiciency
ns=0.85; %efficiency
phi_m=0.82;%efficiency  
Tp=95+273.15; %stagnation primary temperature
Pp=6.04e5; %primary pressure
Ts=8+273.15; %stagnation secondary temperature
Ps=0.40e5; %secondary intial pressure
Mpt=1; %mach number at primary throat
Msy=1; %mach number at secondary throat 
%% Coolprop to find gamma and Rg 
Cp=PropsSI('Cp0mass','P',Pp,'T',Tp,fluid); %J/kg/K (ideal gas mas specific constant pressure specific heat
Cv=PropsSI('Cvmass','P',Pp,'T',Tp,fluid); %couldn't find ideal Cv value 
%Rg=Cp-Cv; %probably incorrect
M=PropsSI('molarmass',fluid); %molar mass 
R=PropsSI('gas_constant', fluid); %univeral gas constant 
Rg=R/M; %specific gas constant 
gamma=Cp/(Cp-Rg); %isentropic factor 
%% Area of the sections 
At=(pi*Dt^2)/4;
Ap1=(pi*D1^2)/4;
A2=(pi*D2^2)/4;

%% Mass of Primary flow
massflow_p=sqrt(gamma/Rg)*(Pp/sqrt(Tp))*Mpt*(1+(gamma-1)*(Mpt^2)/2)^(-(gamma+1)/(2*(gamma-1)))*At*sqrt(np);

%% Section 1-1 
syms mp1 
equa01=(Ap1/At)==(Mpt*(1+((gamma-1)/2)*Mpt^2)^(-(gamma+1)/(2*(gamma-1))))...
    /(mp1*(1+((gamma-1)/2)*mp1^2)^(-(gamma+1)/(2*(gamma-1))));
Mp1=double(vpasolve(equa01,mp1)); 

Pp1=Pp*((1+((gamma-1)/2)*Mp1^2))^(-gamma/(gamma-1));

%% Section y-y 
Psy=Ps*((1+((gamma-1)/2)*Msy^2))^(-gamma/(gamma-1)); %secondary pressure at y section 
Ppy=Psy;

syms mpy
equa02=(Ppy/Pp1)==(((1+((gamma-1)/2)*Mp1^2)^(gamma/(gamma-1)))/((1+((gamma-1)/2)*mpy^2)^(gamma/(gamma-1))));
Mpy=double(vpasolve(equa02,mpy));
Msy=1;
%temperature of both fluids in the section y
Tpy=Tp*(1+(gamma-1)*Mpy^2/2)^-1;
Tsy=Ts*(1+(gamma-1)*Msy^2/2)^-1;
%sound speed 
Cpy=sqrt(gamma*Rg*Tpy);
Csy=sqrt(gamma*Rg*Tsy);
%speed of the both fluids
Vsy=Msy*Csy;
Vpy=Mpy*Cpy;
%area of the both sections
Apy=Ap1*npy*(Mp1*(1+((gamma-1)/2)*Mp1^2)^(-(gamma+1)/(2*(gamma-1))))...
    /(Mpy*(1+((gamma-1)/2)*Mpy^2)^(-(gamma+1)/(2*(gamma-1))));
Asy=A2-Apy;
%mass flow of the secondary flow

massflow_s=sqrt(gamma/Rg)*(Ps/sqrt(Ts))*Msy*(1+(gamma-1)*(Msy^2)/2)^(-(gamma+1)/(2*(gamma-1)))*Asy*sqrt(np);

%% section m-m
Vm=phi_m*(massflow_p*Vpy+massflow_s*Vsy)/(massflow_s+massflow_p)

syms tm
eq03=massflow_p*(Cp*Tpy+(Vpy^2)/2)+massflow_s*(Cp*Tsy+(Vsy^2)/2)==(massflow_s+massflow_p)*(Cp*tm+(Vm^2)/2);
Tm=double(vpasolve(eq03,tm));
Mm=Vm/sqrt(gamma*Rg*Tm);
Pm=Psy;

%%section 2 (After shock)
P2=Pm*(1+(2*gamma)*(Mm^2-1))/(gamma+1);
M2=sqrt((1+((gamma-1)*(Mm^2)/2))/((gamma*Mm^2)-(gamma-1)/2));
T2=Tm;
C2=sqrt(gamma*Rg*T2);
V2=M2*C2;
%outlet section
Pct=P2*(1+(gamma-1)*M2^2/2)^(gamma/(gamma-1));
Vct=0;






%%final calculation and ploting

%creating matrix for the nominal point plot
Pp_mat=[Pp,Pp1,Ppy,Pm,P2,Pct];
matX1=[0,2,3,4,4.1,5];

Ps_mat=[Ps,Psy,Pm,P2,Pct];
matX2=[0,3,4,4.1,5];

Mp_mat=[0,Mpt,Mp1,Mpy, Mm, M2,0 ];
matx3=[0,1,2,3,4,4.1,5];
Ms_mat=[0,Msy,Mm,M2,0];
matx4=[0,3,4,4.1,5];

plot(matX1,Pp_mat)
hold on
plot(matX2,Ps_mat)
legend('Primary pressure','Secondary pressure')
xlabel('0:inlet, 1: prymary throat, 2: after the prymary throat, 3: Secondary throat, 4:After mixing, 4.1: After a shockwave, 5: outlet')
ylabel('Pressure (pa)')


figure
plot(matx3,Mp_mat)
hold on
plot(matx4,Ms_mat)
legend('Primary Mach number','Secondary Mach number')
xlabel('0:inlet, 1: prymary throat, 2: after the prymary throat, 3: Secondary throat, 4:After mixing, 4.1: After a shockwave, 5: outlet')
ylabel('Mach number')

