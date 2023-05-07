%Second exercise of Thermal equipment
%written by Srijan Dasgupta and Soroush Rostami
%this code has been written for question 3
clear all;
close all;
clc;


pyversion
[v,e] = pyversion; system([e,' -m pip uninstall -y CoolProp']);
[v,e] = pyversion; system([e,' -m pip install --user -U CoolProp']);

import py.CoolProp.CoolProp.*
fluid='R141b';
i=1
q=1

for j=1:5

%% Input Data
Dt(j)=2.64e-3 + (j-1)*0.1e-3; %section diameter at throat
D1=4.50e-3; %section diameter for the section 1
D2=6.70e-3; %section diameter for section 2
np=0.95; %efficiency
npy=0.88; %efiiciency
ns=0.85; %efficiency
phi_m=0.82;%efficiency  
Tp(q)=95+273.15+(q-1)*15; %stagnation primary temperature
Pp(i)=6.04e5+(i-1)*0.25*1e5; %primary pressure
Ts=8+273.15; %stagnation secondary temperature
Ps=0.40e5; %secondary intial pressure
Mpt=1; %mach number at primary throat
Msy=1; %mach number at secondary throat 
%% Coolprop to find gamma and Rg 
Cp=PropsSI('Cp0mass','P',Pp(i),'T',Tp(q),fluid); %J/kg/K (ideal gas mas specific constant pressure specific heat
Cv=PropsSI('Cvmass','P',Pp(i),'T',Tp(q),fluid); %couldn't find ideal Cv value 
%Rg=Cp-Cv; %probably incorrect
M=PropsSI('molarmass',fluid); %molar mass 
R=PropsSI('gas_constant', fluid); %univeral gas constant 
Rg=R/M; %specific gas constant 
gamma=Cp/(Cp-Rg); %isentropic factor 
%% Area of the sections 
At=(pi*Dt(j)^2)/4;
Ap1=(pi*D1^2)/4;
A2=(pi*D2^2)/4;

%% Mass of Primary flow
massflow_p=sqrt(gamma/Rg)*(Pp(i)/sqrt(Tp(q)))*Mpt*(1+(gamma-1)*(Mpt^2)/2)^(-(gamma+1)/(2*(gamma-1)))*At*sqrt(np);

%% Section 1-1 
syms mp1 
equa01=(Ap1/At)==(Mpt*(1+((gamma-1)/2)*Mpt^2)^(-(gamma+1)/(2*(gamma-1))))...
    /(mp1*(1+((gamma-1)/2)*mp1^2)^(-(gamma+1)/(2*(gamma-1))));
Mp1=double(vpasolve(equa01,mp1)); 
Pp1=Pp(i)*((1+((gamma-1)/2)*Mp1^2))^(-gamma/(gamma-1));
%% Section y-y 
Psy=Ps*((1+((gamma-1)/2)*Msy^2))^(-gamma/(gamma-1)); %secondary pressure at y section 
Ppy=Psy;

syms mpy
equa02=(Ppy/Pp1)==(((1+((gamma-1)/2)*Mp1^2)^(gamma/(gamma-1)))/((1+((gamma-1)/2)*mpy^2)^(gamma/(gamma-1))));
Mpy=double(vpasolve(equa02,mpy));
Msy=1;
%temperature of both fluids in the section y
Tpy=Tp(q)*(1+(gamma-1)*Mpy^2/2)^-1;
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
%I have doubt about it
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


E_r(j)=massflow_s/massflow_p;
C_r(j)=Pct/Ps;


end

plot(Dt,E_r)
xlabel('Throat Diameter (m)')
ylabel('Entaintment ratio')
figure
plot(Dt,C_r)
xlabel('Throat Diameter (m)')
ylabel('Compression ratio')

