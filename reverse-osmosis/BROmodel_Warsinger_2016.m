function [w_kWhm3_ET_tank,w_kWhm3_PX,w_kWhm3_CCRO,w_kWhm3_RO,w_kWhm3_RO_PX,s_mix_tank,s_mix_CC] = batch_semibatch_continuous(s_f_ppt,RR)
%Emily W. Tow (2016)
%MATLAB code that models specific energy consumption and mixing entropy generation of different
%RO configurations with fixed terminal hydraulic-osmotic pressure difference

%Inputs:
%   s_f_ppt salinity in parts per thousand
%   RR recovery ratio

%Outputs:
%   w_kWhm3_ET_tank specific energy consumption in kWh/m^3 permeate of batch
%       RO configuration with tank, assuming tank is compressed by pumping in
%       permeate
%   w_kWhm3_PX specific energy consumption in kWh/m^3 permeate in batch RO
%       with pressure exchanger
%   w_kWhm3_CCRO specific energy consumption in kWh/m^3 permeate in
%       semi-batch RO
%   w_kWhm3_RO specific energy consumption in kWh/m^3 permeate of RO
%       without pressure exchanger
%   w_kWhm3_RO_PX specific energy consumption in kWh/m^3 permeate or RO
%       with pressure exhanger
%   s_b concentrate salinity in parts per thousand
%   s_mix_tank entopy generation due to mixing in J/K-kg permeate in RO
%       configuration with tank, assuming tank is compressed by pumping in
%       permeate
%   s_mix_CC entopy generation due to mixing in J/K-kg permeate in semi-batch RO

% Required additional functions for NaCl solution properties:
%   NaClDensity(molality, temperature [K], pressure [bar]) outputs solution density in kg/m^3
%   NaClOsmoticPressure(molality, temperature [K], pressure [bar]) outputs osmotic pressure in bar
%   NaClWaterActivity(molality, temperature [K], pressure [bar]) outputs water activity
%   NaClActivityCoefficient(molality, temperature [K], pressure [bar]) outputs mean molal activity coefficient

%inputs; note: also needs to be set in other configurations' functions!
eta_HP=0.8;%high pressure pump efficiency
eta_circ=eta_HP;%circulation pump efficiency
eta_b=eta_HP;%booste pump efficiency
eta_PX=0.96;%pressur exchanger efficiency
RR_n=0.3;%recovery ratio of one pass through the module in batch/semi-batch configurations
DP_l=1;%hydraulic loss in bar of one pass through module [bar]
DP_T=5;%terminal hydraulic-osmotic pressure difference [bar]

%initializing
n=50;%~subdivisions of module
n=2*n+1;%make sure it's oddDP_T=5;%bar
T=20+273.15;%temp. in K
V_p=1;%[m^3] permeate volume (arbitrary)
V_m=(1-RR)/RR*V_p;%module volume
V_t0=V_p;%initial tank volume
W=0;%energy consumption at beginning of cycle = 0
W2=0;%energy consumption at beginning of cycle = 0
S_mix=0;%entropy generation due to mixing at beginning of cycle = 0
s_f=s_f_ppt/1000;%salinity as mass fraction
rho_p=NaClDensity(0,T,0);%permeate density
rho_f=rho_p;%assume constant density for this model
S_t=V_t0*rho_f*s_f;%mass of salt in the tank
w_t=V_t0*rho_f*(1-s_f);%mass of water in the tank
m_t=S_t+w_t;%total tank mass
s_n=ones(1,n)*s_f;%fill module with feed at feed salinity
s_t=s_f;%fill tank with feed
V_n1=V_m/(n-RR_n*(n-1)/2);%calculate volume of first module subdivision
dV_p=V_n1*RR_n;%differential volume of permeate removed from whole module in each step
V_n=ones(1,n)*V_m/n-(-(n-1)/2:1:(n-1)/2)*dV_p/n;%calculate volume of each section of module
m_n=V_n*rho_f;%distribute mass in module
S_n=s_f*m_n;%distribute salt in module
w_n=(1-s_f)*m_n;%distribute water in module
V_p_tot=0;%start with zero permeate produced at beginning of cycle

%modeling system pressure and mixing entropy generation as function of permeate production
%for systems with variable-volume tank (with high pressure tank or pressure exchanger)
while V_p_tot<V_p%for all times in cycle
dm_p=dV_p*rho_p;%mass of permeate produced in step
S_out=S_n(end);%salt mass leaving last section
w_out=w_n(end)-dm_p/n;%water leaving last section
for i=n:-1:2%for all sections of module
    S_n(i)=S_n(i-1);%salt in module moves forward
    w_n(i)=w_n(i-1)-dm_p/n;%water from prev section moves forward (minus permeate from section)
    m_n(i)=S_n(i)+w_n(i);%mass of each section
    s_n(i)=S_n(i)/m_n(i);%salinity of section
end
S_n(1)=s_t*V_n(1)*rho_f;%salt mass in first section comes from tank
w_n(1)=(1-s_t)*V_n(1)*rho_f;%water mass in first section comes from tank
s_n(1)=S_n(1)/(S_n(1)+w_n(1));%mass fraction in first section
S_mix=S_mix+entropy_gen(m_t,w_out+S_out,s_t,S_out/(S_out+w_out));%entropy generation due to mixing of concentrate leaving last section and fluid in tank
S_t=S_t+S_out-S_n(1);%update tank salt mass
w_t=w_t+w_out-w_n(1);%update tank water mass
m_t=S_t+w_t;%update tank mass
s_t=S_t/m_t;%update tank salinity
V_p_tot=V_p_tot+dV_p;%update permeate total
P=NaClOsmoticPressure(MolalityfromSalinity(s_n(end)),T,0)+DP_T;%set pressure to exit osmotic pressure plus driving pressure difference
W=W+P*1e5/eta_HP*dV_p;%update work of cycle [J] in high pressure tank configuration
W2=W2+dV_p*(1-RR_n)/RR_n*(P+DP_l-(eta_PX*(P)))*1e5/eta_b+dV_p/eta_HP*(P+DP_l)*1e5;%update work of cycle [J] in pressure exchanger configuration
end

%total specific energy consumption of cycle (batch processes) [J/m^3 permeate]
w_J=W/V_p_tot+DP_l*1e5/(RR_n*eta_circ)+(1-RR)/RR*DP_l*1e5/eta_circ;%high pressure tank configuration
w_2J=W2/V_p_tot+(1-RR)/RR*DP_l*1e5/eta_HP;%pressure exchanger configuration

%calculating specific energy consumption (SEC) [kWh/m^3 permeate] and total
%mixing entropy generation [kWh/K-kg permeate]
w_kWhm3_ET_tank=w_J/3.6e6;%high pressure tank configuration
w_kWhm3_PX=w_2J/3.6e6;%pressure exchanger configuration
w_kWhm3_RO = ET_RO(s_f_ppt,RR);%continuous, single-stage RO w/o energy recovery
w_kWhm3_RO_PX = ET_RO_PX(s_f_ppt,RR);%continuous, single-stage RO with pressure exchanger
[w_kWhm3_CCRO,s_mix_CC] = CCRO(s_f_ppt,RR);%SEC and entropy generation due to mixing in semi-batch RO (CCRO)
s_mix_tank=S_mix/(V_p_tot*rho_p);%entropy generation due to mixing in batch RO (regardless of tank type)
end

function [w_kWhm3] = ET_RO_PX(s_f_ppt,RR)
%Emily W. Tow (2016)
%Simple continuous RO model with pressure exchanger

%Inputs:
%   s_f_ppt salinity in parts per thousand
%   RR recovery ratio

%Outputs:
%   w_kWhm3 specific energy consumption in kWh/m^3 permeate in RO with
%   pressure exchanger

%inputs; note: also needs to be set in other configurations' functions!
eta_PX=0.96;%pressure exchanger efficiency
DP_T=5;%terminal hydraulic-osmotic pressure difference [bar]
eta_HP=0.8;%high pressure pump efficiency
eta_b=0.8;%circulation pump efficiency
RR_n=0.3;%module recovery ratio
DP_l=1;%pressure drop of one pass through module of batch size [bar]
T=20+273.15;%temp. in K
s_f=s_f_ppt/1000;%salt mass fraction in feed
s_b=s_f/(1-RR);%salt mass fraction of concentrate
passes=log(1-RR)/log(1-RR_n);%equivalent number of passes for pressure drop calculation (N in paper)
pi_b=NaClOsmoticPressure(MolalityfromSalinity(s_b),T,0);%osmotic pressure at brine exit
P=pi_b+DP_T+DP_l*passes;%calculate pressure
w_J=P*1e5/eta_HP+(P-eta_PX*(pi_b+DP_T))*1e5*(1-RR)/RR/eta_b;%SEC
w_kWhm3=w_J/3.6e6;%SEC in kWh/m^3 permeate
end

function [w_kWhm3] = ET_RO(s_f_ppt,RR)
%Emily W. Tow (2016)
%Simple continuous RO model without pressure exchanger

%Inputs:
%   s_f_ppt salinity in parts per thousand
%   RR recovery ratio

%Outputs:
%   w_kWhm3 specific energy consumption in kWh/m^3 permeate in RO without
%   pressure exchanger

DP_T=5;%terminal hydraulic-osmotic pressure difference [bar]
eta_HP=0.8;%high pressure pump efficiency
RR_n=0.3;%module recovery ratio
DP_l=1;%pressure drop of one pass through module of batch size [bar]
T=20+273.15;%temp. in K
s_f=s_f_ppt/1000;%salt mass fraction in feed
s_b=s_f/(1-RR);%salt mass fraction of concentrate
passes=log(1-RR)/log(1-RR_n);%equivalent number of passes for pressure drop calculation (N in paper)
P=NaClOsmoticPressure(MolalityfromSalinity(s_b),T,0)+DP_T+DP_l*passes;%calculate pressure
w_J=P/RR*1e5/eta_HP;%SEC
w_kWhm3=w_J/3.6e6;%SEC in kWh/m^3 permeate
end


function [w_kWhm3,s_mix] = CCRO(s_f_ppt,RR)
%Emily W. Tow (2016)
%Semi-batch RO model
%Inputs:
%   s_f_ppt salinity in parts per thousand
%   RR recovery ratio

%Outputs:
%   w_kWhm3specific energy consumption in kWh/m^3 permeate in
%       semi-batch RO
%   s_mix entopy generation due to mixing in J/K-kg permeate in semi-batch RO

%inputs
DP_T=5;%bar
eta_HP=0.8;
eta_circ=0.8;%high
RR_n=0.3;%.1985;
DP_l=1;%RR_n/0.4;%bar (1 bar at RR_n=0.4)

%initialize
n=50;%~subdivisions of module
n=2*n+1;%make sure it's odd
T=20+273.15;%temp. in K
V_p=1;%[m^3] permeate volume (arbitrary)
V_m=(1-RR)/RR*V_p;%module volume
W=0;%start with zero work for cycle
S_mix=0;%start with zero entopy generated
s_f=s_f_ppt/1000;%salinity as mass fraction
rho_p=NaClDensity(0,T,0);%permeate density
rho_f=rho_p;%approximate density as equal to permeate density (this model)
s_n=ones(1,n)*s_f;%fill module with feed
V_n1=V_m/(n-RR_n*(n-1)/2);%volume of first section of module
dV_p=V_n1*RR_n;%water produced in module in each timestep
V_n=ones(1,n)*V_m/n-(-(n-1)/2:1:(n-1)/2)*dV_p/n;%volume of each module subdivision
m_n=V_n*rho_f;%distribute mass in module
S_n=s_f*m_n;%distribute total salt in module
w_n=(1-s_f)*m_n;%water in module
V_p_tot=0;%zero permeate to start

%calculate pressure and entropy generation
while V_p_tot<V_p%for all time steps
dm_p=dV_p*rho_p;%differential permeate mass
S_out=S_n(end);%salt leaving last section
w_out=w_n(end)-dm_p/n;%water leaving last section
for i=n:-1:2%for all module sections
    S_n(i)=S_n(i-1);%salt in module moves forward
    w_n(i)=w_n(i-1)-dm_p/n;%water from prev section moves forward (minus permeate from section)
    m_n(i)=S_n(i)+w_n(i);%mass of each section
    s_n(i)=S_n(i)/m_n(i);%salinity of section
end
V_f=V_n(1)-(V_n(end)-dV_p/n);%volume of fresh feed added to module during step
S_f=V_f*rho_f*s_f;%mass of salt added to module from feed during step
S_n(1)=S_f+S_out;%total salt mass added to module during step
s_n(1)=S_n(1)/(V_n(1)*rho_f);%salinity of first section
w_n(1)=(1-s_n(1))*V_n(1)*rho_f;%water mass in first section 
S_mix=S_mix+entropy_gen(V_f*rho_f,w_out+S_out,s_f,S_out/(S_out+w_out));%entropy generation due to mixing in first section
V_p_tot=V_p_tot+dV_p;%collect permeate
P=NaClOsmoticPressure(MolalityfromSalinity(s_n(end)),T,0)+DP_T+DP_l;%calculate module pressure
W=W+dV_p*P*1e5/eta_HP;%calculate energy consumption of time step
end

%compute energy consumption and entropy generation
w_J=W/V_p_tot+(1-RR_n)/RR_n*DP_l*1e5/eta_circ+(1-RR)/RR*DP_l*1e5/eta_HP;%SEC in J/m^3
s_mix=S_mix/V_p_tot/rho_p;%entropy generation per unit permeate mass
w_kWhm3=w_J/3.6e6;%SEC in kWh/m^3 permeate
end

function [S_gen] = entropy_gen(m_1,m_2,s_1,s_2)
%Emily W. Tow (2016)
%calculates entropy due to mixing of streams (1 and 2) in tank of bacth processes or
%in tee of semi-batch processes

%inputs
%   m_1,m_2 masses in kg of solution mixed together
%   s_1,s_2 salinities (as mass fraction salt) of solutions mixed together
%output: entropy generation due to mixing in J/K

%constants
R = 8.3145; % Universal gas constant, J/mol-K
M_w = 18.02/1000;%water molar mass kg/mol
M_s=58.44/1000;%NaCl (salt) molar mass kg/mol
T = 293.15;%temp. in K
P=1;%absolute pressure [bar]

%convert to molar
n_s1=m_1*s_1/M_s;%stream 1 moles of salt
n_s2=m_2*s_2/M_s;%stream 2 moles of salt
n_w1=m_1*(1-s_1)/M_w;%stream 1 moles of water
n_w2=m_2*(1-s_2)/M_w;%stream 2 moles of water

%get mixed salinity
m_3=m_1+m_2;%mass of mixed stream
n_s3=n_s1+n_s2;%moles of salt in mixed stream
s_3=(m_1*s_1+m_2*s_2)/m_3;%mass fraction of salt in mixture

%covert to molal
b_1=n_s1/(m_1*(1-s_1));%molality of 1
b_2=n_s2/(m_2*(1-s_2));%molality of 2
b_3=n_s3/(m_3*(1-s_3));%molality of mixture

%activities of salt and water
a_1w=NaClWaterActivity(b_1,T,P);%stream 1
a_2w=NaClWaterActivity(b_2,T,P);%stream 2
a_3w=NaClWaterActivity(b_3,T,P);%mixture
gamma_1s=NaClActivityCoefficient(b_1,T,P);%mean molal activity coefficient of salt in stream 1
gamma_2s=NaClActivityCoefficient(b_2,T,P);%mean molal activity coefficient of salt in stream 2
gamma_3s=NaClActivityCoefficient(b_3,T,P);%mean molal activity coefficient of salt in mixture
a_1s=(gamma_1s*b_1)^2;%activity of salt in stream 1
a_2s=(gamma_2s*b_2)^2;%activity of salt in stream 2
a_3s=(gamma_3s*b_3)^2;%activity of salt in mixture

%compute entropy generation
S_gen=-R*(n_s1*log(a_3s/a_1s)+n_s2*log(a_3s/a_2s)+n_w1*log(a_3w/a_1w)+n_w2*log(a_3w/a_2w));[J/kg-K]
end

function b=MolalityfromSalinity(s)
%input: salinity as a mass fraction
%output: molality
M_s=58.44/1000;%molar mass of NaCl [kg/mol]
b=s/(1-s)/M_s;%molality 
end
