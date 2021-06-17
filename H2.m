function hydrogen_oxygen_combustion(g)

%Hydrogen and oxygen combustion
%Mateusz Krasuski
%MKWS 2021 r.


%input data
oneatm = 100000; %ambient pressure [bar];
p_cham=296.1*oneatm; %chamber pressure [atm];
T0=273; %initial temperature [K];
p_e=oneatm; %pressure at the nozzle exit [atm];
M_ful=2; % fuel molecular mass [g/mol];
M_ox=32.0; % oxidizer molecular mass [g/mol];
R=8314.46; %universal gas constant [J/(kmol*K)];
gravity=9.81; %gravity constant;


%finding methane and oxygen indices
gas = Solution('gri30.yaml');
nsp = nSpecies(gas);
ih2 = speciesIndex(gas,'H2');
io2  = speciesIndex(gas,'O2');

%calculations are conducted for nO_F different O/F ratios beggining at o_f_min
%with constant interval equal to o_f_min
nO_F = 16;
o_f_min=0.8;

%initializing arrays
o_f = linspace(o_f_min, o_f_min*nO_F, nO_F);
tad=zeros(nO_F,1);
xeq(nsp,nO_F) = 0;
M_gas=zeros(nO_F,1);
gamma=zeros(nO_F,1);
Isp=zeros(nO_F,1);
Results_CEA=zeros(2*nO_F,3);
Isp_CEA=zeros(nO_F,1);
Temp_CEA=zeros(nO_F,1);

%equilibrium and ISP calculations using Cantera 
for i =1:1:nO_F
   x = zeros(nsp,1);
   x(ih2,1) = 1.0;
   x(io2,1) = o_f(i)*M_ful/M_ox;
   set(gas,'Temperature',T0,'Pressure',oneatm,'MoleFractions',x);
   equilibrate(gas,'HP');
   tad(i) = temperature(gas);
   xeq(:,i) = moleFractions(gas);
   M_gas(i)=meanMolecularWeight(gas);
   gamma(i)=cp_mass(gas)/cv_mass(gas);
   Isp(i)=sqrt(2*gamma(i)*R*tad(i)/((gamma(i)-1)*M_gas(i))*...
   (1-(p_e/p_cham)^((gamma(i)-1)/gamma(i))));
end

%reading results from CEA - temperature, gamma and molecular weight;
Results_CEA = readmatrix('H2_CEA_data.txt');
for i =1:1:nO_F
Isp_CEA(i)=sqrt(Results_CEA(2*i-1,1)*R/Results_CEA(2*i-1,2)*...
    (2*Results_CEA(2*i-1,3))/(Results_CEA(2*i-1,3)-1)*...
    (1-(p_e/p_cham)^((Results_CEA(2*i-1,3)-1)/Results_CEA(2*i-1,3))));
Temp_CEA(i)=Results_CEA(2*i-1,1);
end

%making plots and printing results
figure(1);
subplot(1,2,1);
plot([o_f_min:o_f_min:o_f_min*nO_F],Isp/gravity);
hold on;
plot([o_f_min:o_f_min:o_f_min*nO_F],Isp_CEA/gravity);
hold off;
grid minor;
xlabel('Oxidizer-Fuel Ratio');
ylabel('Specific impulse [s]');


subplot(1,2,2);
plot([o_f_min:o_f_min:o_f_min*nO_F],tad,'DisplayName','CANTERA');
hold on;
plot([o_f_min:o_f_min:o_f_min*nO_F],Temp_CEA,'DisplayName','CEA');
hold off;
grid minor;
xlabel('Oxidizer-Fuel Ratio');
ylabel('Temperature in combustion chamber [K]');
lgd = legend;
lgd.NumColumns = 1;
