%%  Initialize Variables

%%  Constants
R = 8.3145;                 %[J/Kmol]
g = 9.81;                   %[m/s2]
Patm = 101000;              %[Pa]

%%  Oxidizer System Input Variables
%   The current oxidizer system assumes supply pressure
%   is directly above the injector. So no valve actuation
%   or flow inductance is modeled.

Psup = 5.516e+6;         %[Pa]
%Rval = 1;       % Control Variable
Tsup = 300;              %[K]
%Lpipe = 1;               %[m]
%Apipe = 0.02;            %[m2]
Mox = 0;        %   [kg/s]

%   Injector Properties
Aox = 0.000012566;        %  [m^2]
Cdox = 0.81;    %

%%  Oxidizer System Calculated Parameters
rho_sup = Psup/(R*Tsup);

%%  Combustion Chamber Input Variables
%   Structure
Lcomb = 0.381;  %   [m]
Lpost = 0.0381; %   [m]
Lpre = 0.0254;  %   [m]
Din = 0.0762;   %   [m]
Rg = 0.0127;    %   [m] State Variable
Vcomb = pi()*(Rg^2)*Lcomb+(Lpost+Lpre)*pi()*(Din/4)^2;

%   Fuel Grain Properties
r = 0.0016;             %  [m/s]
rho_fuel = 900;         %  [kg/m^3]
y = 1.4;               %  [J/kg*K] 
Mfuel = 0;  %   [kg/s]

%   Operating Conditions of Combustion Chamber
Pcomb = 0;   %   [Pa]
Tcomb = 1200;              %  [K]


%%   Nozzle Properties and Exhaust Plume
A_thr = 0.0002;             %  [m^2]
exit_theta = 15;            %   [deg]
lamda = 0.5*(1+cosd(exit_theta));
ratio = 70;
Aexit = A_thr * ratio;
Mprop = Mox + Mfuel;

%   Thrust Equations
Pexit = 0;
Texit = 300;
FVac = 0;
FSL = 0;

%%  Model Run Parameters
t = 0.0001;
time = 20;
data = zeros(time/t,12);
%data(1,:) = [i Pcomb Rg Mox Mfuel Mprop r FVac Pexit Texit Vcomb Aburn];


%%  Run Model
for i = 0:t:time
    
    %%  Relation Equations
    Aburn = 2*pi()*Rg*Lcomb;
    Vcomb = pi()*(Rg^2)*Lcomb+(Lpost+Lpre)*pi()*(Din/4)^2;
    Mprop = Mox + Mfuel;
    OF_ratio = Mox/Mfuel;
    Pexit = Pcomb/((1+((y-1)/2)*(t*Mprop)^2)^(y/(y-1)));
    Texit = Tcomb/(1+((y-1)/2)*(t*Mprop)^2);
    FVac = lamda*Mprop*Mprop*t*sqrt(y*Rg*Texit)+Aexit*Pexit;
    %FSL = FVac - Aexit*Patm;
    
    %%  Governing Equations
    dPcomb = ((Aburn*r)/Vcomb)*(rho_fuel*Rg*Tcomb-Pcomb)...
          -Pcomb*((A_thr/Vcomb)*sqrt(y*Rg*Tcomb*(2/(y+1))^((y+1)/(y-1))))...
          +(Rg*Tcomb/Vcomb)*Aox*Cdox*sqrt(2*rho_sup*(Psup-Pcomb));
    dRg = r;
    Mox = Aox*Cdox*sqrt(2*rho_sup*(Psup-Pcomb));
    Mfuel = rho_fuel*Aburn*r;
    
    Pcomb = dPcomb * t + Pcomb;
    Rg = dRg * t + Rg;
    
    if Rg >= Din/2
        r = 0;
        Rg = Din/2;
        Mox = 0;
        Psup = 0;
        Aburn = 0;
    end
    data(uint32((i/t)+1),:) = [i Pcomb Rg Mox Mfuel Mprop r FVac Pexit Texit Vcomb Aburn];
end

%%  Data Plotting
figure(2)
subplot(2,2,1)

plot(data(:,1), data(:,2),'-',data(:,1),data(:,9),'--','Linewidth',2)
title('Combustion Chamber and Exit Pressure', 'Fontsize', 24)
legend('Combustion Chamber', 'Nozzle Exit', 'Supply')
xlabel('Time [s]', 'Fontsize', 24)
ylabel('Pressure [Pascal]', 'Fontsize', 24)

figure(2)
subplot(2,2,2)

plot(data(:,1), data(:,4),'-',data(:,1),data(:,5),'--',data(:,1),data(:,6),':','Linewidth',2)
title('Mass Flow Rates', 'Fontsize', 24)
legend('Oxidizer', 'Fuel', 'Exhaust Gas')
xlabel('Time [s]', 'Fontsize', 24)
ylabel('Mass Flow Rate [kg/s]', 'Fontsize', 24)

figure(2)
subplot(2,2,3)

plot(data(:,1), data(:,3),'-',data(:,1),data(:,12),'--',data(:,1),data(:,11),':','Linewidth',2)
title('Combustion Chamber Volume Parameters', 'Fontsize', 24)
legend('Fuel Regression Radius', 'Burn Area', 'Combustion Chamber Volume')
xlabel('Time [s]', 'Fontsize', 24)
ylabel('Radius [m] | Area [m^2] |  Volume [m^3]', 'Fontsize', 24)

figure(2)
subplot(2,2,4)

plot(data(:,1), data(:,8),'-','Linewidth',2)
title('Thrust', 'Fontsize', 24)
legend('Thrust')
xlabel('Time [s]', 'Fontsize', 24)
ylabel('Thrust [N]', 'Fontsize', 24)