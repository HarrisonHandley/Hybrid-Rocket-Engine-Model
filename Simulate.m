function data = Simulate(Lcomb, Lpre, Lpost, Din, Rg, Psup,...
    Tsup, Patm, Tatm, Tcomb, Aox, Cdox, rho_fuel, y, A_thr,...
    exit_theta, ratio, time, t, r)
%%  Constants
R = 8.3145;                 %[J/Kmol]
g = 9.81;                   %[m/s2]
y_air = 1.4;
%%  Initial Conditions
lamda = 0.5*(1+cosd(exit_theta));
Aexit = A_thr * ratio;
Mox = 0;
Mfuel = 0;
Pcomb = Patm;

%%  Unit Conversion
Psup = Psup*6894.7;
%%  Initial Time Array
data = zeros(time/t,12);

% Psup = zeros(time/t + 1, 1);

for i = 0:t:time
    
    %%  Relation Equations
    rho_sup = Psup/(R*Tsup);
    Aburn = 2*pi()*Rg*Lcomb;
    Vcomb = pi()*(Rg^2)*Lcomb+(Lpost+Lpre)*pi()*(Din/4)^2;
    Mprop = Mox + Mfuel;
    OF_ratio = Mox/Mfuel;
    Pexit = Pcomb/((1+((y-1)/2)*(t*Mprop)^2)^(y/(y-1)));
    Texit = Tcomb/(1+((y-1)/2)*(t*Mprop)^2);
    FVac = lamda*Mprop*Mprop*t*sqrt(y*Rg*Texit)+Aexit*Pexit;
    FSL = FVac - Aexit*Patm;
    
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
        %Mox = 0;
        Psup = 0;
        Aburn = 0;
         Pexit = Patm;
         FVac = lamda*Mprop*Mprop*t*sqrt(y*Rg*Texit)+Aexit*Pexit;
         FSL = FVac - Aexit*Patm;
    end
    data(uint32((i/t)+1),:) = [i Pcomb./6894.7 Rg Mox Mfuel Mprop r FSL Pexit Texit Vcomb Aburn];
end

end