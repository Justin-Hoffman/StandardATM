% Standard Atmosphere Properties
    % 
    % Usage: 
    %    air = standardatm(h,convert)
    %
    % Structureair is a data structure witht the following fields:
    %    air.T     - Temperature
    %    air.P     - Pressure
    %    air.rho   - Density
    %    air.h     - Altitude
    %
    %    air.Tsl   - Sea Level Temperature
    %    air.Psl   - Sea Level Pressure
    %    air.rhosl - Sea Level Density
    %
    %    air.theta - Ratio of T to Tsl
    %    air.delta - Ratio of P to Psl
    %    air.sigma - Ratio of rho to rhosl
    %
    %    air.a     - Temperature Lapse Rate
    %
    %
    % Units: 
    %    Units are Ft, Slugs, Lbf, Rankine by default. An optional second
    %    argument causes the units returned to be in M, Kg, N, Kelvin
    %
    %
function air = standardatm(h,varargin)
    %Units are in US units, and converted before function return if
    %nargin>1
    hb1 = 36100; % Top of Troposhere
    hb2 = 20*3280.84; % Top of lapse rate 1 Stratosphere
    hb3 = 32*3280.84; % Top of lapse rate 2 Stratosphere
    hb4 = 47*3280.84; % Top of lapse rate 3 Stratosphere
    hb5 = 51*3280.84; % Top of lapse rate 1 Mesosphere
    if nargin > 1
       h = h*3.28084;
    end
    
    g = 32.2;
    R = 1716;
    Tsl = 518.67;
    Psl = 14.697*144;
    rhosl = 23.77*10^-4;
    gma = 1.4;
    
    for i = 1:length(h)
        if h(i)<hb1
            [T(i) P(i) rho(i) a(i)] = troposphere(h(i));
        elseif h(i) <= hb2
            [T(i) P(i) rho(i) a(i)] = stratosphere1(h(i));
        elseif h(i) <= hb3
            [T(i) P(i) rho(i) a(i)] = stratosphere2(h(i));
        elseif h(i) <= hb4
            [T(i) P(i) rho(i) a(i)] = stratosphere3(h(i));
        elseif h(i) <= hb5
            [T(i) P(i) rho(i) a(i)] = mesosphere(h(i));
        end
    end
    %Assemble Structure
    air.T = T;
    air.Tsl = Tsl;
    air.P = P;
    air.Psl = Psl;
    air.rho = rho;
    air.rhosl = rhosl;
    air.a = sqrt(gma*R*air.T);
    air.sigma = air.rho/rhosl;
    air.delta = air.P/Psl;
    air.theta = air.T/Tsl;
    air.Ta = a;
    air.R = R;
    air.h = h;
    
    %Convert Units if Nargin > 1
    if nargin > 1
       air.T = air.T*5/9;
       air.P = air.P*47.8803;
       air.Tsl = air.Tsl*5/9;
       air.Psl = air.Psl*47.8803;
       air.a = air.a/3.28084;
       air.rho = air.rho*515.378818;
       air.rhosl = air.rhosl*515.378818;
       air.Ta = air.Ta*5/9*3.28084;
       air.R = 286.9;
       air.h = air.h/3.28084;
    end
    
    %Standard Atmosphere Equations
    function [T P rho a] = troposphere(h)
       a = -.00356;
       T = Tsl+a*h;
       P = Psl*(T/Tsl)^(-g/(a*R));
       rho = rhosl.*(T/Tsl)^(-1*(1+g/(a*R)));
    end

    function [T P rho a] = stratosphere1(h)
       a = 0;
       [TT PT rhoT] = troposphere(hb1);
       T = TT;
       P = PT*exp(-g*(h-hb1)/(R*T));
       rho = rhoT*exp(-g*(h-hb1)/(R*T));
    end
    function [T P rho a] = stratosphere2(h)
       a = 1.0/1000/(5/9*3.28084);
       [TS PS rhoS] = stratosphere1(hb2);
       T = TS+a*h;
       P = PS*(T/TS)^(-g/(a*R));
       rho = rhoS.*(T/TS)^(-1*(1+g/(a*R)));
    end
    function [T P rho a] = stratosphere3(h)
       a = 1.0/1000/(5/9*3.28084);
       [TS PS rhoS] = stratosphere2(hb3);
       T = TS+a*h;
       P = PS*(T/TS)^(-g/(a*R));
       rho = rhoS.*(T/TS)^(-1*(1+g/(a*R)));
    end
    function [T P rho a] = mesosphere(h)
       a = 0;
       [TT PT rhoT] = stratosphere3(hb4);
       T = TT;
       P = PT*exp(-g*(h-hb1)/(R*T));
       rho = rhoT*exp(-g*(h-hb1)/(R*T));
    end
end