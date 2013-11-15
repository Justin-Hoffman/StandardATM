% Standard Atmosphere
    % 
    % Usage: 
    %    air = standardatm(h,convert)
    %
    % Structureair is a data structure with the following fields:
    %    air.T     - Temperature
    %    air.P     - Pressure
    %    air.rho   - Density
    %    air.a     - Speed of Sound
    %    air.h     - Altitude
    %    air.mu    - Dynamic Viscosity
    %
    %    air.Tsl   - Sea Level Temperature
    %    air.Psl   - Sea Level Pressure
    %    air.rhosl - Sea Level Density
    %
    %    air.theta - Ratio of T to Tsl
    %    air.delta - Ratio of P to Psl
    %    air.sigma - Ratio of rho to rhosl
    %
    %    air.Ta     - Temperature Lapse Rate
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
    hb1 = 36089; % Top of Troposhere
    hb2 = 65617; % Top of lapse rate 1 Stratosphere
    hb3 = 104987; % Top of lapse rate 2 Stratosphere
    hb4 = 154199; % Top of lapse rate 3 Stratosphere
    hb5 = 167323; % Top of lapse rate 1 Mesosphere
    hb6 = 232940;
    if nargin > 1
       h = h*3.28084;
    end
    
    % Error Checking
    if min(h) > hb6
        exception = MException('VerifyOutput:OutOfBounds',...
            'Altitude outside the allowable limits');
        throw(exception);
    end
    if max(h) > hb6
       if min(h) <= hb6
            display('WARNING: Max altitude out of range, truncated results returned');
       else 
            exception = MException('VerifyOutput:OutOfBounds',...
             'Altitude outside the allowable limits');
            throw(exception);
       end  
    end
    
    g = 32.2;
    R = 1716;
    Tsl = 518.67;
    Psl = 14.697*144;
    rhosl = 23.77*10^-4;
    gma = 1.4;
    muinterp = [3.8360    3.7370    3.6370    3.5340    3.4300    3.3240    3.2170    3.1070    2.9950    2.9690    2.9690    2.9690    2.9690    2.9840 3.0180    3.0520    3.0870    3.5110    3.2790    2.8460]*10^-7;
    hinterp = [ -5000           0        5000       10000       15000       20000       25000       30000       35000       40000       45000       50000 60000       70000       80000       90000      100000      150000      200000      250000];
    
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
            [T(i) P(i) rho(i) a(i)] = mesosphere1(h(i));
        elseif h(i) <= hb6
            [T(i) P(i) rho(i) a(i)] = mesosphere2(h(i));
        end
        mu(i) = interp1(hinterp,muinterp,h(i));
    end

    
    %Assemble Structure
    air.T = T;
    air.Tsl = Tsl;
    air.P = P;
    air.Psl = Psl;
    air.rho = rho;
    air.rhosl = rhosl;
    air.mu = mu;
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
       air.mu = air.mu*4.448222/(0.092903);
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
       T = TS+a*(h-hb2);
       P = PS*(T/TS)^(-g/(a*R));
       rho = rhoS.*(T/TS)^(-1*(1+g/(a*R)));
    end
    function [T P rho a] = stratosphere3(h)
       a = 2.8/1000/(5/9*3.28084);
       [TS PS rhoS] = stratosphere2(hb3);
       T = TS+a*(h-hb3);
       P = PS*(T/TS)^(-g/(a*R));
       rho = rhoS.*(T/TS)^(-1*(1+g/(a*R)));
    end
    function [T P rho a] = mesosphere1(h)
       a = 0;
       [TT PT rhoT] = stratosphere3(hb4);
       T = TT;
       P = PT*exp(-g*(h-hb4)/(R*T));
       rho = rhoT*exp(-g*(h-hb4)/(R*T));
    end
    function [T P rho a] = mesosphere2(h)
       a = -0.00085344;
       [TS PS rhoS] = mesosphere1(hb5);
       T = TS+a*(h-hb5);
       P = PS*(T/TS)^(-g/(a*R));
       rho = rhoS.*(T/TS)^(-1*(1+g/(a*R)));
    end
end
