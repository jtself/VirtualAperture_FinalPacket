function [T, P, rho] = stdatm_SelfJustin(h)
%stdatm_SelfJustin ARDC Std Atmosphere
% %   Write a MATLAB function that calculates the air properties at a given 
% altitude using the 1959 ARDC standard atmosphere model.  

%{
INPUT
h in meters

OUTPUTS
T in K
P in Pa
rho in kg/m3
%}

% Constants
g0 = 9.80665; % m/s^2
R = 287; % Gas constant dry air, J/(kg*K)

% Set up standard atmosphere parameters
h0 = 0; % Sea level altitude, m 
T0 = 288.16; % Sea level temperature, K
P0 = 101325; % Sea level pressure Pa
rho0 = 1.2250; % Sea level density kg/m^3 
 

% Layer 1, Gradient layer
a1 = -0.0065; % K/m 
h1 = 11000; % Top of Layer 1, m
T1 = T0 + a1*(h1-h0);   % ** make sure this is h1, T1, etc **
P1 = P0* (T1/T0)^(-g0/(a1*R));
rho1 = rho0 * (T1/T0)^-((g0/(a1*R))+1);

% Layer 2, Isothermal layer (no lapse rate, for iso layers)
a2 = 0;
h2 = 25000; % Top of Layer 2, m
T2 = T1;
P2 = P1 * exp(-g0/(R*T1)*(h2-h1));
rho2 = rho1 * exp(-g0/(R*T1)*(h2-h1));

% Layer 3, Gradient layer  
a3 = 0.0030; % K/m 
h3 = 47000; % Top of Layer 3, m
T3 = T2 + a3*(h3-h2);  
P3 = P2* (T3/T2)^(-g0/(a3*R));
rho3 = rho2 * (T3/T2)^-((g0/(a3*R))+1);

% Layer 4, Isothermal layer
a4 = 0;         % 47 km to 53 km
h4 = 53000;     % top of layer 4, m
T4 = T3;
P4 = P3 * exp(-g0/(R*T3)*(h4-h3));
rho4 = rho3 * exp(-g0/(R*T3)*(h4-h3));

%  Layer 5, Gradient layer  
a5 = -0.0045;   % K/m up to 79 km        
h5 = 79000;     % top of layer 5, m
T5 = T4 + a5*(h5-h4);   
P5 = P4* (T5/T4)^(-g0/(a5*R));
rho5 = rho4 * (T5/T4)^-((g0/(a5*R))+1); 

%  Layer 6, Isothermal layer            
a6 = 0;         % 79 km to 90 km
h6 = 90000;     % top of layer 6, m
T6 = T5;
P6 = P5 * exp(-g0/(R*T5)*(h6-h5));
rho6 = rho5 * exp(-g0/(R*T5)*(h6-h5));


% Layer 7, Gradient layer  
a7 = 0.0040; % K/m up to 100 km
h7 = 100000;
T7 = T6 + a7*(h7-h6);   
P7 = P6 * (T7/T6)^(-g0/(a7*R));
rho7 = rho6 * (T7/T6)^-((g0/(a7*R))+1); 


% based on altitude h, solve for P, rho, and T at that altitude
if h < h1
    % Layer 1: gradient layer 1
    % DONE
    T = T0 + a1*(h-h0);
    P = P0 * (T/T0)^(-g0/(a1*R));
    rho = rho0 * (T/T0)^-((g0/(a1*R))+1);

elseif h < h2
    % Layer 2: isothermal layer 1
    % DONE
    T = T1;
    P = P1 * exp(-g0/(R*T)*(h-h1));
    rho = rho1 * exp((-g0/(R*T))*(h-h1));

elseif h < h3
    % Layer 3: gradient layer 2
    % DONE
    T = T2 + a3*(h-h2);
    P = P2 * (T/T2)^(-g0/(a3*R));
    rho = rho2 * (T/T2)^-((g0/(a3*R))+1);
    
elseif h < h4 
    % Layer 4: isothermal layer 2
    % DONE
    T = T4;
    P = P4 * exp(-g0/(R*T)*(h-h4));
    rho = rho4 * exp((-g0/(R*T))*(h-h4));

elseif h < h5
    % Layer gradient layer 3
    T = T4 + a5*(h-h4);                
    P = P4 * (T/T4)^(-g0/(a5*R));
    rho = rho4 * (T/T4)^-((g0/(a5*R))+1);

elseif h < h6
    % isothermal layer 3 
    T = T6;
    P = P6 * exp(-g0/(R*T)*(h-h6));
    rho = rho6 * exp((-g0/(R*T))*(h-h6));
   

elseif h <= h7
    % gradient layer 4
    T = T6 + a7*(h-h6);                
    P = P6 * (T/T6)^(-g0/(a7*R));
    rho = rho6 * (T/T6)^-((g0/(a7*R))+1);

    end

end

