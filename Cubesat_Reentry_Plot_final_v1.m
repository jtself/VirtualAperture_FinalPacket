% Self, Justin Thomas
% Reentry Map for CubeSat Deorbit
% California Polytechnic State University, San Luis Obispo
% August 23, 2023

% DESCRIPTION
% This code should produce a relationship between altitude and velocity
% for a ballistic object upon reentry into Earth's atmosphere.

tic
close all; clear all; clc;

%% GOVERNING EQUATION

% define velocity (t) equation from Anderson Hypersonics and Gas Dynamics...textbook
% original equation for delta v for flight path during reentry
% -1/g * dV/dt = ((W/CdS)^-1)*(rho*v^2)/2

%% Running the function

% Initial conditions and defining variables
tspan = [0 51*60]; % hits the earth around here! Velocity doesn't reach 0 above the surface. 

g = 9.80665; % m/s2

% 1-sigma plot ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Ballistic coefficient = Beta = W/CdS.
W = 1.33; % kg; typical
Cd = 2.22; % typical value
r = sqrt(0.05); % hypot of 5x5 cm triangle (10cm x 10cm 1U CS face)
S = 4*pi*r^2; % treat as sphere for tumbling (worst case)
Beta = W/(Cd*S); % N/m2
h0 = 100e3; % m; starting height
v0 = 7.5e3; % m/s; reentry velocity (LEO)

% Initial conditions
%[T0, P0, rho0] = stdatm_SelfJustin(h0);
state = [h0;v0];

% Options for ode45
options = odeset('RelTol',1e-8, 'AbsTol', 1e-8);

% NOMINAL (MEAN VALUE) for Cd
% Call ode45 function
[time,newstateVect] = ode45(@reentry_velocity, tspan,state,options,Beta);
newstate_km1 = [newstateVect(:,1)./1000, newstateVect(:,2)./1000]; % km and km/s

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 3-SIGMA for Cd
uncertain = 0.2; % per paper
Cd3 = Cd + uncertain*Cd;
Beta3 = W/(Cd3*S); % N/m2
Cd3minus = Cd - uncertain*Cd;
Beta3_minus = W/(Cd3minus*S); % N/m2

% Call ode45 function
% Plus Cd3
[time3,newstateVect3] = ode45(@reentry_velocity, tspan,state, options,Beta3);
newstate_km3 = [newstateVect3(:,1)./1000, newstateVect3(:,2)./1000]; % converts m/s to km/s

% Minus Cd3
[time3_minus,newstateVect3_minus] = ode45(@reentry_velocity, tspan,state,options, Beta3_minus);
newstate_km3_minus = [newstateVect3_minus(:,1)./1000, newstateVect3_minus(:,2)./1000]; % converts m/s to km/s

%% Total Pressure and Density Using Normal Shock Eqns

% Calculate Mach vector at each altitude.
hVect = newstateVect(:,1);

% Extract Ttotal, Ptotal, and rho from (Freestream) from earlier calcs
% Preallocate
T = zeros(length(hVect),1);
P = T;
rho = T;

for i = 1:length(newstateVect)
[T(i), P(i), rho(i)] = stdatm_SelfJustin(hVect(i));
end
T = T';
P = P';
rho = rho';

% Preallocate loop outputs
a = zeros(length(hVect),1);
gammaNew = zeros(length(hVect),1);
M = zeros(length(hVect),1);
Pt = zeros(length(hVect),1);
Pt0 = zeros(length(hVect),1);
Tt0 = zeros(length(hVect),1);

Pt_postshock = zeros(length(hVect),1);
P_postshock = zeros(length(hVect),1);
Tt_postshock = zeros(length(hVect),1);
T_postshock = zeros(length(hVect),1);
M_postshock = zeros(length(hVect),1);

Pt_postshock_altitude = zeros(length(hVect),1);
P_postshock_altitude = zeros(length(hVect),1);
Tt_postshock_altitude = zeros(length(hVect),1);
T_postshock_altitude = zeros(length(hVect),1);
rho_postshock_altitude = zeros(length(hVect),1);
M_postshock_altitude = zeros(length(hVect),1);

for i = 1:length(hVect)
    [a(i),gammaNew(i)] = speed_of_sound(hVect(i),T(i));
    M(i) = newstateVect(i,2) / a(i); % unitless; units check out

    % Using the Mach vector we can now calculate stagnation parameters at each altitude
    Tt0(i,1)= T(i).*(1 + (gammaNew(i)-1)*0.5*M(i).^2); %freestream isentropic relations
    Pt0(i,1)= P(i).*(Tt0(i,1)/T(i)).^(gammaNew(i)./(gammaNew(i)-1)); %freestream isentropic relations
   
    % gamma convergence after normal shock
    gammaPostShock = gammaNew(i);
    for kk=1:500
        [Pt_postshock(kk), P_postshock(kk), Tt_postshock(kk), T_postshock(kk), rho_postshock(kk), M_postshock(kk)] = normalshock(M(i),gammaPostShock,Pt0(i),P(i),Tt0(i),T(i),rho(i));
        
        % Constants
        R = 287; % J/kg*K
        theta = 3056; % K (thermal constant)
        gammaPerf = 1.4; % specific heat ratio for calorically perfect (ideal) gas

        % CALORICALLY IMPERFECT GAS
        denom = 1 + (gammaPerf - 1)*( (theta/T_postshock(kk)).^2 * ((exp(theta/T_postshock(kk))) / (exp(theta/T_postshock(kk)) - 1).^2  )  );
        gammaPostShock = 1 + (gammaPerf - 1) / denom; 

    end
    Pt_postshock_altitude(i,1) = Pt_postshock(kk);
    P_postshock_altitude(i,1) = P_postshock(kk);
    Tt_postshock_altitude(i,1) = Tt_postshock(kk);
    T_postshock_altitude(i,1) = T_postshock(kk);
    rho_postshock_altitude(i,1) = rho_postshock(kk);
    M_postshock_altitude(i,1) = M_postshock(kk);

    %[Pt_postshock(i), P_postshock(i), Tt_postshock(i), T_postshock(i), rho_postshock(i), M_postshock(i)] = [Pt_postshock(kk), P_postshock(kk), Tt_postshock(kk), T_postshock(kk), rho_postshock(kk), M_postshock(kk)];
end

% Plot results
%{
figure()
plot(a,hVect)

% Graph pretty 
ylim padded 
xlim tight 
xLab = xlabel('a','Interpreter','latex'); 
yLab = ylabel('alt','Interpreter','latex'); 
plotTitle = title('Speed of sound by altitude','interpreter','latex'); 
set(plotTitle,'FontSize',14,'FontWeight','bold') 
set(gca,'FontName','Palatino Linotype') 
set([xLab, yLab],'FontName','Palatino Linotype') 
set(gca,'FontSize', 9) 
set([xLab, yLab],'FontSize', 14) 
grid on 

figure()
plot(M,hVect)

% Graph pretty 
ylim padded 
xlim tight 
xLab = xlabel('Mach','Interpreter','latex'); 
yLab = ylabel('Alt','Interpreter','latex'); 
plotTitle = title('Mach by altitude','interpreter','latex'); 
set(plotTitle,'FontSize',14,'FontWeight','bold') 
set(gca,'FontName','Palatino Linotype') 
set([xLab, yLab],'FontName','Palatino Linotype') 
set(gca,'FontSize', 9) 
set([xLab, yLab],'FontSize', 14) 
grid on 


figure()
plot(Pt_postshock_altitude(1:238),hVect(1:238)) 

% Graph pretty 
ylim padded 
xlim tight 
xLab = xlabel('Stagnation Pressure','Interpreter','latex'); 
yLab = ylabel('altitude','Interpreter','latex'); 
plotTitle = title('Stagnation Pressure by Altitude','interpreter','latex'); 
set(plotTitle,'FontSize',14,'FontWeight','bold') 
set(gca,'FontName','Palatino Linotype') 
set([xLab, yLab],'FontName','Palatino Linotype') 
set(gca,'FontSize', 9) 
set([xLab, yLab],'FontSize', 14) 
grid on 

figure()
plot(Tt_postshock_altitude(1:238),hVect(1:238))

% Graph pretty 
ylim padded 
xlim tight 
xLab = xlabel('Stagnation Temperature','Interpreter','latex'); 
yLab = ylabel('Altitude','Interpreter','latex'); 
plotTitle = title('Total Temperature by Altitude','interpreter','latex'); 
set(plotTitle,'FontSize',14,'FontWeight','bold') 
set(gca,'FontName','Palatino Linotype') 
set([xLab, yLab],'FontName','Palatino Linotype') 
set(gca,'FontSize', 9) 
set([xLab, yLab],'FontSize', 14) 
grid on 

figure()
%yyaxis right
%plot(rho_postshock_altitude(1:238),hVect(1:238))
%hold on
%yyaxis left
%plot(rho,hVect)

t = tiledlayout(1,1);
ax1 = axes(t);
plot(ax1,rho_postshock_altitude(1:238),hVect(1:238),'-r')
ax1.XColor = 'r';
ax1.YColor = 'r';

ax2 = axes(t);
plot(ax2,rho,hVect,'-k')
ax2.XAxisLocation = 'top';
ax2.YAxisLocation = 'right';
ax2.Color = 'none';
ax1.Box = 'off';
ax2.Box = 'off';

% Graph pretty 
ylim padded 
xlim tight 
xLab = xlabel('rho','Interpreter','latex'); 
yLab = ylabel('altitude','Interpreter','latex'); 
plotTitle = title('Rho before and after normal shock','interpreter','latex'); 
set(plotTitle,'FontSize',14,'FontWeight','bold') 
set(gca,'FontName','Palatino Linotype') 
set([xLab, yLab],'FontName','Palatino Linotype') 
set(gca,'FontSize', 9) 
set([xLab, yLab],'FontSize', 14) 
grid on 
legend('rho after shock','rho freestream', 'interpreter','latex','Location', 'best')
%}

%% Fig 1

% Figure 1: Combination of velocity-altitude map with Mach altitude (two
% x-axes)
for fig = 1

% FIGURE 1
figure()
t = tiledlayout(1,1);
ax1 = axes(t);
% Plot mean 
plot(ax1,newstate_km1(:,2), newstate_km1(:,1),'--')
ax1.XColor = 'r';
ax1.YColor = 'r';
hold on

% Interpolation for colorfill plot (shade between two curves)
% since ode45 outputs different x-vectors each time!!
for colorfill = 1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% SAVE THIS CODE
% Bottom
x1 = newstate_km3(:,2)';
y1 = newstate_km3(:,1)';

% Top
% x2 = newstate_km3_minus(1:857,2)';
% y2 = newstate_km3_minus(1:857,1)';

x2 = newstate_km3_minus(:,2)';
y2 = newstate_km3_minus(:,1)';

% Plot shaded area between two
xm = [x1(1:length(x2)); x2];                                                                  % Create 1X: Matrix
xl = [min(xm(:))  max(xm(:))];
xu = linspace(xl(1), xl(2), 1000);                                                  % Create Sorted Interpolation Vector
ym = [y1(1:length(y2)); y2]; 
yi = zeros(2,1000);

% Create ‘Y’ Matrix
    for k1 = 1:size(xm,1)
        yi(k1,:) = interp1(xm(k1,:), ym(k1,:), xu, 'linear','extrap');                  % Interpolate To Provide Common (x,y) Values
    end
ymin = min(yi);
ymax = max(yi);

% Color fill
patch([xu(:)' fliplr(xu(:)')], [ymin fliplr(ymax)], 'r','FaceAlpha',0.3) % Plot Filled Background

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% SAVE THIS CODE

end % colorfill

% % Plot upper and lower bounds
% plot(ax1,newstate_km3(:,2), newstate_km3(:,1),'--b')
% plot(ax1,newstate_km3_minus(:,2), newstate_km3_minus(:,1),'--b')

% Graph pretty for Colorfill
ylim padded 
xlim tight 
xLab = xlabel('Velocity [km/s]','Interpreter','latex'); 
yLab = ylabel('Altitude [km]','Interpreter','latex');  
set(gca,'FontName','Palatino Linotype') 
set([xLab, yLab],'FontName','Palatino Linotype') 
set(gca,'FontSize', 10) 
set([xLab, yLab],'FontSize', 12) 
grid on 
legend('Mean $C_d$','$3\sigma$ $C_d$','interpreter','latex','Location', 'south')

% Plot Mach against altitude
ax2 = axes(t);
plot(ax2,M,hVect./1000,'-k','Linewidth',2)  % PRE SHOCK, NOT POST SHOCK
ax2.XAxisLocation = 'top';
ax2.YAxisLocation = 'right';
ax2.Color = 'none';
ax1.Box = 'off';
ax2.Box = 'off';

% Graph pretty 
ylim padded 
xlim tight 
xLab = xlabel('$M_\infty$','Interpreter','latex'); 
yLab = ylabel('Altitude [km]','Interpreter','latex'); 
%plotTitle = title('Velocity and Mach by Altitude','interpreter','latex'); 
%set(plotTitle,'FontSize',16,'FontWeight','bold') 
set(gca,'FontName','Palatino Linotype') 
set([xLab, yLab],'FontName','Palatino Linotype') 
set(gca,'FontSize', 10) 
set([xLab, yLab],'FontSize', 12) 
grid off
legend('$M_\infty$', 'interpreter','latex','Location', 'southeast')

end % fig 1

%% Fig. 2
for fig = 2
%{
Similarly, show Figure 5 Stagnation Temperature on the bottom axis, 
Figure 4 Stagnation Pressure on the top x-axis; altitude y-axis. 
%}

% No rho, no gamma, etc. 

% Fig. 2
figure()
%yyaxis right
%plot(rho_postshock_altitude(1:238),hVect(1:238))
%hold on
%yyaxis left
%plot(rho,hVect)

% NOTE: limit vector plot elements to 1:233 since Mach goes <1 at 234.
% *********Normal shock assumption fails @ (234:end); about 70km up. 

t = tiledlayout(1,1);
ax1 = axes(t);
plot(ax1,Tt_postshock_altitude(1:233),hVect(1:233)/1000,'r','Linewidth',2)
ax1.XColor = 'r';
ax1.YColor = 'r';

% Graph pretty 
ylim padded 
xlim tight 
xLab = xlabel('Stagnation Temperature [K]','Interpreter','latex'); 
yLab = ylabel('Altitude [km]','Interpreter','latex'); 
%plotTitle = title('Stagnation Temperature and Pressure by Altitude','interpreter','latex'); 
%set(plotTitle,'FontSize',14,'FontWeight','bold') 
set(gca,'FontName','Palatino Linotype') 
set([xLab, yLab],'FontName','Palatino Linotype') 
set(gca,'FontSize', 10) 
set([xLab, yLab],'FontSize', 12) 
grid on 
hold on
%legend('Stagnation Temperature', 'interpreter','latex','Location', 'northeast')

ax2 = axes(t);
plot(ax2,Pt_postshock_altitude(1:233),hVect(1:233)/1000,'--k','Linewidth',2)
ax2.XAxisLocation = 'top';
ax2.YAxisLocation = 'right';
ax2.Color = 'none';
ax1.Box = 'off';
ax2.Box = 'off';

% Graph pretty 
ylim padded 
xlim tight 
xLab = xlabel('Stagnation Pressure [Pa]','Interpreter','latex'); 
yLab = ylabel('Altitude [km]','Interpreter','latex'); 
%plotTitle = title('Stagnation Temperature and Pressure by Altitude','interpreter','latex'); 
%set(plotTitle,'FontSize',16,'FontWeight','bold') 
set(gca,'FontName','Palatino Linotype') 
set([xLab, yLab],'FontName','Palatino Linotype') 
set(gca,'FontSize', 10) 
set([xLab, yLab],'FontSize', 12) 
grid on 
%legend('Stagnation Pressure', 'interpreter','latex','Location', 'best')
end % fig 2

%% Mag plot 
%{
% From 2022 AMS paper work... Eduardo Cardenas wrote this originally, I
% modified it. 

% Parameters
 R1 = [458.37 305.58 229.18 183.34];   % [m]
 S1 = 89.5*10^3:-5*10^3:0;   % [m]
 
% Equations 
sprime = @(s,R) ((1/s)+(2/R))^-1;
m = @(sp,s) -1*sp/s;

% Data
SP_1 = zeros(length(S1),length(R1));
m_1 = zeros(length(S1),length(R1));
for ii = 1:length(R1)
   Rii = R1(ii);
   for i=1:length(S1)
      si = S1(i);
      SP_1(i,ii) = sprime(si,Rii);
      m_1(i,ii) = m(SP_1(i),si);
   end 
end 

x = 1;
y = 2;

figure()
plot(abs(m_1(:,1)),S1*10^-3,abs(m_1(:,2)),S1*10^-3,abs(m_1(:,3)),S1*10^-3,'LineWidth',2)

% Old code
%{
%legend('R1','R2','R3')
ax = gca;
%ax.XDir = 'reverse';
ax.XMinorTick = 'on';
ax.YMinorTick = 'on';
ylabel('Distance [km]')
xlabel('Magnification')
grid on 
grid minor
%} 

% Graph pretty 
ylim padded 
xlim tight 
xLab = xlabel('Magnification','Interpreter','latex'); 
yLab = ylabel('Distance [km]','Interpreter','latex'); 
%plotTitle = title('Optical Magnification','interpreter','latex'); 
%set(plotTitle,'FontSize',16,'FontWeight','bold') 
set(gca,'FontName','Palatino Linotype') 
set([xLab, yLab],'FontName','Palatino Linotype') 
set(gca,'FontSize', 10) 
set([xLab, yLab],'FontSize', 12) 
grid on 

%}
toc

%% Functions
function dstate = reentry_velocity(time, state,ballistic_coefficient)

% Self, Justin
% Reentry Velocity Function using numerical integration

% Description
% This function uses the velocity equation (anderson) for ballistic reentry and
% outputs a state vector in which the COLUMNS ARE: 
    % (1) height above earth's surface
    % (2) velocity (m/s)

% We obtain density at every time step using standard atmosphere model code developed in A215 and optimized in AERO302. 

alt = state(1); % height part of the state vector

% Changing gravity with altitude
g = gravity(alt);

% ~~~~~~~~~~ Replace with standard atmosphere solver of your choice ~~~~
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
[~,~,rho] = stdatm_SelfJustin(alt); % get rho at height
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
dx = -state(2); % velocity is the second column of the state vector

ddx = -(rho * g) * 0.5* (ballistic_coefficient^(-1)) * (dx^2); % anderson 
dstate = [dx; ddx]; % final state is [altitude; velocity] 

end

function [g] = gravity(h)

% https://en.wikipedia.org/wiki/Gravity_of_Earth

g0 = 9.80665; % m/s2
re = 6378*1000; % km

g = g0*((re/(re+h))^2); % m/s

end

