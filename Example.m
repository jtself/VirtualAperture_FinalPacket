%{

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Justin Self, Cal Poly, Winter 2024: AERO 407 | Reentry Aerodynamics
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Deliverable #1 ... Problem #2
%}

% housekeeping
clear all; close all; clc;

%..... ode45 options
options = odeset('RelTol', 1e-8, 'AbsTol',1e-8); 

disp("Deliverable #2: Reentry Mechanics")
disp(" ")

% Define knowns
re = 6378; % km; Radius of earth
tspan = [0 2e3]; % sec
alt0 = (6478 - re)*1000;
v0 = 11e3; % km/s
m = 1e4; % kg
S = pi*(3.9/2)^2; % m2
Cd = 1.55; % not 2.22 for spaceflight; this is from 
beta0 = m/(Cd*S);
state = [alt0;v0];

for i = 0:9
    beta = beta0 - (i/10)*beta0;
    [timenew,statenew] = ode45(@reentry_velocity,tspan,state,options,beta);
    time{1,i+1} = timenew;
    stateNEW{1,i+1} = statenew;
end

% Extract data
state1 = stateNEW{1,1}./1000;
state2 = stateNEW{1,2}./1000;
state3 = stateNEW{1,3}./1000;
state4 = stateNEW{1,4}./1000;
state5 = stateNEW{1,5}./1000;
state6 = stateNEW{1,6}./1000;
state7 = stateNEW{1,7}./1000;
state8 = stateNEW{1,8}./1000;
state9 = stateNEW{1,9}./1000;
state10 = stateNEW{1,10}./1000;

figure()
plot(state1(:,2),state1(:,1),'LineWidth',2) % nominal Beta
hold on
plot(state2(:,2),state2(:,1),'LineWidth',1) % 90% of Beta
plot(state3(:,2),state3(:,1),'LineWidth',1) % 80% 
plot(state4(:,2),state4(:,1),'LineWidth',1) % 70%
plot(state5(:,2),state5(:,1),'LineWidth',1) % 60%
plot(state6(:,2),state6(:,1),'LineWidth',1) % 50%
plot(state7(:,2),state7(:,1),'LineWidth',1) % 40%
plot(state8(:,2),state8(:,1),'LineWidth',1) % 30%
plot(state9(:,2),state9(:,1),'LineWidth',1) % 20%
plot(state10(:,2),state10(:,1),'LineWidth',1) % 10%
yline(0) % show surface of earth

% Graph pretty 
ylim padded 
xlim tight 
xLab = xlabel('Velocity [km/s]','Interpreter','latex'); 
yLab = ylabel('Altitude [km]','Interpreter','latex'); 
plotTitle = title('Velocity altitude map for changing ballistic coefficient','interpreter','latex'); 
set(plotTitle,'FontSize',14,'FontWeight','bold') 
set(gca,'FontName','Palatino Linotype') 
set([xLab, yLab],'FontName','Palatino Linotype') 
set(gca,'FontSize', 9) 
set([xLab, yLab],'FontSize', 12) 
grid on 
legend('Nominal $\beta$','$0.9 \beta$','$0.8 \beta$','$0.7 \beta$',...
    '$0.6 \beta$','$0.5 \beta$','$0.4 \beta$','$0.3 \beta$','$0.2 \beta$','$0.1 \beta$','interpreter','latex','Location', 'best')

%% Flight path angle changes
gamma0 = [-90,-75,-60,-45,-30,-15,-1];
distx = 0;
h = alt0;
alty = alt0;
z = 0; % not using this
heatshield = "no";

dragPercentChange = 0; % not needed here since heatshield = 'no'

for i = 1:7
    fpa0 = gamma0(i);
    vx0 = v0*cosd(fpa0);
    vy0 = v0*sind(fpa0);
    state = [distx;alty;v0;gamma0(i);h];
    [~,statenew] = ode45(@reentry_velocity_gammaVh,tspan,state,options,beta,m,z,dragPercentChange,heatshield);
    stateNEW{1,i} = statenew;
end


% Extract data (change to km and km/s)
fpa90 = stateNEW{1,1}./1000;
fpa75 = stateNEW{1,2}./1000;
fpa60 = stateNEW{1,3}./1000;
fpa45 = stateNEW{1,4}./1000;
fpa30 = stateNEW{1,5}./1000;
fpa15 = stateNEW{1,6}./1000;
fpa1 = stateNEW{1,7}./1000; % literally orbit (not realistic)

%% Plot flight path angle against y-velocity

%{
indices:
1. dist x (km)
2. dist y (alt), km
3. velocity (mag), km/s
4. gamma
5. height (altitude, y)
%}

figure() % plot y(veloc, alt)
plot(fpa90(:,3),fpa90(:,2),'LineWidth',2) % should be just like nominal
hold on
% Plot flight path angle changes
plot(fpa75(:,3),fpa75(:,2),'LineWidth',2)
plot(fpa60(:,3),fpa60(:,2),'LineWidth',2)
plot(fpa45(:,3),fpa45(:,2),'LineWidth',2)
plot(fpa30(:,3),fpa30(:,2),'LineWidth',2)
plot(fpa15(:,3),fpa15(:,2),'LineWidth',2)
plot(fpa1(:,3),fpa1(:,2),'LineWidth',2)
yline(0) % show surface of earth

% Graph pretty 
ylim padded 
xlim tight 
xLab = xlabel('Velocity (magnitude) [km/s]','Interpreter','latex'); 
yLab = ylabel('Altitude [km]','Interpreter','latex'); 
plotTitle = title('Velocity altitude map for changing flight path angle $\in [-90,0]^\circ$','interpreter','latex'); 
set(plotTitle,'FontSize',14,'FontWeight','bold') 
set(gca,'FontName','Palatino Linotype') 
set([xLab, yLab],'FontName','Palatino Linotype') 
set(gca,'FontSize', 9) 
set([xLab, yLab],'FontSize', 12) 
grid on 

legend('$\gamma = -90^\circ$','$\gamma = -75^\circ$','$\gamma = -60^\circ$',...
    '$\gamma = -45^\circ$','$\gamma = -30^\circ$','$\gamma = -15^\circ$','$\gamma = -1^\circ$',...
    'interpreter','latex','Location', 'best')


%% Plot flight path angles vs x-distance (ground footprint)
%{
indices:
1. dist x
2. dist y (alt)
3. velocity (mag)
4. gamma
5. height (altitude, y)
%}

figure()
plot(fpa90(:,1),fpa90(:,2),'LineWidth',2) % should be just like nominal
hold on
% Plot x, alt
plot(fpa75(:,1),fpa75(:,2),'LineWidth',2)
plot(fpa60(:,1),fpa60(:,2),'LineWidth',2)
plot(fpa45(:,1),fpa45(:,2),'LineWidth',2)
plot(fpa30(:,1),fpa30(:,2),'LineWidth',2)
plot(fpa15(:,1),fpa15(:,2),'LineWidth',2)
plot(fpa1(:,1),fpa1(:,2),'LineWidth',2)
yline(0) % show surface of earth

% Graph pretty 
ylim padded 
xlim tight 
xLab = xlabel('Ground footprint (x-distance) [km]','Interpreter','latex'); 
yLab = ylabel('Altitude [km]','Interpreter','latex'); 
plotTitle = title('Ground footprint for changing flight path angle $\in [-90,0]^\circ$','interpreter','latex'); 
set(plotTitle,'FontSize',14,'FontWeight','bold') 
set(gca,'FontName','Palatino Linotype') 
set([xLab, yLab],'FontName','Palatino Linotype') 
set(gca,'FontSize', 9) 
set([xLab, yLab],'FontSize', 12) 
grid on 

legend('$\gamma = -90^\circ$','$\gamma = -75^\circ$','$\gamma = -60^\circ$',...
    '$\gamma = -45^\circ$','$\gamma = -30^\circ$','$\gamma = -15^\circ$','$\gamma = -1^\circ$',...
    'interpreter','latex','Location', 'best')

%% Model ground footprint with heatshield deployment

%{
IDEA: 
Show the impact on ground footprint for a sudden increase in drag coefficient during re-entry created by
opening of an umbrella looking heat shield at a particular altitude for at least three altitudes of your choice, and
for at least three Cd values. Suggestion: Consider 10% to 100% increase in base Cd value keeping constant base
values for W and S

Assume flight path angle is: -1 deg; (realistic)
% CIS LUNAR TRAJECTORY REENTRY FLIGHT PATH ANGLE WAS -23.7 degrees.

PLAN: 
    - if statement inside ode funct that causes Cd to instantaneously
    change from nominal to new value
    - new values are: Cd*1.1; Cd*1.5; Cd*2 (10%, 50%, and 100% increase)
    - Do that again for different altitude values (3 total)
%}

% S/c and reentry parameters (update for this part)
Cd = 1.4;
beta = m/(Cd*S);
gamma = -23.7; % flight path angle; degrees for cis-lunar drop in

% Compute cases
z1 = 80e3; % meters 
z2 = 60e3;
z3 = 40e3;
dragPercentChange = [1.1;1.5;2]; % that is, 1.1*Cd = 10% increase in Cd

distx = 0; % original distance in the x direction
alty = alt0;
vx0 = v0*cosd(gamma); % x-velocity
vy0 = v0*sind(gamma); % y-velocity

% Turn heatshield 'on'
heatshield = "yes";
h = alt0;
state = [distx;alty;v0;gamma;h];
tspan = [0 500]; % to truncate unneeded data at end

% .... Call ode for each configuration

%%%%%%%% CAN BE OPTIMIZED (to do)

% 80 km (z1d1 = 80 km alt, drag coeff (1)
[~,statenew_z1d1] = ode45(@reentry_velocity_gammaVh,tspan,state,options,beta,m,z1,dragPercentChange(1),heatshield);
[~,statenew_z1d2] = ode45(@reentry_velocity_gammaVh,tspan,state,options,beta,m,z1,dragPercentChange(2),heatshield);
[~,statenew_z1d3] = ode45(@reentry_velocity_gammaVh,tspan,state,options,beta,m,z1,dragPercentChange(3),heatshield);

% 60 km
[~,statenew_z2d1] = ode45(@reentry_velocity_gammaVh,tspan,state,options,beta,m,z2,dragPercentChange(1),heatshield);
[~,statenew_z2d2] = ode45(@reentry_velocity_gammaVh,tspan,state,options,beta,m,z2,dragPercentChange(2),heatshield);
[~,statenew_z2d3] = ode45(@reentry_velocity_gammaVh,tspan,state,options,beta,m,z2,dragPercentChange(3),heatshield);

% 40 km
[~,statenew_z3d1] = ode45(@reentry_velocity_gammaVh,tspan,state,options,beta,m,z3,dragPercentChange(1),heatshield);
[~,statenew_z3d2] = ode45(@reentry_velocity_gammaVh,tspan,state,options,beta,m,z3,dragPercentChange(2),heatshield);
[~,statenew_z3d3] = ode45(@reentry_velocity_gammaVh,tspan,state,options,beta,m,z3,dragPercentChange(3),heatshield);

% How fast was the capsule going at impact? (Just for fun)
vf.z1d1 = statenew_z1d1(end,3); % m/s
vf.z1d2 = statenew_z1d2(end,3); % m/s
vf.z1d3 = statenew_z1d3(end,3); % m/s
vf.z2d1 = statenew_z2d1(end,3); % m/s
vf.z2d2 = statenew_z2d2(end,3); % m/s
vf.z2d3 = statenew_z2d3(end,3); % m/s
vf.z3d1 = statenew_z3d1(end,3); % m/s
vf.z3d2 = statenew_z3d2(end,3); % m/s
vf.z3d3 = statenew_z3d3(end,3); % m/s

v_impact = [vf.z1d1;vf.z1d2;vf.z1d3];
figure
x = ["10% increase","50% increase","100% increase"];
bar(x,v_impact)
grid on

% Plot ground footprint (x) vs altitude for 80km case
figure()
% 80 km, all three drag coeffs.
plot(statenew_z1d1(:,1)/1000,statenew_z1d1(:,2)/1000, "LineWidth",2) % x-footprint vs height
hold on
plot(statenew_z1d2(:,1)/1000,statenew_z1d2(:,2)/1000, "LineWidth",2)
plot(statenew_z1d3(:,1)/1000,statenew_z1d3(:,2)/1000, "LineWidth",2)

% 60 km
plot(statenew_z2d1(:,1)/1000,statenew_z2d1(:,2)/1000, "LineWidth",2)
plot(statenew_z2d2(:,1)/1000,statenew_z2d2(:,2)/1000, "LineWidth",2)
plot(statenew_z2d3(:,1)/1000,statenew_z2d3(:,2)/1000, "LineWidth",2)

% 40 km
plot(statenew_z3d1(:,1)/1000,statenew_z3d1(:,2)/1000, "LineWidth",2)
plot(statenew_z3d2(:,1)/1000,statenew_z3d2(:,2)/1000, "LineWidth",2)
plot(statenew_z3d3(:,1)/1000,statenew_z3d3(:,2)/1000, "LineWidth",2)

yline(0) % show surface of earth

% ************** NOTE: THEY LOOK OVERLAID, BUT THEY ARE DIFFERENT WHEN YOU ZOOM IN. 

% Graph pretty 
ylim padded 
xlim tight 
xLab = xlabel('Ground footprint (x-distance) [km]','Interpreter','latex'); 
yLab = ylabel('Altitude [km]','Interpreter','latex'); 
plotTitle = title('Ground footprint changes by heatshield deployment altitude','interpreter','latex'); 
set(plotTitle,'FontSize',14,'FontWeight','bold') 
set(gca,'FontName','Palatino Linotype') 
set([xLab, yLab],'FontName','Palatino Linotype') 
set(gca,'FontSize', 9) 
set([xLab, yLab],'FontSize', 12) 
grid on 
legend('80 km, $C_{D1}$','80 km, $C_{D2}$','80 km, $C_{D3}$','60 km, $C_{D1}$','60 km, $C_{D2}$',...
    '60 km, $C_{D3}$','40 km, $C_{D1}$','40 km, $C_{D2}$','40 km, $C_{D3}$', 'interpreter','latex','Location', 'best')

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

[~,~,rho] = stdatm_SelfJustin(alt); % get rho at height
dx = -state(2); % velocity is the second column of the state vector

ddx = -(rho * g) * 0.5* (ballistic_coefficient^(-1)) * (dx^2); % anderson 
dstate = [dx; ddx]; % final state is [altitude; velocity] 

end

function dstate = reentry_velocity_gammaVh(time,state,ballistic_coefficient,m,z,dragPercentChange,heatshield)

% Self, Justin
% Reentry Velocity Function using numerical integration using heatshield
% that pops out (changes CD)

% Description
% This function uses the velocity equation (anderson) for ballistic reentry and
% outputs a state vector in which the COLUMNS ARE: 
    % (1) height above earth's surface
    % (2) velocity (m/s)
    % GAMMA IN DEGREES
    % dragPercentChange is the desired change in CD value (inside ball.
    % coeff beta) for each heat shield configuration. 
    % i.e. dragPercentChange = 1.1 for a new Beta including CD(new) =
    % 1.1*Cd(old); increase of 10%

distx = state(1); % not used
alt = state(2);
v = state(3);

%vx = state(3);
%vy = state(4);
gamma = state(4);
h = state(5); 
% Changing gravity with altitude
g = gravity(alt);

% get rho at height
[~,~,rho] = stdatm_SelfJustin(alt); 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Check if user wants "heatshield" part of problem
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Take user input for heatshield deployment configuration or not
if strcmpi(heatshield,'yes')

    % Pop out heat shield at prescribed altitude
    if alt < z % until here, Cd = Cd normal
        ballistic_coefficient = ballistic_coefficient * dragPercentChange^-1;
    end
end

L = 0; % no lift
R = h + 6378e3; % radius of flight from center of earth (m)
vc = sqrt(g*R);

% Compute drag force
D = 0.5*rho*v^2*m/ballistic_coefficient;
W = m*g;

% Hankey equations 2.119 - 2.121
vdot = -g*(D/W + sind(gamma)); % change in x veloc
dgamma = g/v * (L/W - (1 - v^2/vc^2)*cosd(gamma) ); % change in gamma
dh = v*sind(gamma); % change in height
vx = cosd(gamma)*v;
vy = sind(gamma)*v;
dstate = [vx; vy; vdot; dgamma; dh]; % final state is [altitude; velocity]
end
