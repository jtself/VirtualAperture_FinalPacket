% Groundtrack testing for 1U CubeSat (quick estimate)
% Justin Self
% July 6 2023

clear all; close all; clc;

% Goal: Plot groundtrack of orbital flightpath using real satellite data
% NOTE: Groundtrack function courtesy of Tamas Kis
    % Copyright © 2021 Tamas Kis
    % Last Update: 2022-07-06
    % Website: https://tamaskis.github.io
    % Contact: tamas.a.kis@outlook.comTamas Kis

%{
for PAYLOAD G: AAU CubeSat 
https://en.wikipedia.org/wiki/AAU_CubeSat
Sat ID: 27846

From HEAVENS - ABOVE
Epoch (UTC):	05 July 2023 16:12:59
Eccentricity:	0.0010351
inclination:	98.6854°
perigee height:	809 km
apogee height:	823 km
right ascension of ascending node:	192.7829°
argument of perigee:	51.7474°
revolutions per day:	14.22726258
mean anomaly at epoch:	308.4634°
%}

% ORBITAL PARAMETERS
epoch = [2023,07, 06, 16, 12, 59]; % <-- user enters this
mu = 398600; % km3/s2
revsperday = 14.22726258; % revs
T = revsperday * 86400; % seconds
a = (T*sqrt(mu) / (2*pi)) ^ (2/3); % km
ecc = 0.0010351; 
inc = 98.6854; % deg
alt_peri = 809; % km
alt_apo = 823; % km
RAAN = 192.7829; % deg
w = 51.7474; % deg
Me = 308.4634; % deg

re = 6378; % radius of earth; km

% Orbital radius (magnitude)
r = re + 0.5*(alt_apo + alt_peri); % km; simplifying assumption for circular

% Orbital speed (assume circular)
v1 = sqrt(mu/r); % km/s

% rVect, vVect in Perifocal (orbital plane)
h = sqrt(r*mu*(1 + ecc* cosd(inc)));
r = (h^2/mu) * (1 + ecc * cosd(inc))^-1 *  [cosd(inc); sind(inc); 0];
v = (mu/h) * [-sind(inc); ecc + cosd(inc); 0];

% Convert from orbital plane to ECI
C_ECI_perifocal = peri2geo_313_rotation_matrix(RAAN,inc,w); % in degrees
    
% Transform r, v from Perifocal into ECI
r_eci = C_ECI_perifocal*r;
v_eci = C_ECI_perifocal*v;

rmag_eci = norm(r_eci);
h_eci = cross(r_eci,v_eci);
hmag_eci = norm(h_eci);

% Now we have orbital parameters in ECI.

% Find lla from eci
% UTC of lookup: 05 July 2023 16:12:59
% output of eci2lla = lat, long, alt in geodetic

% Propogate orbit forward using ode45 and orbit equation (no perts)
start = 0;
stop = 10200;
tspan = [start stop];
options = odeset('RelTol', 1e-8, 'AbsTol',1e-8);
state = [r_eci; v_eci]; 

% call ode here (no perts)
[time, statenew] = ode45(@non_impulsive_COAST,tspan,state,options,mu);

% Now we have a vector of r,v.
n = length(statenew);
lla  = zeros(length(statenew),3);
delta_t = 1; % one day

% Build one DAY of time vector (all in UTC)
t1 = datetime(epoch); % start time 
t2 = datetime(epoch(1),epoch(2), epoch(3) + delta_t, epoch(4), epoch(5), epoch(6)); % stop time (plus one day)

t = t1:seconds(86400/n):t2;
DateVector = datevec(t);

for i = 1:n
lla(i,:) = eci2lla([statenew(i,1),statenew(i,2),statenew(i,3)] .* 1e3,DateVector(i,:));
end

lat = lla(:,1);
lon = lla(:,2);

% Call plotting functions
% initialize figure
figure('Position',[300,300,1000,500]);

% Plot groundtrack
ground_track(lat,lon,'Earth')

%% Functions
function [Q] = peri2geo_313_rotation_matrix(RAAN,inc,w)

% This function produces a transformation matrix, Q, from the PERIFOCAL to
% the GEOCENTRIC frames.

% Self, Justin
% Fall 2022

% INPUTS MUST BE IN DEGREES

Cz_RAAN =        [  cosd(RAAN)     sind(RAAN)     0;
                   -sind(RAAN)     cosd(RAAN)     0;
                    0             0             1]; % 3

Cx_inc =        [   1             0             0;
                    0             cosd(inc)      sind(inc);
                    0             -sind(inc)     cosd(inc)]; % 1

Cz_w =        [     cosd(w)        sind(w)        0;
                    -sind(w)       cosd(w)        0;
                    0             0             1]; % 3

Q = Cz_RAAN * Cx_inc * Cz_w;
Q = Q';
end % 313 funct

function [dstate] = non_impulsive_COAST(time,state,mu)
%{ 

*FOR COAST PHASE ONLY*

**This function assumes that thrust is exactly parallel to initial velocity
direction.

This function is for numerical integration (i.e., plug into ode45).
NOTE: only use this function for the BURN PHASE of the non-impulsive
maneuver. For the BURN phase, use "non_impulsive_BURN"

INPUTS: 
    INITIAL CONDITIONS
    state vector(1:3) = position (rx,ry,rz)
    state vector(4:6) = velocity (vx,vy,vz)

OUTPUT: New state vector after integration.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Justin Self, Cal Poly, Fall 2022; Introduction to Orbital Mechanics.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%}

g0 = 9.807; % m/s2
g0 = g0/1000; % to correct the units; km/s2

rx = state(1);
ry = state(2);
rz = state(3);
vx = state(4);
vy = state(5);
vz = state(6);

r_vect = [rx ry rz];
r = norm(r_vect);

v_vect = [vx vy vz];
v = norm(v_vect);

xdot = vx;
ydot = vy;
zdot = vz;
xddot = (-mu*rx)/(r^3);
yddot = (-mu*ry)/(r^3);
zddot = (-mu*rz)/(r^3);

dstate = [xdot; ydot; zdot; xddot; yddot; zddot];

end % coast

%==========================================================================
%
% ground_track  Plots the ground track of an orbit.
%
%   ground_track(lat,lon)
%   ground_track(lat,lon,opts)
%   ground_track(__,planet)
%
% See also planet3D.
%
% Copyright © 2021 Tamas Kis
% Last Update: 2022-07-06
% Website: https://tamaskis.github.io
% Contact: tamas.a.kis@outlook.com
%
% IMAGE SOURCES:
% https://tamaskis.github.io/files/Ground_Track_Image_Sources.pdf
%
%--------------------------------------------------------------------------
%
% ------
% INPUT:
% ------
%   lat     - (1D double array) planetodetic latitude [deg]
%   lon     - (1D double array) planetodetic longitude [deg]
%   opts    - (OPTIONAL) (1×1 struct) plot options
%       • Color     - (char or 1×3 double) line color
%                       --> can be specified as a name, short name, or RGB 
%                           triplet [rgb]
%       • LineWidth - (1×1 double) line width
%       • LineStyle - (char) line style
%   planet  - (OPTIONAL) (char) 'Blank', 'Sun', 'Moon', 'Mercury', 'Venus',
%             'Earth', 'Earth Cloudy', 'Earth Coastlines', 'Earth Night', 
%             'Earth Night Cloudy', 'Mars', 'Jupiter', 'Saturn', 'Uranus',
%             'Neptune', or 'Pluto'
%
%==========================================================================
function ground_track(lat,lon,opts,planet)
    
    % ------------------------------------
    % Sets (or defaults) plotting options.
    % ------------------------------------
    
    % sets line color (defaults to default MATLAB color)
    if (nargin < 3) || isempty(opts) || ~isfield(opts,'Color')
        Color = [0,0.4470,0.7410];
    else
        Color = opts.Color;
    end
    
    % sets line style (defaults to solid line)
    if (nargin < 3) || isempty(opts) || ~isfield(opts,'LineStyle')
        LineStyle = '-';
    else
        LineStyle = opts.LineStyle;
    end
    
    % sets line width (defaults to 1.5)
    if (nargin < 3) || isempty(opts) || ~isfield(opts,'LineWidth')
        LineWidth = 1.5;
    else
        LineWidth = opts.LineWidth;
    end
    
    % sets background (defaults to Earth coastlines)
    if (nargin < 4) || isempty(planet)
        planet = 'Earth Coastlines';
    end
    
    % --------------------------------------
    % Draws background of ground track plot.
    % --------------------------------------
    
    % sets background of ground track plot
    if strcmpi(planet,'Earth Coastlines')
        
        % loads Earth topographic data
        load('topo.mat','topo');
        
        % rearranges Earth topopgrahic data so it goes from -180 to 180 
        % deg longitude
        topoplot = [topo(:,181:360),topo(:,1:180)];
        
        % plots Earth map by making a contour plot of topographic data at
        % elevation of 0
        contour(-180:179,-90:89,topoplot,[0,0],'black');
        
    elseif ~strcmpi(planet,'Blank')
        
        % reads in image
        if strcmpi(planet,'Earth Cloudy')
            cdata = imread('earth.png')+imread('clouds.png');
        elseif strcmpi(planet,'Earth Night')
            cdata = imread('earthnight.png');
        elseif strcmpi(planet,'Earth Night Cloudy')
            cdata = imread('earthnight.png')+0.1*...
                imread('clouds.png');
        else
            cdata = imread(strcat(lower(planet),'.png'));
        end
        
        % plots background
        image('CData',flipud(cdata),'XData',[-180,180],'YData',[-90,90]);
        
        % determines grid style based on background
        if strcmpi(planet,'Sun')
            GridLineWidth = 1;
            GridColor = 'k';
        elseif strcmpi(planet,'Moon')
            GridLineWidth = 0.75;
            GridColor = 'y';
        elseif strcmpi(planet,'Mercury')
            GridLineWidth = 0.75;
            GridColor = 'y';
        elseif strcmpi(planet,'Venus')
            GridLineWidth = 0.5;
            GridColor = 'k';
        elseif strcmpi(planet,'Earth')
            GridLineWidth = 0.5;
            GridColor = 'y';
        elseif strcmpi(planet,'Earth Cloudy')
            GridLineWidth = 1.5;
            GridColor = 'y';
        elseif strcmpi(planet,'Earth Night')
            GridLineWidth = 0.5;
            GridColor = [0.5,0.5,0.5];
        elseif strcmpi(planet,'Earth Night Cloudy')
            GridLineWidth = 0.5;
            GridColor = [0.5,0.5,0.5];
        elseif strcmpi(planet,'Mars')
            GridLineWidth = 0.75;
            GridColor = 'y';
        elseif strcmpi(planet,'Jupiter')
            GridLineWidth = 1;
            GridColor = 'k';
        elseif strcmpi(planet,'Saturn')
            GridLineWidth = 0.65;
            GridColor = 'k';
        elseif strcmpi(planet,'Uranus')
            GridLineWidth = 0.5;
            GridColor = 'k';
        elseif strcmpi(planet,'Neptune')
            GridLineWidth = 0.5;
            GridColor = 'w';
        elseif strcmpi(planet,'Pluto')
            GridLineWidth = 0.65;
            GridColor = 'g';
        end
        
        % manually adds grid
        hold on;
        for i = 1:11
            plot([-180+i*30,-180+i*30],[-90,90],'Color',GridColor,...
                'LineWidth',GridLineWidth,'LineStyle',':');
        end
        for i = 1:5
            plot([-180,180],[-90+i*30,-90+i*30],'color',GridColor,...
                'LineWidth',GridLineWidth,'LineStyle',':');
        end
        
    end
    
    % ----------------------------------------
    % Plotting ground track / axis formatting.
    % ----------------------------------------
    
    % determines indices where ground track crosses figure border (i.e.
    % from 180 to -180 or -180 to 180) to avoid "jumps" in the plot
    j = [];
    for i = 2:length(lon)
        if ((lon(i) > 170) && (lon(i-1) < -170)) || ((lon(i) < -170) &&...
                (lon(i-1) > 170))
            j = [j,i];
        end
    end
    
    % adds last index to "j" in order to plot ground track between last
    % figure border crossing and the last input longitude
    j = [j,length(lon)];
    
    % plots groundtrack (starts new plot every time ground track crosses
    % left border or right border)
    hold on;
    plot(lon(1:(j(1)-1)),lat(1:(j(1)-1)),'Color',Color,'LineStyle',...
        LineStyle,'LineWidth',LineWidth);
    for i = 1:(length(j)-1)
        plot(lon(j(i):(j(i+1)-1)),lat(j(i):(j(i+1)-1)),'Color',Color,...
           'LineStyle',LineStyle,'LineWidth',LineWidth,...
           'HandleVisibility','off');
    end
    
    % axis formatting
    axis equal
    grid on
    if strcmpi(planet,'Earth Coastlines') || strcmpi(planet,...
            'Blank')
        ax = gca;
        ax.GridColor = [0.35,0.35,0.35];
        ax.GridAlpha = 1;
        ax.GridLineStyle = ':';
    end
    xlim([-180,180]);
    xticks(-180:30:180);
    ylim([-90,90]);
    yticks(-90:30:90);
    xlabel('Longitude, $\lambda\;[^{\circ}]$','Interpreter','latex',...
        'FontSize',18);
    ylabel('Latitude, $\phi\;[^{\circ}]$','Interpreter','latex',...
        'FontSize',18);
    
end

