The purpose of this packet is to provide the (next) researcher all derived code from my time (2021-2023) working with Dr. Nandeesh Hiremath on the "Virtual Aperture Multispectral Imaging for Atmospheric Reentry Studies using High Altitude Reflective Arrays" project. The .m files included in this repository constitute the state of the art of this project as of Spring 2023. 
The "groundtrack testing" file takes user-defined TLE data and produces (using Tamas Kis groundtrack() code) groundtrack plot for simple and quick visualizations.
The "CubeSat Reentry Plot" file walks through a simulated 1U CubeSat reentry, making several general assumptions. The scenario initial conditions are: 
W = 1.33; % kg; typical
Cd = 2.22; % typical value
r = sqrt(0.05); % hypot of 5x5 cm triangle (10cm x 10cm 1U CS face)
S = 4*pi*r^2; % treat as sphere for tumbling (worst case)
Beta = W/(Cd*S); % N/m2
h0 = 100e3; % m; starting height
v0 = 7.5e3; % m/s; reentry velocity (LEO)

The "new_EOMs" file will be useful for implementing a more robust reentry model that takes flight path angle (and its derivative) into account. 
The "Example.m" file shows a useful example for plotting various ballistic coefficients (tumbling, etc) along with varying initial flight path angles for reentry. 
  This simulation was used as a homework assignment for Dr. Hiremath's "Reentry Aerodynamics" course in 2024. The problem includes the simulated deployment of a heatshield and characterizes the flight trajectory (2D) from 100 km to earth's surface.
