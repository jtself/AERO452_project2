%{
AERO452 | SPACEFLIGHT DYNAMICS II
Group Project #2
Authors: 
    Travis Bouck
    Justin Self
Due date: 
    Dec. 8, 2023
%}

% Right now this function takes 0.0001069 s to run for 1 year...good or
% bad? idk 

% Housekeeping
clear all; close all; clc;

% Add Path
addpath("Functions\")

% Governing Constants
re = 6378;                 % km
mu = 398600;               % km3/s2
wEarth = [0;0;72.9211e-6]; % rad/s
muSun = 132.712e9;         % km3/s2

%% Section 1

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TLE Input (RANDOM TLE- GEO from other Proj)
SC.init.UTC = "02/10/2023 05:38:01"; % UTC
SC.init.ecc = 0.000001016; % actual ecc but we must change it to zero
SC.init.inc = deg2rad(0.0009);               % rad
SC.init.rp = 35782 + re;                 % km
SC.init.ra = 35791 + re;                 % km
SC.init.raan = deg2rad(246.2601);            % rad
SC.init.w = deg2rad(273.8346);               % rad
SC.init.revperday = 1.00269599;              % rev/day
SC.init.Me = deg2rad(199.9092);              % rad

% Calculate Additional COEs
SC.init.a = 0.5*(SC.init.rp + SC.init.ra);             % km
SC.init.T = ( (2*pi) / (sqrt(mu)) ) * SC.init.a^(3/2); % sec
SC.init.n = SC.init.revperday*(2*pi/(24*60*60));       % rad/sec
SC.init.jd = juliandate(datetime(SC.init.UTC,"Format","MM/dd/uuuu HH:mm:ss"));

SC.init.E = newtonsKepler(SC.init.Me, SC.init.ecc);
SC.init.TA = 2*atan((tan(SC.init.E/2))...
    /(sqrt((1 - SC.init.ecc)/(1 + SC.init.ecc))));
SC.init.h = findh(SC.init.a, mu, SC.init.ecc, SC.init.TA);

% Determine epoch r and v
[SC.init.rVect,SC.init.vVect] = COES2RandV(SC.init.h,SC.init.ecc,...
    SC.init.inc,SC.init.raan,SC.init.w,SC.init.TA,mu);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Need to model EXP Drag, J2-J6, Sun as extra Body, and SRP

% RANDOM PARAMETERS
diameter = 1;          % m
mass = 100;            % kg
area = (1/(1000^2))*pi*(diameter/2)^2; % km
Cd = 2.2;
tfinal = 365*24*60*60; % RANDOM 

tspan = [0 tfinal]; 
tic 
options = odeset('RelTol', 1e-13, 'AbsTol',1e-13); % No event function needed cause were not deorbiting
init = [SC.init.h; SC.init.ecc; SC.init.TA; SC.init.raan ;SC.init.inc ;SC.init.w]; 
[time, state] = ode45(@vop_ODE, tspan, init, options,wEarth, re, mu, muSun, Cd, area, mass, SC.init.jd); 
tocTest1 = toc(tic);

disp("Part 1 took: " + tocTest1 + " sec")

% Find rVector
r = zeros(length(state),3);
for i = 1:length(state)
    [r_temp,~] = COES2RandV(state(i,1),state(i,2),state(i,5),state(i,4),state(i,6),state(i,3),mu);
    r(i,1:3) = r_temp;
end

plot3(r(:,1),r(:,2),r(:,3))


%% Section 2