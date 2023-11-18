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

%% Section 1: HammerSat

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TLE Input (HammerSat)
SC.init.UTC = "11/18/2023 06:41:13"; % UTC
SC.init.ecc = 0.0001978; % actual ecc but we must change it to zero
SC.init.inc = deg2rad(51.6407);               % rad
SC.init.rp = 411 + re;                 % km
SC.init.ra = 414 + re;                 % km
SC.init.raan = deg2rad(285.1807);            % rad
SC.init.w = deg2rad(121.7946);               % rad
SC.init.revperday = 15.51474585;              % rev/day
SC.init.Me = deg2rad(238.3237);              % rad

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
options = odeset('RelTol', 1e-8, 'AbsTol',1e-8,'Events',@eventDeOrbit);
init = [SC.init.h; SC.init.ecc; SC.init.TA; SC.init.raan ;SC.init.inc ;SC.init.w]; 
[time, state] = ode45(@vop_ODE, tspan, init, options,wEarth, re, mu, muSun, Cd, area, mass, SC.init.jd); 
tocTest1 = toc(tic);

disp("Part 1 took: " + tocTest1 + " sec")


time = time/(24*3600);
% Find rVector
r = zeros(length(state),3);
for i = 1:length(state)
    [r_temp,~] = COES2RandV(state(i,1),state(i,2),state(i,5),state(i,4),state(i,6),state(i,3),mu);
    r(i,1:3) = r_temp;
    posNorm(i) = norm(r_temp);
end

[~, apogeeIndex] = findpeaks(posNorm);
[~,perigeeIndex] = findpeaks(-posNorm);
%apogeeIndex = apogeeIndex(1:(length(perigeeIndex)));

for i = 1:length(apogeeIndex)
    apogee(i) = posNorm(apogeeIndex(i));
    perigee(i) = posNorm(perigeeIndex(i));
    timeA(i) = time(apogeeIndex(i));
    timeP(i) = time(perigeeIndex(i));
end

apogee = apogee - re;
perigee = perigee - re;

figure
h1 = gca;
earth_sphere(h1)
hold on
plot3(r(:,1),r(:,2),r(:,3))


figure
plot(timeA,apogee,'LineWidth',2)
hold on
plot(timeP,perigee,'LineWidth',2)
grid on
legend("Apogee","Perigee",'Location','best')
title("HammerSAT Orbital Path")
ylabel("Altitude [km]")
xlabel("Time [Days]")


%% Section 2