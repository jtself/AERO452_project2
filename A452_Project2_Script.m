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

% For hammerSat d~0.5 and m~50 but that makes the pert lower 
% RANDOM PARAMETERS
diameter = 1;          % m
mass = 100;            % kg
area = (1/(1000^2))*pi*(diameter/2)^2; % km
Cd = 2.2;
Cr = 1.2;
Psr = 4.57*10^-6;

tfinal = 365*24*60*60; % RANDOM 
tspan = [0 tfinal]; 
ticStart = tic;
options = odeset('RelTol', 1e-9, 'AbsTol',1e-9,'Events',@eventDeOrbit);
init = [SC.init.h; SC.init.ecc; SC.init.TA; SC.init.raan ;SC.init.inc ;SC.init.w]; 
[time, state] = ode45(@vop_ODE, tspan, init, options,wEarth, re, mu, muSun, Cd, area, mass, SC.init.jd, Cr, Psr); 


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

tocEnd = toc(ticStart);
disp("HammerSat took: " + tocEnd + " sec to run")

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

%% Adding Lamberts
% This will run but is incorrect

days = 5;
timeVector = linspace(0,days*24*60*60,days);
tm = 1;
mu = 398600;
dt = 3600; % 1 minute lambert
tspanValue = 24*60*60;
deltaV = 0;

init_vop = [SC.init.h; SC.init.ecc; SC.init.TA; SC.init.raan ;SC.init.inc ;SC.init.w]; 
init_coast =[SC.init.rVect;SC.init.vVect];
options = odeset('RelTol', 1e-8, 'AbsTol',1e-8);

rStorage = []; 
tStorage = [];
for i = 1:length(timeVector)

tspan = [tspanValue*i - tspanValue; tspanValue*i]; % run for 1 day
[time_vop, state_vop] = ode45(@vop_ODE, tspan, init_vop, options, wEarth, re, mu, muSun, Cd, area, mass, SC.init.jd, Cr, Psr); 
[time_coast, state_coast] = ode45(@coast_ODE, tspan, init_coast, options, mu);

[r0vec,v0vec] = COES2RandV(state_vop(end,1),state_vop(end,2),state_vop(end,5),state_vop(end,4),state_vop(end,6),state_vop(end,3),mu);

r1vec = r0vec;
r2vec = state_coast(end,1:3)';

[v1vec, v2vec] = lambert(r1vec, r2vec, dt, tm, mu);

deltaV = deltaV + norm(v0vec-v1vec);

rinit = [state_coast(end,1);state_coast(end,2);state_coast(end,3)];
vinit = [state_coast(end,4);state_coast(end,5);state_coast(end,6)];

[h, inc, RAAN, ecc, w, theta, epsilon, a, T] = OrbitalElements(rinit,vinit,mu);
init_vop = [h; ecc; theta; RAAN; inc; w];
init_coast = [rinit;vinit];

rStorage = [rStorage;state_vop];
tStorage = [tStorage;time_vop];

end
%%
for i = 1:length(rStorage)
    [r_temp,~] = COES2RandV(rStorage(i,1),rStorage(i,2),rStorage(i,5),rStorage(i,4),rStorage(i,6),rStorage(i,3),mu);
    r(i,1:3) = r_temp;
    posNorm(i) = norm(r_temp);
end

%%

figure
plot(tStorage,rStorage(:,5))


%% Section 2