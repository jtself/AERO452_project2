%% AERO452 | SPACEFLIGHT DYNAMICS II

%{
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
% SC.init.UTC = "11/18/2023 06:41:13"; % UTC
% SC.init.ecc = 0.0001978; % actual ecc but we must change it to zero
% SC.init.inc = deg2rad(51.6407);               % rad
% SC.init.rp = 411 + re;                 % km
% SC.init.ra = 414 + re;                 % km
% SC.init.raan = deg2rad(285.1807);            % rad
% SC.init.w = deg2rad(121.7946);               % rad
% SC.init.revperday = 15.51474585;              % rev/day
% SC.init.Me = deg2rad(238.3237);              % rad

% TLE Input (Outfit1)
SC.init.UTC = "11/29/2023 13:39:51"; % UTC
SC.init.ecc = 0.0063209; % actual ecc but we must change it to zero
SC.init.inc = deg2rad(98.1528);               % rad
SC.init.rp = 363 + re;                 % km
SC.init.ra = 439 + re;                 % km
SC.init.raan = deg2rad(196.9400);            % rad
SC.init.w = deg2rad(224.8522);               % rad
SC.init.revperday = 15.53527537;              % rev/day
SC.init.Me = deg2rad(134.7611);              % rad

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
area = (1/(1000^2))*pi*(diameter/2)^2; % km2
Cd = 2.2;
Cr = 1.2;
Psr = 4.57*10^-6;

tfinal = 10*24*60*60; % RANDOM 
tspan = [0 tfinal]; 
ticStart = tic;
options = odeset('RelTol', 1e-12, 'AbsTol',1e-12,'Events',@eventDeOrbit);
init = [SC.init.h; SC.init.ecc; SC.init.TA; SC.init.raan ;SC.init.inc ;SC.init.w]; 
[time, state] = ode45(@vop_ODE, tspan, init, options,wEarth, re, mu, muSun, Cd, area, mass, SC.init.jd, Cr, Psr); 


time = time/(24*3600);
% Find rVector
r = zeros(length(state),3);
for i = 1:length(state)
    [r_temp,v_temp] = COES2RandV(state(i,1),state(i,2),state(i,5),state(i,4),state(i,6),state(i,3),mu);
    r(i,1:3) = r_temp;
    v(i,1:3) = v_temp;
    h = state(i,1);
    ecc = state(i,2);
    a = (h^2)/(mu*(1-ecc^2));
    ra(i) = a + a*ecc;
    rp(i) = 2*a - ra(i);
end

tocEnd = toc(ticStart);
disp("HammerSat took: " + tocEnd + " sec to run")

figure
h1 = gca;
earth_sphere(h1)
hold on
plot3(r(:,1),r(:,2),r(:,3),'.')

figure
plot(time,ra-re,'LineWidth',2)
hold on
plot(time,rp-re,'LineWidth',2)
grid on
legend("Apogee","Perigee",'Location','best')
title("HammerSAT Orbital Path")
ylabel("Altitude [km]")
xlabel("Time [Days]")

figure
plot(time,state(:,1)- SC.init.h)
title("h")
figure
plot(time,state(:,2)-SC.init.ecc)
title("ecc")
figure
plot(time,rad2deg(state(:,3)-SC.init.TA))
title("theta")
figure
plot(time,rad2deg(state(:,4)-SC.init.raan))
title("raan")
figure
plot(time,deg2rad(state(:,5)-SC.init.inc))
title("inc")
figure
plot(time,deg2rad(state(:,6)-SC.init.w))
title("w")


%% Time for Lambert Back to Orginal Spot 

r1vec = [r(end,1);r(end,2);r(end,3)];
vend = [v(end,1);v(end,2);v(end,3)];
dt = 3600; % 60 sec to return
tm = 1;

tspan = [0 tfinal]; 
init_coast = [SC.init.rVect;SC.init.vVect];
options = odeset('RelTol', 1e-9, 'AbsTol',1e-9);
[time_coast, state_coast] = ode45(@coast_ODE, tspan, init_coast, options, mu);
deltaV = zeros(1,length(state_coast));


for i = 1:length(state_coast)
    rf = [state_coast(i,1);state_coast(i,2);state_coast(i,3)];
    [v1vec, v2vec] = lambert(r1vec, rf, dt, tm, mu);
    vf = [state_coast(i,4);state_coast(i,5);state_coast(i,6)];
    deltaV(i) = norm(v1vec - vend) + norm(v2vec - vf);
end


[value,index] = min(deltaV);
rf = [state_coast(index,1);state_coast(index,2);state_coast(index,3)];
[v1vec, v2vec] = lambert(r1vec, rf, dt, tm, mu);
vf = [state_coast(index,4);state_coast(index,5);state_coast(index,6)];

tspan = [0,dt];
init_coast = [r1vec;v1vec];
[time_lambert, state_lambert] = ode45(@coast_ODE, tspan, init_coast, options, mu);

figure
h1 = gca;
earth_sphere(h1)
hold on
plot3(state_lambert(:,1),state_lambert(:,2),state_lambert(:,3))
plot3(state_coast(:,1),state_coast(:,2),state_coast(:,3))
plot3(state_lambert(end,1),state_lambert(end,2),state_lambert(end,3),'*','LineWidth',2)
plot3(state_lambert(1,1),state_lambert(1,2),state_lambert(1,3),'*','LineWidth',2)
plot3(r(:,1),r(:,2),r(:,3),'.')


