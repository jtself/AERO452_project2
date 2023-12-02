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

%% Section 1: Outfit1 Calculations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TLE Input (HammerSat)
% SC.init.UTC = "11/18/2023 06:41:13"; % UTC
% SC.init.ecc = 0.0001978;
% SC.init.inc = deg2rad(51.6407);               % rad
% SC.init.rp = 411 + re;                        % km
% SC.init.ra = 414 + re;                        % km
% SC.init.raan = deg2rad(285.1807);             % rad
% SC.init.w = deg2rad(121.7946);                % rad
% SC.init.revperday = 15.51474585;              % rev/day
% SC.init.Me = deg2rad(238.3237);               % rad

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
SC.init.h = findh(SC.init.a, mu, SC.init.ecc);

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

tfinal = (1)*24*3600; % RANDOM 
tspan = [0 tfinal]; 
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

figure
h1 = gca;
earth_sphere(h1)
hold on
plot3(r(:,1),r(:,2),r(:,3),'.')
plot3(r(1,1),r(1,2),r(1,3),'*','LineWidth',3)
plot3(r(end,1),r(end,2),r(end,3),'*','LineWidth',3)
legend("Earth","Orbital Path","Start Position","End Position")
xlabel("X [Km]")
ylabel("Y [Km]")
zlabel("Z [Km]")



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
xlabel("Time [days]")
ylabel("")
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


%% Section 1: Outfit1 Lamberts

r1vec = [r(end,1);r(end,2);r(end,3)];
vend = [v(end,1);v(end,2);v(end,3)];
dt = 500; % 60 sec to return
tm = -1;

tspan = [0 SC.init.T]; 
init_coast = [SC.init.rVect;SC.init.vVect];
options = odeset('RelTol', 1e-8, 'AbsTol',1e-8);
[time_coast, state_coast] = ode45(@coast_ODE, tspan, init_coast, options, mu);
deltaV = zeros(1,length(state_coast));

for i = 1:length(state_coast)
    rf = [state_coast(i,1);state_coast(i,2);state_coast(i,3)];
    [v1vec, v2vec] = lambert(r1vec, rf, dt, tm, mu);
    vf = [state_coast(i,4);state_coast(i,5);state_coast(i,6)];
    deltaV(i) = norm(v1vec - vend) + norm(v2vec - vf);%%%%%%%%%%%%%%%%%%%%%%%%%%% NEED TO CHECK %%%%%%%%%%%%%%%%%%%
end

[value,index] = min(deltaV);
rf = [state_coast(index,1);state_coast(index,2);state_coast(index,3)];
[v1vec, v2vec] = lambert(r1vec, rf, dt, tm, mu);
vf = [state_coast(index,4);state_coast(index,5);state_coast(index,6)];

tspan = [0,dt];
init_coast = [r1vec;v1vec];
[time_lambert, state_lambert] = ode45(@coast_ODE, tspan, init_coast, options, mu);
close all

%
figure
h1 = gca;
%earth_sphere(h1)
%hold on
plot3(r(16000:end,1),r(16000:end,2),r(16000:end,3),'LineWidth',1)
hold on
plot3(state_lambert(:,1),state_lambert(:,2),state_lambert(:,3))
plot3(state_coast(:,1),state_coast(:,2),state_coast(:,3),'LineWidth',1)
plot3(state_lambert(end,1),state_lambert(end,2),state_lambert(end,3),'*','LineWidth',2)
plot3(state_lambert(1,1),state_lambert(1,2),state_lambert(1,3),'*','LineWidth',2)

%% Section 1: Falcon 9 R/B Calculations

% Housekeeping
clear all; close all; clc;

% Governing Constants
re = 6378;                 % km
mu = 398600;               % km3/s2
wEarth = [0;0;72.9211e-6]; % rad/s
muSun = 132.712e9;         % km3/s2

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TLE Input (Falcon 9 R/B)
RB.init.UTC = "11/29/2023 06:19:47"; % UTC
RB.init.ecc = 0.7883697; 
RB.init.inc = deg2rad(16.3434);      % rad
RB.init.rp = 313 + re;               % km
RB.init.ra = 50162 + re;             % km
RB.init.raan = deg2rad(73.9241);     % rad
RB.init.w = deg2rad(293.5891);       % rad
RB.init.revperday = 1.54435022;      % rev/day
RB.init.Me = deg2rad(6.2424);        % rad

% Calculate Additional COEs
RB.init.a = 0.5*(RB.init.rp + RB.init.ra);             % km
RB.init.T = ( (2*pi) / (sqrt(mu)) ) * RB.init.a^(3/2); % sec
RB.init.n = RB.init.revperday*(2*pi/(24*60*60));       % rad/sec
RB.init.jd = juliandate(datetime(RB.init.UTC,"Format","MM/dd/uuuu HH:mm:ss"));
RB.init.E = newtonsKepler(RB.init.Me, RB.init.ecc);
RB.init.TA = 2*atan((tan(RB.init.E/2))...
    /(sqrt((1 - RB.init.ecc)/(1 + RB.init.ecc))));
RB.init.h = findh(RB.init.a, mu, RB.init.ecc);

% Determine epoch r and v
[RB.init.rVect,RB.init.vVect] = COES2RandV(RB.init.h,RB.init.ecc,...
    RB.init.inc,RB.init.raan,RB.init.w,RB.init.TA,mu);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Need to model EXP Drag, J2-J6, Sun as extra Body, and SRP

% RANDOM PARAMETERS
% Height	13.8 m / 45.3 ft
% Diameter	3.7 m / 12.1 ft
% Empty Mass	3,900 kg / 8,598 lb
diameter = 13.8;        % m
mass = 3900;            % kg
area = (1/(1000^2))*pi*(diameter/2)^2; % km2
Cd = 2.2;
Cr = 1.2;
Psr = 4.57*10^-6;

tfinal = (10)*24*3600; % RANDOM 
tspan = [0 tfinal]; 
ticStart = tic;
options = odeset('RelTol', 1e-12, 'AbsTol',1e-12,'Events',@eventDeOrbit);
init = [RB.init.h; RB.init.ecc; RB.init.TA; RB.init.raan ;RB.init.inc ;RB.init.w]; 
[time, state] = ode45(@vop_ODE, tspan, init, options,wEarth, re, mu, muSun, Cd, area, mass, RB.init.jd, Cr, Psr); 


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
disp("Falcon 9 R/B took: " + tocEnd + " sec to run")

figure
h1 = gca;
earth_sphere(h1)
hold on
plot3(r(:,1),r(:,2),r(:,3),'.')
plot3(r(end,1),r(end,2),r(end,3),'*','LineWidth',5)

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
plot(time,state(:,1)- RB.init.h)
title("h")
figure
plot(time,state(:,2)-RB.init.ecc)
title("ecc")
figure
plot(time,rad2deg(state(:,3)-RB.init.TA))
title("theta")
figure
plot(time,rad2deg(state(:,4)-RB.init.raan))
title("raan")
figure
plot(time,deg2rad(state(:,5)-RB.init.inc))
title("inc")
figure
plot(time,deg2rad(state(:,6)-RB.init.w))
title("w")


%% Section 1: Outfit1 Lamberts
clc

r1vec = [r(end,1);r(end,2);r(end,3)];
vend = [v(end,1);v(end,2);v(end,3)];
dt = 500; 
tm = 1;

tspan = [0 RB.init.T]; 
init_coast = [RB.init.rVect;RB.init.vVect];
options = odeset('RelTol', 1e-10, 'AbsTol',1e-10);
[time_coast, state_coast] = ode45(@coast_ODE, tspan, init_coast, options, mu);
deltaV = zeros(1,length(state_coast));

%for i = 1:length(state_coast)
% 600-615
    i = 608;
    rf = [state_coast(i,1);state_coast(i,2);state_coast(i,3)];
    [v1vec, v2vec] = lambert(r1vec, rf, dt, tm, mu);
    vf = [state_coast(i,4);state_coast(i,5);state_coast(i,6)];
    deltaV(i) = norm(v1vec - vend) + norm(v2vec - vf);%%%%%%%%%%%%%%%%%%%%%%%%%%% NEED TO CHECK %%%%%%%%%%%%%%%%%%%
%end

disp(nonzeros(deltaV))
%[value,index] = min(deltaV);
index = i;
rf = [state_coast(index,1);state_coast(index,2);state_coast(index,3)];
[v1vec, v2vec] = lambert(r1vec, rf, dt, tm, mu);
vf = [state_coast(index,4);state_coast(index,5);state_coast(index,6)];

tspan = [0,dt];
init_coast = [r1vec;v1vec];
[time_lambert, state_lambert] = ode45(@coast_ODE, tspan, init_coast, options, mu);
close all

%
figure
h1 = gca;
%earth_sphere(h1)
%hold on
plot3(r(16000:end,1),r(16000:end,2),r(16000:end,3),'LineWidth',1)
hold on
plot3(state_lambert(:,1),state_lambert(:,2),state_lambert(:,3))
plot3(state_coast(:,1),state_coast(:,2),state_coast(:,3),'LineWidth',1)
plot3(state_lambert(end,1),state_lambert(end,2),state_lambert(end,3),'*','LineWidth',2)
plot3(state_lambert(1,1),state_lambert(1,2),state_lambert(1,3),'*','LineWidth',2)
