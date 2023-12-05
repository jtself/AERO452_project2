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
SC.init.TA = 2*atan((tan(SC.init.E/2))/(sqrt((1 - SC.init.ecc)/(1 + SC.init.ecc))));
SC.init.h = findh(SC.init.a, mu, SC.init.ecc);
SC.init.h = sqrt(SC.init.rp*mu*(1 + SC.init.ecc* cos(SC.init.TA)));


% Determine epoch r and v
[SC.init.rVect,SC.init.vVect] = COES2RandV(SC.init.h,SC.init.ecc,...
    SC.init.inc,SC.init.raan,SC.init.w,SC.init.TA,mu);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Need to model EXP Drag, J2-J6, Sun as extra Body, and SRP

% This is a 1U Cubesat
% 10x10x10 cm cube
% mass = 2 kg
diameter = sqrt(0.1^2 + 0.1^2);        % m
mass = 2;            % kg
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
plot3(r(1,1),r(1,2),r(1,3),'*','LineWidth',5)
plot3(r(end,1),r(end,2),r(end,3),'*','LineWidth',5)
lgd = legend("Earth","Orbital Path","Start Position","End Position",'Location','southoutside');
lgd.NumColumns = 2;
xlabel("X [Km]")
ylabel("Y [Km]")
zlabel("Z [Km]")

figure
plot(time,ra-re,'LineWidth',2)
hold on
plot(time,rp-re,'LineWidth',2)
% Graph pretty 
ylim padded 
xlim tight 
xLab = xlabel('Time [days]','Interpreter','latex'); 
yLab = ylabel('Altitude [km]','Interpreter','latex'); 
%plotTitle = title('title $x$','interpreter','latex'); 
%set(plotTitle,'FontSize',14,'FontWeight','bold') 
set(gca,'FontName','Palatino Linotype') 
set([xLab, yLab],'FontName','Palatino Linotype') 
set(gca,'FontSize', 9) 
set([xLab, yLab],'FontSize', 9) 
grid on 
legend('Apogee','Perigee', 'interpreter','latex','Location', 'best')


figure()
tiledlayout("vertical")
nexttile
plot(time,state(:,1)- SC.init.h)
% Graph pretty 
ylim padded 
xlim tight  
yLab = ylabel('$h - h_0$','Interpreter','latex');  
set(gca,'FontName','Palatino Linotype') 
set(yLab,'FontName','Palatino Linotype') 
set(gca,'FontSize', 9) 
set(yLab,'FontSize', 9) 
grid on 

nexttile
plot(time,state(:,2)-SC.init.ecc)
% Graph pretty 
ylim padded 
xlim tight  
yLab = ylabel('$e - e_0$','Interpreter','latex');  
set(gca,'FontName','Palatino Linotype') 
set(yLab,'FontName','Palatino Linotype') 
set(gca,'FontSize', 9) 
set(yLab,'FontSize', 9) 
grid on 


%nexttile DONT PLOT THETA
%plot(time,rad2deg(state(:,3)-SC.init.TA))
% Graph pretty 
%ylim padded 
%xlim tight 
%xLab = xlabel('Time [days]','Interpreter','latex'); 
%yLab = ylabel('$\theta - \omega_0$ [degs]','Interpreter','latex'); 
%plotTitle = title('title $x$','interpreter','latex'); 
%set(plotTitle,'FontSize',14,'FontWeight','bold') 
%set(gca,'FontName','Palatino Linotype') 
%set([xLab, yLab],'FontName','Palatino Linotype') 
%set(gca,'FontSize', 9) 
%set([xLab, yLab],'FontSize', 14) 
%grid on 



nexttile
plot(time,rad2deg(state(:,4)-SC.init.raan))
% Graph pretty 
ylim padded 
xlim tight 
yLab = ylabel('$\Omega - \Omega_0$','Interpreter','latex'); 
set(gca,'FontName','Palatino Linotype') 
set(yLab,'FontName','Palatino Linotype') 
set(gca,'FontSize', 9) 
set(yLab,'FontSize', 9) 
grid on 

nexttile
plot(time,state(:,5)-SC.init.inc)
% Graph pretty 
ylim padded 
xlim tight 
yLab = ylabel('$inc - inc_0$','Interpreter','latex');  
set(gca,'FontName','Palatino Linotype') 
set(yLab,'FontName','Palatino Linotype') 
set(gca,'FontSize', 9) 
set(yLab,'FontSize', 9) 
grid on 

nexttile
plot(time,rad2deg(state(:,6)-SC.init.w))
% Graph pretty 
ylim padded 
xlim tight 
xLab = xlabel('Time [days]','Interpreter','latex'); 
yLab = ylabel('$\omega - \omega_0$','Interpreter','latex'); 
set(gca,'FontName','Palatino Linotype') 
set([xLab, yLab],'FontName','Palatino Linotype') 
set(gca,'FontSize', 9) 
set([xLab, yLab],'FontSize', 9) 
grid on 

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
    deltaV(i) = norm(v1vec - vend) + norm(v2vec - vf);
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
earth_sphere(h1)
hold on
plot3(r(:,1),r(:,2),r(:,3),'LineWidth',1)
plot3(state_coast(:,1),state_coast(:,2),state_coast(:,3),'LineWidth',1)
plot3(state_lambert(:,1),state_lambert(:,2),state_lambert(:,3),'LineWidth',4)
plot3(state_lambert(1,1),state_lambert(1,2),state_lambert(1,3),'*','LineWidth',5)
plot3(state_lambert(end,1),state_lambert(end,2),state_lambert(end,3),'*','LineWidth',5)
lgd = legend("Earth","Perturbed Orbit","Osculating Orbit","Lambert Trajectory","Initial Position","Final Position",'Location','southoutside');
lgd.NumColumns = 2;

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
diameter = 3.7;      % m
height = 13.8; % m
mass = 3900;            % kg
area = (1/(1000^2))*(diameter*height); % km2
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
plot3(r(1,1),r(1,2),r(1,3),'*','LineWidth',5)
plot3(r(end,1),r(end,2),r(end,3),'*','LineWidth',5)
lgd = legend("Earth","Orbital Path","Start Position","End Position",'Location','southoutside');
lgd.NumColumns = 2;
xlabel("X [Km]")
ylabel("Y [Km]")
zlabel("Z [Km]")

figure
plot(time,ra-re,'LineWidth',2)
hold on
plot(time,rp-re,'LineWidth',2)

% Graph pretty 
ylim padded 
xlim tight 
xLab = xlabel('Time [days]','Interpreter','latex'); 
yLab = ylabel('Altitude [km]','Interpreter','latex'); 
%plotTitle = title('title $x$','interpreter','latex'); 
%set(plotTitle,'FontSize',14,'FontWeight','bold') 
set(gca,'FontName','Palatino Linotype') 
set([xLab, yLab],'FontName','Palatino Linotype') 
set(gca,'FontSize', 9) 
set([xLab, yLab],'FontSize', 9) 
grid on 
legend('Apogee','Perigee', 'interpreter','latex','Location', 'best')



% Plot COEs
figure()
tiledlayout("vertical")
nexttile
plot(time,state(:,1)- RB.init.h)
% Graph pretty 
ylim padded 
xlim tight  
yLab = ylabel('$h - h_0$','Interpreter','latex');  
set(gca,'FontName','Palatino Linotype') 
set(yLab,'FontName','Palatino Linotype') 
set(gca,'FontSize', 9) 
set(yLab,'FontSize', 9) 
grid on 

nexttile
plot(time,state(:,2)-RB.init.ecc)
% Graph pretty 
ylim padded 
xlim tight  
yLab = ylabel('$e - e_0$','Interpreter','latex');  
set(gca,'FontName','Palatino Linotype') 
set(yLab,'FontName','Palatino Linotype') 
set(gca,'FontSize', 9) 
set(yLab,'FontSize', 9) 
grid on 

nexttile
plot(time,rad2deg(state(:,4)-RB.init.raan))
% Graph pretty 
ylim padded 
xlim tight 
yLab = ylabel('$\Omega - \Omega_0$','Interpreter','latex'); 
set(gca,'FontName','Palatino Linotype') 
set(yLab,'FontName','Palatino Linotype') 
set(gca,'FontSize', 9) 
set(yLab,'FontSize', 9) 
grid on 

nexttile
plot(time,state(:,5)-RB.init.inc)
% Graph pretty 
ylim padded 
xlim tight 
yLab = ylabel('$inc - inc_0$','Interpreter','latex');  
set(gca,'FontName','Palatino Linotype') 
set(yLab,'FontName','Palatino Linotype') 
set(gca,'FontSize', 9) 
set(yLab,'FontSize', 9) 
grid on 

nexttile
plot(time,rad2deg(state(:,6)-RB.init.w))
% Graph pretty 
ylim padded 
xlim tight 
xLab = xlabel('Time [days]','Interpreter','latex'); 
yLab = ylabel('$\omega - \omega_0$','Interpreter','latex'); 
set(gca,'FontName','Palatino Linotype') 
set([xLab, yLab],'FontName','Palatino Linotype') 
set(gca,'FontSize', 9) 
set([xLab, yLab],'FontSize', 9) 
grid on 

%% Section 1: Falcon 9 R/B Lamberts
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
    deltaV(i) = norm(v1vec - vend) + norm(v2vec - vf);
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
earth_sphere(h1)
hold on
plot3(r(:,1),r(:,2),r(:,3),'LineWidth',1)
plot3(state_coast(:,1),state_coast(:,2),state_coast(:,3),'LineWidth',1)
plot3(state_lambert(:,1),state_lambert(:,2),state_lambert(:,3),'LineWidth',4)
plot3(state_lambert(1,1),state_lambert(1,2),state_lambert(1,3),'*','LineWidth',5)
plot3(state_lambert(end,1),state_lambert(end,2),state_lambert(end,3),'*','LineWidth',5)
lgd = legend("Earth","Perturbed Orbit","Osculating Orbit","Lambert Trajectory","Initial Position","Final Position",'Location','southoutside');
lgd.NumColumns = 2;
%% Section 2: Continous Corrections for Outfti 1 

clear all
close all
clc
% Governing Constants
re = 6378;                 % km
mu = 398600;               % km3/s2
wEarth = [0;0;72.9211e-6]; % rad/s
muSun = 132.712e9;         % km3/s2

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
SC.init.TA = 2*atan((tan(SC.init.E/2))/(sqrt((1 - SC.init.ecc)/(1 + SC.init.ecc))));
SC.init.h = findh(SC.init.a, mu, SC.init.ecc);

% Determine epoch r and v
[SC.init.rVect,SC.init.vVect] = COES2RandV(SC.init.h,SC.init.ecc,...
    SC.init.inc,SC.init.raan,SC.init.w,SC.init.TA,mu);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Need to model EXP Drag, J2-J6, Sun as extra Body, and SRP

% This is a 1U Cubesat
% 10x10x10 cm cube
% mass = 2 kg
diameter = sqrt(0.1^2 + 0.1^2);        % m
mass = 2;            % kg
area = (1/(1000^2))*pi*(diameter/2)^2; % km2
Cd = 2.2;
Cr = 1.2;
Psr = 4.57*10^-6;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
h = SC.init.h;
ecc = SC.init.ecc;
TA = SC.init.TA;
raan = SC.init.raan;
inc = SC.init.inc;
w = SC.init.w;
tm = -1;
overallTime = [];
overallState = [];
overallDeltaV = [];
initalPoints = [];
finalPoints = [];
initialTime = [];
finalTime = [];
init_pert = [h; ecc; TA; raan ; inc ;w]; 
init_noPert = [h; ecc; TA; raan ;inc ;w];
init = [init_pert;init_noPert];
options = odeset('RelTol', 1e-12, 'AbsTol',1e-12,'Events',@eventDeOrbit);


perTime = 1.35*24*3600;
tfinal = perTime;

tspan = [0 tfinal]; 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for i = 1:5

[time, state] = ode45(@vop_ODE_Part2, tspan, init, options,wEarth, re, mu, muSun, Cd, area, mass, SC.init.jd, Cr, Psr); 


for j = 1:length(state)
    if   time(end)-time(length(state)-j) > 500
        dt = time(end)-time(length(state)-j);
        break
    end
end


[r1vec,v1] = COES2RandV(state((length(state)-j),1)...
    ,state((length(state)-j),2)...
    ,state((length(state)-j),5)...
    ,state((length(state)-j),4)...
    ,state((length(state)-j),6)...
    ,state((length(state)-j),3)...
    ,mu);

[r2vec,v2] = COES2RandV(state(end,7),state(end,8),state(end,11),state(end,10),state(end,12),state(end,9),mu);

[v1vec, v2vec] = lambert(r1vec, r2vec, dt, tm, mu);

deltaV = norm(v1-v1vec) + norm(v2-v2vec);

time = time + perTime*i - perTime;

overallDeltaV = [overallDeltaV;deltaV];
overallTime   = [overallTime;time];
overallState  = [overallState;state];

init_pert   = [state(end,7);state(end,8);state(end,9);state(end,10);state(end,11);state(end,12)];
init_noPert = [state(end,7);state(end,8);state(end,9);state(end,10);state(end,11);state(end,12)]; 


finalPertLocation = [state(end,1);state(end,2);state(end,3);state(end,4);state(end,5);state(end,6)];
initPertLocation  = [state(1,1);state(1,2);state(1,3);state(1,4);state(1,5);state(1,6)];

init         = [init_pert;init_noPert];
initalPoints = [initalPoints;initPertLocation'];
initialTime  = [initialTime;time(1)];
finalPoints  = [finalPoints;finalPertLocation'];
finalTime    = [finalTime;time(end)];


end

overallTime = overallTime/24/3600;
initialTime = initialTime/24/3600;
finalTime= finalTime/24/3600;   
%%% PRETTY PLOTS
figure
plot(overallTime,overallState(:,1)- SC.init.h)
hold on
plot(overallTime,overallState(:,7)- SC.init.h,'LineWidth',2)
plot(initialTime,initalPoints(:,1)-SC.init.h,'*','LineWidth',3)
plot(finalTime,finalPoints(:,1)-SC.init.h,'*','LineWidth',3)
% Graph pretty 
ylim padded 
xlim tight 
xLab = xlabel('Time [days]','Interpreter','latex'); 
yLab = ylabel('$h - h_0$ [km2/s]','Interpreter','latex'); 
%plotTitle = title('title $x$','interpreter','latex'); 
%set(plotTitle,'FontSize',14,'FontWeight','bold') 
set(gca,'FontName','Palatino Linotype') 
set([xLab, yLab],'FontName','Palatino Linotype') 
set(gca,'FontSize', 9) 
set([xLab, yLab],'FontSize', 10) 
grid on 
legend("Perturbed Orbit","Osculating Orbit","Initial change in $h$","Final change in $h$", 'interpreter','latex','Location', 'best')



figure
plot(overallTime,overallState(:,2)-SC.init.ecc)
hold on
plot(overallTime,overallState(:,8)-SC.init.ecc,'LineWidth',2)
plot(initialTime,initalPoints(:,2)-SC.init.ecc,'*','LineWidth',3)
plot(finalTime,finalPoints(:,2)-SC.init.ecc,'*','LineWidth',3)
% Graph pretty 
ylim padded 
xlim tight 
xLab = xlabel('Time [days]','Interpreter','latex'); 
yLab = ylabel('$ecc-ecc_0$','Interpreter','latex'); 
%plotTitle = title('title $x$','interpreter','latex'); 
%set(plotTitle,'FontSize',14,'FontWeight','bold') 
set(gca,'FontName','Palatino Linotype') 
set([xLab, yLab],'FontName','Palatino Linotype') 
set(gca,'FontSize', 9) 
set([xLab, yLab],'FontSize', 10) 
grid on 
legend("Perturbed Orbit","Osculating Orbit","Initial change in $ecc$","Final change in $ecc$", 'interpreter','latex','Location', 'best')

figure
plot(overallTime,rad2deg(overallState(:,4)-SC.init.raan))
hold on
plot(overallTime,rad2deg(overallState(:,10)-SC.init.raan),'LineWidth',2)
plot(initialTime,rad2deg(initalPoints(:,4)-SC.init.raan),'*','LineWidth',3)
plot(finalTime,rad2deg(finalPoints(:,4)-SC.init.raan),'*','LineWidth',3)
% Graph pretty 
ylim padded 
xlim tight 
xLab = xlabel('Time [days]','Interpreter','latex'); 
yLab = ylabel('$\Omega - \Omega_0$ [deg]','Interpreter','latex'); 
%plotTitle = title('title $x$','interpreter','latex'); 
%set(plotTitle,'FontSize',14,'FontWeight','bold') 
set(gca,'FontName','Palatino Linotype') 
set([xLab, yLab],'FontName','Palatino Linotype') 
set(gca,'FontSize', 9) 
set([xLab, yLab],'FontSize', 10) 
grid on 
legend("Perturbed Orbit","Osculating Orbit","Initial change in $\Omega$","Final change in $\Omega$", 'interpreter','latex','Location', 'best')


figure
plot(overallTime,rad2deg(overallState(:,5)-SC.init.inc))
hold on
plot(overallTime,overallState(:,11)-SC.init.inc,'LineWidth',2)
plot(initialTime,rad2deg(initalPoints(:,5)-SC.init.inc),'*','LineWidth',3)
plot(finalTime,rad2deg(finalPoints(:,5)-SC.init.inc),'*','LineWidth',3)
% Graph pretty 
ylim padded 
xlim tight 
xLab = xlabel('Time [days]','Interpreter','latex'); 
yLab = ylabel('$inc - inc_0$ [deg]','Interpreter','latex'); 
%plotTitle = title('title $x$','interpreter','latex'); 
%set(plotTitle,'FontSize',14,'FontWeight','bold') 
set(gca,'FontName','Palatino Linotype') 
set([xLab, yLab],'FontName','Palatino Linotype') 
set(gca,'FontSize', 9) 
set([xLab, yLab],'FontSize', 10) 
grid on 
legend("Perturbed Orbit","Osculating Orbit","Initial change in $inc$","Final change in $inc$", 'interpreter','latex','Location', 'best')

figure
plot(overallTime,rad2deg(overallState(:,6)-SC.init.w))
hold on
plot(overallTime,rad2deg(overallState(:,12)-SC.init.w),'LineWidth',2)
plot(initialTime,rad2deg(initalPoints(:,6)-SC.init.w),'*','LineWidth',3)
plot(finalTime,rad2deg(finalPoints(:,6)-SC.init.w),'*','LineWidth',3)
% Graph pretty 
ylim padded 
xlim tight 
xLab = xlabel('Time [days]','Interpreter','latex'); 
yLab = ylabel('$\omega - \omega_0$ [deg]','Interpreter','latex'); 
%plotTitle = title('title $x$','interpreter','latex'); 
%set(plotTitle,'FontSize',14,'FontWeight','bold') 
set(gca,'FontName','Palatino Linotype') 
set([xLab, yLab],'FontName','Palatino Linotype') 
set(gca,'FontSize', 9) 
set([xLab, yLab],'FontSize', 10) 
grid on 
legend("Perturbed Orbit","Osculating Orbit","Initial change in $\omega$","Final change in $\omega$", 'interpreter','latex','Location', 'best')

breakpoint = 1;


%%%%% END PRETTY PLOTS 

for i = 1:length(overallState)
    [r_temp,v_temp] = COES2RandV(overallState(i,1),overallState(i,2),overallState(i,5),overallState(i,4),overallState(i,6),overallState(i,3),mu);
    r(i,1:3) = r_temp;
    v(i,1:3) = v_temp;
    h = overallState(i,1);
    ecc = overallState(i,2);
    a = (h^2)/(mu*(1-ecc^2));
    ra(i) = a + a*ecc;
    rp(i) = 2*a - ra(i);
end

figure
plot(overallTime,ra-re,'LineWidth',2)
hold on
plot(overallTime,rp-re,'LineWidth',2)
% Graph pretty 
ylim padded 
xlim tight 
xLab = xlabel('Time [days]','Interpreter','latex'); 
yLab = ylabel('Altitude [km]','Interpreter','latex'); 
%plotTitle = title('title $x$','interpreter','latex'); 
%set(plotTitle,'FontSize',14,'FontWeight','bold') 
set(gca,'FontName','Palatino Linotype') 
set([xLab, yLab],'FontName','Palatino Linotype') 
set(gca,'FontSize', 9) 
set([xLab, yLab],'FontSize', 9) 
grid on 
legend('Apogee','Perigee', 'interpreter','latex','Location', 'best')
%% Section 3: Special Orbit


% TLE Input (Envisat)
SC.init.UTC = "12/04/2023 12:16:40"; % UTC
SC.init.ecc = 0.0001316; % actual ecc but we must change it to zero
SC.init.inc = deg2rad(98.2559);               % rad
SC.init.rp = 762 + re;                 % km
SC.init.ra = 764 + re;                 % km
SC.init.raan = deg2rad(301.7273);            % rad
SC.init.w = deg2rad(83.6424);               % rad
SC.init.revperday = 14.38460028;              % rev/day
SC.init.Me = deg2rad(29.2718);              % rad

% Calculate Additional COEs
SC.init.a = 0.5*(SC.init.rp + SC.init.ra);             % km
SC.init.T = ( (2*pi) / (sqrt(mu)) ) * SC.init.a^(3/2); % sec
SC.init.n = SC.init.revperday*(2*pi/(24*60*60));       % rad/sec
SC.init.jd = juliandate(datetime(SC.init.UTC,"Format","MM/dd/uuuu HH:mm:ss"));
SC.init.E = newtonsKepler(SC.init.Me, SC.init.ecc);
SC.init.TA = 2*atan((tan(SC.init.E/2))/(sqrt((1 - SC.init.ecc)/(1 + SC.init.ecc))));
SC.init.h = findh(SC.init.a, mu, SC.init.ecc);
SC.init.h = sqrt(SC.init.rp*mu*(1 + SC.init.ecc* cos(SC.init.TA)));


% Determine epoch r and v
[SC.init.rVect,SC.init.vVect] = COES2RandV(SC.init.h,SC.init.ecc,...
    SC.init.inc,SC.init.raan,SC.init.w,SC.init.TA,mu);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Need to model EXP Drag, J2-J6, Sun as extra Body, and SRP

diameter = sqrt(0.1^2 + 0.1^2);        % m
mass = 8211;            % kg
area = (1/(1000^2))*26*10; % km2
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
plot3(r(1,1),r(1,2),r(1,3),'*','LineWidth',5)
plot3(r(end,1),r(end,2),r(end,3),'*','LineWidth',5)
lgd = legend("Earth","Orbital Path","Start Position","End Position",'Location','southoutside');
lgd.NumColumns = 2;
xlabel("X [Km]")
ylabel("Y [Km]")
zlabel("Z [Km]")

figure
% something is busted on the line below:
plot(time,ra-re,'LineWidth',2)
hold on
plot(time,rp-re,'LineWidth',2)
% Graph pretty 
ylim padded 
xlim tight 
xLab = xlabel('Time [days]','Interpreter','latex'); 
yLab = ylabel('Altitude [km]','Interpreter','latex'); 
%plotTitle = title('title $x$','interpreter','latex'); 
%set(plotTitle,'FontSize',14,'FontWeight','bold') 
set(gca,'FontName','Palatino Linotype') 
set([xLab, yLab],'FontName','Palatino Linotype') 
set(gca,'FontSize', 9) 
set([xLab, yLab],'FontSize', 9) 
grid on 
legend('Apogee','Perigee', 'interpreter','latex','Location', 'best')

figure
plot(time,state(:,1)- SC.init.h)
ylabel("h-h0 [km2/s]")
xlabel("Time [days]")
grid on
figure
plot(time,state(:,2)-SC.init.ecc)
ylabel("ecc-ecc0")
xlabel("Time [days]")
grid on
figure
plot(time,rad2deg(state(:,3)-SC.init.TA))
ylabel("theta-theta0 [degs]")
xlabel("Time [days]")
grid on
figure
plot(time,rad2deg(state(:,4)-SC.init.raan))
ylabel("raan-raan0 [degs]")
xlabel("Time [days]")
grid on
figure
plot(time,state(:,5)-SC.init.inc)
ylabel("inc-inc0 [degs]")
xlabel("Time [days]")
grid on
figure
plot(time,rad2deg(state(:,6)-SC.init.w))
ylabel("w-w0 [degs]")
xlabel("Time [days]")
grid on
