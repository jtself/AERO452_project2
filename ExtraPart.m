% Housekeeping
clear all; close all; clc;

%{
        |
       / \
      / _ \
     |.o '.|
     |'._.'|
     |     |
   ,'|  |  |`.
  /  |  |  |  \
  |,-'--|--'-.| 

%}

% Add Path
addpath("Functions\")

% Governing Constants
rm = 3396;                         % km
mu = 42828;                        % km3/s2
wMars = [0;0;0.0000708907088544];   % rad/s
muSun = 132.712e9;                 % km3/s2


UTC = "2/04/2024 00:00:00"; % UTC
jd_epoch = juliandate(datetime(UTC,"Format","MM/dd/uuuu HH:mm:ss"));

mass = 1031;                  % kg
area = (1/(1000^2))*11.4*2.3; % km2
Cd = 2.2;
Cr = 1.2;
Psr = 586.2/(2.998e8);

tfinal = 50*24*3600; % RANDOM 
tspan = [0 tfinal]; 

h0 = 1.4892e4;
ecc0 = 0.4604;
inc0 = deg2rad(75);
TA0 = 0;
raan0 = 0;
w0 =  0;

options = odeset('RelTol', 1e-12, 'AbsTol',1e-12,'Events',@eventDeOrbit_Mars);
init = [h0; ecc0; TA0; raan0; inc0; w0]; 
[time, state] = ode45(@vop_mars_ode, tspan, init, options,wMars, rm, mu, muSun, Cd, area, mass, jd_epoch, Cr, Psr); 

time = time/(24*3600);

for i = 1:length(state)
    [r_temp,v_temp] = COES2RandV(state(i,1),state(i,2),state(i,5),state(i,4),state(i,6),state(i,3),mu);
    r(i,1:3) = r_temp;
    v(i,1:3) = v_temp;
    h = state(i,1);
    ecc = state(i,2);
    a = (h^2)/(mu*(1-ecc^2));
    ra(i) = a + a*ecc;
    rp(i) = 2*a - ra(i);
    
    % find eclipse
    jd = jd_epoch + time(i)/86400;
    r_S = planetEphemeris(jd,'Mars','Sun');
    r_S = r_S';
    nu(i) = checkeclipse(r_temp,r_S,rm);
end
nu = nu';
% should give us a nu vect (all the eclipse indices)


% eclipse times

% Batch all eclipse indices together
for i = 1:(length(nu)-1)
    if nu(i) - nu(i+1) > 0
        eclipseDuration.index(i) = (i+1);
    elseif nu(i) - nu(i+1) < 0
        eclipseDuration.index(i) = (i);
    end % if
end % for
%%
% All eclipse enter and exit times
x = nonzeros(eclipseDuration.index);
A = zeros(length(x)-1);
for i = 1:2:(length(x)-1)
    A(i) = time(x(i+1)) - time(x(i));
end
% Time s/c is in eclipse
A = sum(nonzeros(A));
disp("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
disp("Time s/c is in eclipse: " + A + " days <------------------------") % days
disp("Percent of time s/c is in eclipse (of 50 day mission): " + 100*A/50 + " percent <------------------------")
%%
% Plot eclipse times visually (dots)
figure()
plot(time,nu,'*', 'LineWidth',1) % days
% Graph pretty 
ylim padded 
xlim tight 
xLab = xlabel('Mission time [days]','Interpreter','latex'); 
yLab = ylabel('Eclipse (1) or not (0)','Interpreter','latex'); 
plotTitle = title('Times s/c is in eclipse','interpreter','latex'); 
set(plotTitle,'FontSize',14,'FontWeight','bold') 
set(gca,'FontName','Palatino Linotype') 
set([xLab, yLab],'FontName','Palatino Linotype') 
set([xLab, yLab],'FontSize', 14) 
grid on 

% plot eclipse dots in MCI
% Pretty MARS PLOT
figure
background('Milky Way')
hold on
ax = gca;
opts.Units = 'km';
planet3D('Mars',opts)
plot3(r(:,1),r(:,2),r(:,3),'--')
plot3(r(1,1),r(1,2),r(1,3),'*','LineWidth',5)
plot3(r(end,1),r(end,2),r(end,3),'*','LineWidth',5)

% ECLIPSE LOCATION
eclipseYes = find(~nu);
p3 = plot3(r(eclipseYes,1),r(eclipseYes,2),r(eclipseYes,3),'*','Linewidth',2);
p3.Color='g';

lgd = legend("Mars","Orbital Path","Start Position","End Position","ECLIPSE",'Location','southoutside');
lgd.NumColumns = 2;
xlabel("X [Km]")
ylabel("Y [Km]")
zlabel("Z [Km]")

figure
plot(time,ra-rm,'LineWidth',2)
hold on
plot(time,rp-rm,'LineWidth',2)
grid on
legend("Apogee","Perigee",'Location','best')
ylabel("Altitude [km]")
xlabel("Time [Days]")

figure
plot(time,state(:,1)- h0)
ylabel("h-h0 [km2/s]")
xlabel("Time [days]")
grid on
figure
plot(time,state(:,2)-ecc0)
ylabel("ecc-ecc0")
xlabel("Time [days]")
grid on
figure
plot(time,rad2deg(state(:,3)-TA0))
ylabel("theta-theta0 [degs]")
xlabel("Time [days]")
grid on
figure
plot(time,rad2deg(state(:,4)-raan0))
ylabel("raan-raan0 [degs]")
xlabel("Time [days]")
grid on
figure
plot(time,state(:,5)-inc0)
ylabel("inc-inc0 [degs]")
xlabel("Time [days]")
grid on
figure
plot(time,rad2deg(state(:,6)-w0))
ylabel("w-w0 [degs]")
xlabel("Time [days]")
grid on
