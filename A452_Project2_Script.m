%{
AERO452 | SPACEFLIGHT DYNAMICS II
Group Project #2
Authors: 
    Travis Bouck
    Justin Self
Due date: 
    Dec. 8, 2023
%}

% Housekeeping
clear all; close all; clc;

addpath("Functions\")

%% Section 1

% Governing Constants
rearth = 6378; % km
mu = 398600; % km3/s2

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TLE Input (RANDOM TLE- GEO from other Proj)
SC.init.UTC = "02 October 2023, 05:38:01"; % UTC
SC.init.ecc = 0.000001016; % actual ecc but we must change it to zero
SC.init.inc = deg2rad(0.0009);               % rad
SC.init.rp = 35782 + rearth;                 % km
SC.init.ra = 35791 + rearth;                 % km
SC.init.raan = deg2rad(246.2601);            % rad
SC.init.w = deg2rad(273.8346);               % rad
SC.init.revperday = 1.00269599;              % rev/day
SC.init.Me = deg2rad(199.9092);              % rad

% Calculate Additional COEs
SC.init.a = 0.5*(SC.init.rp + SC.init.ra);             % km
SC.init.T = ( (2*pi) / (sqrt(mu)) ) * SC.init.a^(3/2); % sec
SC.init.n = SC.init.revperday*(2*pi/(24*60*60));       % rad/sec

SC.init.E = newtonsKepler(SC.init.Me, SC.init.ecc);
SC.init.TA = 2*atan((tan(SC.init.E/2))...
    /(sqrt((1 - SC.init.ecc)/(1 + SC.init.ecc))));
SC.init.h = findh(SC.init.a, mu, SC.init.ecc, SC.init.TA);

% Determine epoch r and v
[SC.init.rVect,SC.init.vVect] = COES2RandV(SC.init.h,SC.init.ecc,...
    SC.init.inc,SC.init.raan,SC.init.w,SC.init.TA,mu);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%






%% Section 2