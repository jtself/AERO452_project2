function [dstate] = vop_ODE(time, state, wEarth, re, mu, muSun, Cd, area, mass, jd_epoch)

%% Find COES
h = state(1);
ecc = state(2);
theta = state(3);
raan = state(4);
inc = state(5);
omega = state(6);

%% COES2RandV

% Solve cos and sin
u = omega + theta;
s_raan = sin(raan);
c_rran = cos(raan);
s_u = sin(u);
c_u = cos(u);
s_inc = sin(inc);
c_inc = cos(inc);
s_w = sin(omega);
c_w = cos(omega);
c_theta = cos(theta);
s_theta = sin(theta);

% From COES into perifocal frame
rperifocal = (h^2/mu) * (1/(1 + ecc * c_theta)) *  [c_theta; s_theta; 0]; % p, q, w coordinate vects
vperifocal = (mu/h) * [-s_theta; ecc + c_theta; 0];                      % p, q, w coordinate vects
C_ECI_PERI = [c_w*c_rran - c_inc*s_w*s_raan,c_w*s_raan + c_inc*c_rran*s_w,s_inc*s_w;
    - c_rran*s_w - c_inc*c_w *s_raan, c_inc*c_w*c_rran - s_w*s_raan, c_w*s_inc;
    s_inc*s_raan, -c_rran*s_inc,c_inc];

% Get r and v in ECI
rVect = C_ECI_PERI' * rperifocal;
vVect = C_ECI_PERI' * vperifocal;
r = norm(rVect);

%% Build Qxx Rotaion Matrix (ECI --> LVLH)
R_hat = rVect/r;
N_hat = cross(rVect,vVect)/norm(cross(rVect,vVect));
T_hat = cross(N_hat,R_hat);
Qxr = [R_hat T_hat N_hat];
Qxr = Qxr';

%% Find Drag Pertubation
alt = r-6378;
rho = (1e9)*expDrag(alt);
Vrel = vVect - cross(wEarth,rVect);
acc_drag = -0.5*Cd*area*(1/mass).*rho.*(norm(Vrel)).*(Vrel);
acc_drag = Qxr*acc_drag;

%% Find Solar Gravity Pertubation
jd = jd_epoch + time/(60*60*24);
[~, ~, r_S] = solar_position(jd);
r_S_SC = r_S - rVect;
q = dot(rVect,(2*r_S - rVect))/(norm(r_S)^2);
F = q*(q^2 - 3*q + 3)/(1+ ((1-q)^(1.5)))';
acc_SolarGravity = muSun*(F*r_S - rVect)/(norm(r_S_SC)^3);
acc_SolarGravity = Qxr*acc_SolarGravity;

%% Find J2-J6 Pertubation

J2 = 1.08262668355e-3;
J3 =  -2.53265648533e-6;
J4 = -1.61962159137e-6;
J5 = -2.27296082869e-7;
J6 = 5.40681239107e-7;

aI_J2 = ((-3*J2*mu*(re^2)*rVect(1))/(2*(r^5)))*(1 - (5*(rVect(3)^2))/(r^2));
aJ_J2 = ((-3*J2*mu*(re^2)*rVect(2))/(2*(r^5)))*(1 - (5*(rVect(3)^2))/(r^2));
aK_J2 = ((-3*J2*mu*(re^2)*rVect(3))/(2*(r^5)))*(3 - (5*(rVect(3)^2))/(r^2));

aI_J3 = ((-5*J3*mu*(re^3)*rVect(1))/(2*(r^7)))*(3*rVect(3) - (7*(rVect(3)^3))/(r^2));
aJ_J3 = ((-5*J3*mu*(re^3)*rVect(2))/(2*(r^7)))*(3*rVect(3) - (7*(rVect(3)^3))/(r^2));
aK_J3 = ((-5*J3*mu*(re^3))/(2*(r^7)))*(6*(rVect(3)^2) - (7*(rVect(3)^4))/(r^2) - (3/5)*r^2);

aI_J4 = ((15*J4*mu*(re^4)*rVect(1))/(8*(r^7)))*(1- ((14*(rVect(3)^2))/(r^2)) + ((21*(rVect(3)^4))/(r^4)));
aJ_J4 = ((15*J4*mu*(re^4)*rVect(2))/(8*(r^7)))*(1- ((14*(rVect(3)^2))/(r^2)) + ((21*(rVect(3)^4))/(r^4)));
aK_J4 = ((15*J4*mu*(re^4)*rVect(3))/(8*(r^7)))*(5- ((70*(rVect(3)^2))/(3*(r^2))) + ((21*(rVect(3)^4))/(r^4)));

aI_J5 = ((3*J5*mu*(re^5)*rVect(1)*rVect(3))/(8*(r^9)))*(35- 210*((rVect(3)^2)/(r^2))+ 231*((rVect(3)^4)/(r^4)));
aJ_J5 = ((3*J5*mu*(re^5)*rVect(2)*rVect(3))/(8*(r^9)))*(35- 210*((rVect(3)^2)/(r^2))+ 231*((rVect(3)^4)/(r^4)));
aK_J5 = ((3*J5*mu*(re^5)*(rVect(3)^2))/(8*(r^9)))*(105- 315*((rVect(3)^2)/(r^2))+ 231*((rVect(3)^4)/(r^4))) - (15*J5*mu*(re^5))/(8*(r^7));

aI_J6 = ((-J6*mu*(re^6)*rVect(1))/(16*(r^9)))*(35 - 945*((rVect(3)^2)/(r^2)) + 3465*((rVect(3)^4)/(r^4)) -3003*((rVect(3)^6)/(r^6)));
aJ_J6 = ((-J6*mu*(re^6)*rVect(2))/(16*(r^9)))*(35 - 945*((rVect(3)^2)/(r^2)) + 3465*((rVect(3)^4)/(r^4)) -3003*((rVect(3)^6)/(r^6)));
aK_J6 = ((-J6*mu*(re^6)*rVect(3))/(16*(r^9)))*(245 - 2205*((rVect(3)^2)/(r^2)) + 4851*((rVect(3)^4)/(r^4)) -3003*((rVect(3)^6)/(r^6)));

aI = aI_J2 + aI_J3 + aI_J4 + aI_J5 + aI_J6;
aJ = aJ_J2 + aJ_J3 + aJ_J4 + aJ_J5 + aJ_J6;
aK = aK_J2 + aK_J3 + aK_J4 + aK_J5 + aK_J6;

acc_J2_6 = [aI; aJ; aK];
acc_J2_6 = Qxr*acc_J2_6;

%% Sum All Pertubations

p = acc_drag + acc_SolarGravity + acc_J2_6; %  + acc_SRP 
pr = p(1);
ps = p(2);
pw = p(3);

%% Solve VoP Gauss Formulation

dh = r*ps;

decc = (h/mu)*s_theta*pr + (1/(mu*h))*((h^2 + mu*r)*c_theta + mu*ecc*r)*ps;

% theta
twobodymotion= h/r^2;
dtheta_pert = (1/(ecc*h))*((h^2/mu)*c_theta*pr - (r + (h^2/mu))*s_theta*ps);
dtheta =  twobodymotion + dtheta_pert;

draan = (r/(h*s_inc))*s_u*pw;

dinc = (r/h)*c_u*pw;

domega = -dtheta_pert- ((r*s_u)/(h*tan(inc)))*pw;

%% Build dState
dstate = [dh; decc; dtheta; draan; dinc; domega];

end