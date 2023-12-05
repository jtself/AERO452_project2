function [dstate] = vop_mars_ode(time, state, wMars, rm, mu, muSun, Cd, area, mass, jd_epoch, Cr, Psr)

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

%% Find Current JD and Solar Position

jd = jd_epoch + time/(60*60*24);
r_S = planetEphemeris(jd,'Mars','Sun');
r_S = r_S';


%% Build Qxx Rotaion Matrix (ECI --> LVLH)
R_hat = rVect/r;
N_hat = cross(rVect,vVect)/norm(cross(rVect,vVect));
T_hat = cross(N_hat,R_hat);
Qxr = [R_hat T_hat N_hat];
Qxr = Qxr';

%% Find Drag Pertubation
alt = r-rm;
rho = (1e9)*marsDrag(alt);
Vrel = vVect - cross(wMars,rVect);
acc_drag = -0.5*Cd*area*(1/mass).*rho.*(norm(Vrel)).*(Vrel);
acc_drag = Qxr*acc_drag;

%% Find Solar Gravity Pertubation
r_S_SC = r_S - rVect;
q = dot(rVect,(2*r_S - rVect))/(norm(r_S)^2);
F = q*(q^2 - 3*q + 3)/(1+ ((1-q)^(1.5)))';
acc_SolarGravity = muSun*(F*r_S - rVect)/(norm(r_S_SC)^3);
acc_SolarGravity = Qxr*acc_SolarGravity;

%% Find SRP Pertubation

thetaB = acos(rm/r);
thetaA = acos(rm/r_S);
theta = acos(dot(r_S,rVect)/(norm(r_S)*r));

if (thetaA + thetaB) < theta
    F = 0;
else
    F = 1;
end

acc_SRP = -(1000)*Psr*Cr*area*(1/mass)*(r_S)*(1/norm(r_S))*F;
acc_SRP = Qxr*acc_SRP;


%% Find J2-J3 Pertubation

J2 = 1.9555e-3;
J3 =  3.14498e-5;

aI_J2 = ((-3*J2*mu*(rm^2)*rVect(1))/(2*(r^5)))*(1 - (5*(rVect(3)^2))/(r^2));
aJ_J2 = ((-3*J2*mu*(rm^2)*rVect(2))/(2*(r^5)))*(1 - (5*(rVect(3)^2))/(r^2));
aK_J2 = ((-3*J2*mu*(rm^2)*rVect(3))/(2*(r^5)))*(3 - (5*(rVect(3)^2))/(r^2));

rx = rVect(1);
ry = rVect(2);
rz = rVect(3);
aI_J3 = ( (-5*J3*mu*rm^3*rx) / (2*r^7) ) * (3*rz - (7*rz^3)/(r^2) );
aJ_J3= ( (-5*J3*mu*rm^3*ry) / (2*r^7) ) * (3*rz - (7*rz^3)/(r^2) );
aK_J3 = ( (-5*J3*mu*rm^3) / (2*r^7) ) * (6*rz^2 - (7*rz^4)/(r^2) - (3/5)*r^2 );

aJ2 = [aI_J2;aJ_J2;aK_J2];
aJ2 = Qxr*aJ2;

aJ3 = [aI_J3;aJ_J3;aK_J3];
aJ3 = Qxr*aJ3;

acc_J2_3 = aJ2 + aJ3;

%% Sum All Pertubations

p = acc_drag + acc_SolarGravity + acc_J2_3 + acc_SRP;
pr = p(1);
ps = p(2);
pw = p(3);

%% Solve VoP Gauss Formulation

dh = r*ps;

decc = (h/mu)*s_theta*pr + (1/(mu*h))*((h^2 + mu*r)*c_theta + mu*ecc*r)*ps;

twobodymotion= h/r^2;
dtheta_pert = (1/(ecc*h))*((h^2/mu)*c_theta*pr - (r + (h^2/mu))*s_theta*ps);
dtheta =  twobodymotion + dtheta_pert;

draan = (r/(h*s_inc))*s_u*pw;

dinc = (r/h)*c_u*pw;

domega = -dtheta_pert- ((r*s_u)/(h*tan(inc)))*pw;

%% Build dState
dstate = [dh; decc; dtheta; draan; dinc; domega];

end