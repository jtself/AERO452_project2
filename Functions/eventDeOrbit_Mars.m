function [value,isterminal,direction] = eventDeOrbit_Mars(time, state, wMars, rm, mu, muSun, Cd, area, mass, jd_epoch, Cr, Psr)

% Find COES
h = state(1);
ecc = state(2);
theta = state(3);
raan = state(4);
inc = state(5);
omega = state(6);

% Solve cos and sin
s_raan = sin(raan);
c_rran = cos(raan);
s_inc = sin(inc);
c_inc = cos(inc);
s_w = sin(omega);
c_w = cos(omega);
c_theta = cos(theta);
s_theta = sin(theta);

% From COES into perifocal frame
rperifocal = (h^2/mu) * (1 + ecc * c_theta)^-1 *  [c_theta; s_theta; 0]; % p, q, w coordinate vects
C_ECI_PERI = [c_w*c_rran - c_inc*s_w*s_raan,c_w*s_raan + c_inc*c_rran*s_w,s_inc*s_w;
    - c_rran*s_w - c_inc*c_w *s_raan, c_inc*c_w*c_rran - s_w*s_raan, c_w*s_inc;
    s_inc*s_raan, -c_rran*s_inc,c_inc];

% Get r and v in ECI
rVect = C_ECI_PERI' * rperifocal;
r = norm(rVect);

ALT = r - rm;

value = ALT - 125;
isterminal = 1;
direction = -1;

end