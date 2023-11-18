function [r,v] = COES2RandV(h,ecc,inc,raan,omega,theta,mu)
% Function to find r and v from given COES

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Provide h in km2/s
% Provide inc, raan, omega, and theta in rad
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% From COES into perifocal frame
rperifocal = (h^2/mu) * (1 + ecc * cos(theta))^-1 *  [cos(theta); sin(theta); 0]; % p, q, w coordinate vects
vperifocal = (mu/h) * [-sin(theta); ecc + cos(theta); 0];                         % p, q, w coordinate vects

% Find rotation matrix to go from PERI to ECI
Cz_omega = [cos(omega)  sin(omega) 0;
       -sin(omega)  cos(omega) 0
        0       0      1];

Cx_inc  = [1   0           0
           0   cos(inc)  sin(inc) 
           0  -sin(inc)  cos(inc)];

Cz_RAAN = [cos(raan) sin(raan) 0;
          -sin(raan) cos(raan) 0
           0         0         1];

C_ECI_PERI = Cz_omega*Cx_inc*Cz_RAAN;

% Get r and v in ECI
r = C_ECI_PERI' * rperifocal;
v = C_ECI_PERI' * vperifocal;

end