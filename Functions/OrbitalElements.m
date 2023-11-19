function [h, inc, RAAN, ecc, w, theta, epsilon, a, T] = OrbitalElements(r0Vect,v0Vect,mu)
% Finds the Six Classical Orbital Elements

r = norm(r0Vect);
v = norm(v0Vect);
vr = dot(r0Vect,v0Vect)/r;


% h, specific angular momentum
hVect = cross(r0Vect,v0Vect);
h = norm(hVect);

% inc, inclination
inc = acos(hVect(3)/h);

while inc > pi
    inc = inc - pi;
end
while inc < 0 
    inc = inc + pi;
end

inc = inc*(180/pi);

% omega, right ascension of ascending node (RAAN)
NVect = cross([0;0;1],hVect);
N = norm(NVect);

if NVect(2) < 0
    RAAN = 2*pi - acos(NVect(1)/N);
else
    RAAN = acos(NVect(1)/N);
end

RAAN = RAAN*(180/pi);


% ecc, eccentricity
eccVect = (1/mu).*((v^2-(mu/r)).*r0Vect-r*vr.*v0Vect);
ecc = norm(eccVect);

% w, argument of perigee

if eccVect(3) < 0
    w = 2*pi - acos(dot(NVect,eccVect)/(N*ecc));
else
    w = acos(dot(NVect,eccVect)/(N*ecc));
end

w = w*(180/pi);

% theta, true anomaly
if vr >= 0
    theta = acos(dot(eccVect,r0Vect)/(ecc*r));
else
    theta = 2*pi - acos(dot(eccVect,r0Vect)/(ecc*r));
end

theta = theta*(180/pi);

% epsilon, specific energy
epsilon = (0.5*v^2) - (mu/r);
if epsilon > 0
    epsilon = nan;
end

% a, semi-major axis
a = -mu/(2*epsilon);


% T, period
T = 2*pi*sqrt((a^3)/mu);



end