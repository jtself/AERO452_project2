function [v1vec, v2vec] = lambert(r1vec, r2vec, dt, tm, mu)
% Lamberts: Universial Variable Solution
% r1: initial position vector [km]
% r2: final position vector [km]
% t: transfer time [s]
% tm: Direction between the points (1 for short way, -1 for long)

% assume prograde
% dt in sec

r1 = norm(r1vec);
r2 = norm(r2vec);
C21 = cross(r1vec,r2vec);
delta_theta = acos( dot(r1vec,r2vec) / (r1*r2) );


if tm == 1
    if C21(3) < 0
        delta_theta = 2*pi - delta_theta;
    end
else
    if C21(3) >=0
        delta_theta = 2*pi - delta_theta;
    end
end


% find A
A = ((sqrt(r1*r2))*sin(delta_theta))/(sqrt(1-cos(delta_theta)));

% y(z)
y = @(z) r1 + r2 + A*( (z*S(z)-1) / (sqrt(C(z)) ) );

% Chi X
X = @(z) sqrt( y(z)/C(z) );

% dtloop

dtloop = @(z) (((X(z)^3)*S(z))/sqrt(mu)) + (A*sqrt(y(z))/sqrt(mu));

% z
z = 0;
zL = -4*pi^2;
zU = 4*pi^2;

TOL = 10^(-5);
while abs(dtloop(z) - dt) > TOL
    
    if dtloop(z) <= dt
        zL = z;
    else
        zU = z;
    end

z = (zU+zL)/2;
end

% lagrange coefficients
f = @(z) 1 - (y(z)/r1);

g = @(z) A*sqrt(y(z)/mu);

fdot = @(z) (sqrt(mu)/(r1*r2)) * sqrt( (y(z))/(C(z)) ) * (z*S(z)-1);

gdot = @(z) 1 - (y(z)/r2);

v1vec = (1/g(z))*(r2vec - f(z)*r1vec);
v2vec = (1/g(z))*(gdot(z)*r2vec-r1vec);

% nested functions
    % stumpff
    function S = S(z)
    % for S
        if z > 0
            S = (sqrt(z) - sin(sqrt(z))) / ((sqrt(z))^3);
        elseif z < 0
            S = ( sinh(sqrt(-z)) - sqrt(-z) ) / ((sqrt(-z))^3);
        else 
            S = (1/6);
        end
    end
    
    function C = C(z)
           % for C
        if z > 0
            C = (1-cos(sqrt(z)))/z;
        elseif z < 0
            C = ( (cosh(sqrt(-z))) - 1)/(-z);
        else
            C = (1/2);
        end
    end

    
end
