function [output] = C(alpha,X)
% C(z) Stumpff Function
% Defined by infinte series

z = alpha*(X^2);
output = (1/2)-(z/24)+((z^2)/720)-((z^3)/40320);

end