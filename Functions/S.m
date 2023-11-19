function [output] = S(alpha,X)
% S(z) Stumpff Function
% Defined by infinte series

z = alpha*(X^2);
output = (1/6)-(z/120)+((z^2)/5040)-((z^3)/362880);

end