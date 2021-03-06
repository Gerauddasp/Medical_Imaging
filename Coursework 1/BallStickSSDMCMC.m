function [sumRes, resJ] = BallStickSSDMCMC(x, Avox, bvals, qhat)
% Extract the parameters
%x = real(x);
S0 = abs(x(1));
diff = abs(x(2));
f = sin(x(3))^2;
theta = x(4);
phi = x(5);


% Synthesize the signals
fibdir = [cos(phi)*sin(theta) sin(phi)*sin(theta) cos(theta)];
fibdotgrad = sum(qhat.*repmat(fibdir, [length(qhat), 1])');
S = S0*(f*exp(-bvals*diff.*(fibdotgrad.^2)) + (1-f)*exp(-bvals*diff));

% Compute the sum of square differences
sumRes = sum((Avox - S').^2);
resJ = S';

end