% Question 2
clear all; close all;

fid=fopen('dwi.Bfloat', 'r', 'b');
dwis = fread(fid, 'float');
fclose(fid);


dwis = reshape(dwis, 33, 112, 112, 50);

qhat = load('grad_dirs.txt')';

bvals = 1000*sum(qhat.*qhat);

Avox = dwis(:,52,62,25);

% Define a starting point for the non-linear fit

startx = [250000, 1E-3, 0.5, 0, 0];
startx(3) = asin(sqrt(startx(3)));

h=optimset('MaxFunEvals',20000,...
   'Algorithm','levenberg-marquardt',...
   'LargeScale','off',...
   'Display','off',...
   'TolX',1e-1000,...
   'TolFun',1e-1000);

% Now run the fitting
[parameter_hat,RESNORM,EXITFLAG,OUTPUT]=fminunc('BallStickSSD',startx,h,Avox,bvals,qhat);
parameter_hat(1) = abs(parameter_hat(1));
parameter_hat(2) = abs(parameter_hat(2));
parameter_hat(3) = sin(parameter_hat(3))^2;

% Here we try to compare the results we obtain:
x = parameter_hat;
S0 = x(1);
diff = x(2);
f = x(3);
theta = x(4);
phi = x(5);

% Synthesize the signals
fibdir = [cos(phi)*sin(theta) sin(phi)*sin(theta) cos(theta)];
fibdotgrad = sum(qhat.*repmat(fibdir, [length(qhat), 1])');
model = S0*(f*exp(-bvals*diff.*(fibdotgrad.^2)) + (1-f)*exp(-bvals*diff));

figure;
% Plot the actual data points!
plot(Avox, ' bs', 'MarkerSize', 16, 'LineWidth', 4);
hold on;
% Add the predictions to the plot!
plot(model, ' rx', 'MarkerSize', 16, 'LineWidth', 4);
% Add labels and legend.!
set(gca, 'FontName', 'Times');
set(gca, 'FontSize', 26);
xlabel('\bf{q} index');
ylabel('S');
legend('Data', 'Model');