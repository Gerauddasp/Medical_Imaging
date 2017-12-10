% here we only have to compute p(AIx):


close all; clear all;

datafile = 'challengeOPEN.txt'; fid = fopen(datafile, 'r', 'b');
% Read in the header
A = fgetl(fid);
% Read in the data
A = fscanf(fid, '%f', [8, inf]); fclose(fid);
% Create the protocol
meas = A(1,:)'; 
grad_dirs = A(2:4,:); 
G = A(5,:);
delta = A(6,:); 
smalldel = A(7,:);
TE = A(8,:);
GAMMA = 2.675987E8;
bvals = ((GAMMA*smalldel.*G).^2).*(delta-smalldel/3);

Avox = meas;
qhat = grad_dirs;

% first model:

sigma = 0.04;

x = [0.9791, 1.3491e-09, 0.7709, 1.5968, -0.0287];    
S0 = abs(x(1));
diff = abs(x(2));
f = sin(x(3))^2;
theta = x(4);
phi = x(5);


% Synthesize the signals
fibdir = [cos(phi)*sin(theta) sin(phi)*sin(theta) cos(theta)];
fibdotgrad = sum(qhat.*repmat(fibdir, [length(qhat), 1])');
S1 = S0*(f*exp(-bvals*diff.*(fibdotgrad.^2)) + (1-f)*exp(-bvals*diff));
%sigma = mean(Avox-S1');
pBSx = prod(-log((2*pi*sigma^2)) - (-((Avox-S1').^2/(2*sigma)^2)));

% Now for the second model:
x = [1.0051, 1.1825e-09, 0.9093, -17.2758, 0.0393, 33.7591, 0.0367];
S0 = abs(x(1));
diff = abs(x(2));
f = sin(x(3))^2;
theta = x(4);
phi = x(5);
lambda1 = abs(x(6));
lambda2 = abs(x(7));

if lambda1 < lambda2
    temp = lambda1;
    lambda1 = lambda2;
    lambda2 = temp;
    clear('temp');
end

% Synthesize the signals
fibdir = [cos(phi)*sin(theta) sin(phi)*sin(theta) cos(theta)];
fibdotgrad = sum(qhat.*repmat(fibdir, [length(qhat), 1])');


S2 = S0*(f*exp(-bvals*diff.*(fibdotgrad.^2)) + (1-f)*exp(-bvals.*(lambda2 + (lambda1 - lambda2)*(fibdotgrad.^2))));

%sigma = mean(Avox-S2');
pZSx = prod((2*pi*sigma^2)^-1 * exp(-((Avox-S2').^2/(2*sigma)^2)))

% And for the last one
x = [1.0087, 1.1692e-09, 0.9080, -1.5425, 0.0105, 486.9675];
S0 = abs(x(1));
diff = abs(x(2));
f = sin(x(3))^2;
theta = x(4);
phi = x(5);
lambda1 = abs(x(6));
lambda2 = (1-f) * lambda1;


% Synthesize the signals
fibdir = [cos(phi)*sin(theta) sin(phi)*sin(theta) cos(theta)];
fibdotgrad = sum(qhat.*repmat(fibdir, [length(qhat), 1])');

S3 = S0*(f*exp(-bvals*diff.*(fibdotgrad.^2)) + (1-f)*exp(-bvals.*(lambda2 + (lambda1 - lambda2)*(fibdotgrad.^2))));

%sigma = mean(Avox-S3');
pZSTx = prod((2*pi*sigma^2)^-1 * exp(-((Avox-S3').^2/(2*sigma)^2)))