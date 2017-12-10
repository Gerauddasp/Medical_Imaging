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

x = [0.9791, 1.3491e-09, 0.7709, 1.5968, -0.0287];

%%%%%%%% FisherMat = funciton(x meas, bvals, grad_dirs)
bd = exp(-2* bvals' * x(2));

value1 = sum(x(1) * bd);
value2 = sum (x(1)* -bvals'.* bd);
value3 = sum (x(1)* -bvals'.* bd);
value4 = sum (x(1)^2* -(bvals').^2.* bd);

FisherMat = [value1, value2;...
    value3, value4];


