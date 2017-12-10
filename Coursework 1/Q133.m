
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


sigma = 0.04;
K = length(meas);

xBS = [0.9791, 1.3491e-09, 0.7709, 1.5968, -0.0287];
SSDBS = BallStickSSD(xBS, meas, bvals, grad_dirs);

AICBS = 2*5 + sigma^-2 * SSDBS
BICBS = 5 * log(K) + sigma^-2 * SSDBS

%%%%%%%%
xZS = [1.0051, 1.1825e-09, 0.9093, -17.2758, 0.0393, 33.7591, 0.0367];
SSDZS = BallStickSSD(xZS, meas, bvals, grad_dirs);

AICZS = 2*5 + sigma^-2 * SSDZS
BICZS = 5 * log(K) + sigma^-2 * SSDZS
%%%%%%%%%%%%%%
xZST = [1.0087, 1.1692e-09, 0.9080, -1.5425, 0.0105, 486.9675];
SSDZST = BallStickSSD(xZST, meas, bvals, grad_dirs);

AICZST = 2*5 + sigma^-2 * SSDZST
BICZST = 5 * log(K) + sigma^-2 * SSDZST

