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

q = grad_dirs';
b = bvals';

y = [ones(size(b)), - b .* q(:,1).^2, -2 * b .* q(:,1) .* q(:,2), - 2 * b .* ...
    q(:,1) .* q(:,3), - b .* q(:,2).^2, - 2 * b .* q(:,2) .* q(:,3), - b .* q(:,3).^2]; 

x = y\log(meas);

S0init = exp(x(1));
Dinit = mean([x(2), x(5), x(7)]);
startx(1) = S0init;
startx(2) = Dinit;

D = [ x(2), x(3), x(4);...
    x(3), x(5), x(6);...
    x(4), x(6), x(7)];

lambda = eig(D);
FA = (3/2)*sqrt((sum((lambda - startx(2)).^2)/sum(lambda)));
