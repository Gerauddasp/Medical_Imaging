clear all; close all;
rng(1);
%%PART 1
%% Question 1.a

Set1 = 1 + randn(50,1);
M1 = mean(Set1); S1 = std(Set1);
Set2 = 1.25 + randn(50,1);
M2 = mean(Set2); S2 = std(Set2);

%% Question 1.b

[h, pval, ci, stats] = ttest2(Set1,Set2);

%% Question 1.c

% i
X = [ [ones(50,1), zeros(50,1)]; [zeros(50,1), ones(50,1)] ];

%X2 = [Set1 ; Set2];
%Y = [ ones(50,1) ; zeros(50,1) ];
Y = [Set1; Set2];
% ii
checkRank = rank(X);


% iii 

beta = inv(X'*X) * X' * Y;

Yhat = X*beta;

% iv
error = Y - Yhat;

sigmasq = (norm(error))^2 / (length(Y) - 2);

% v

cos_angle = ( acos( error'* Yhat ) / ( norm(error) + norm(Yhat) ) / pi ) * 180;

% vi
covMat = sigmasq * inv(X'*X);

% vii
C = [1, -1];

% viii

t2 = (C * beta) / sqrt((C * covMat * C'));

%% Question 1.d
% adding a constant
X = [ones(length(X),1) , X];
oldX = X;
%check the rank again
CheckRank2 = rank(X);
% compute beta
beta2 = pinv(X'*X) * X' * Y;
% estimate the Y value with the beta
Yhat2 = X * beta2;
% the error vector
error2 = Y - Yhat2;
% total error
sigmasq2 = (error2' * error2) / (length(error2) - 2);
% covariance matrix
covMat2 = sigmasq2 * pinv(X'*X);
% confidence vector
C2 = [0, 1, -1];
% ttest
t22 = (C2 * beta2) ./ sqrt(C2 * covMat2 * C2');
%% Question 1.e

% We remove X2
X = X(:,1:2);
%CheckRank3 = rank(X);
% estimate new beta
beta3 = pinv(X'*X) * X' * Y;
% the new Yhat
Yhat3 = X * beta3;
% and epsilon
error3 = Y - Yhat3;
% the total error
sigmasq3 = (error3' * error3) / (length(error3) - 1);
% covariance matrix
covMat3 = sigmasq3 * pinv(X'*X);
% hypothesis
C3 = [0, 1];
% ttest
t23 = (C3 * beta3) / (C3 * covMat3 * C3');


%% Question 2.a
[h2, p2, ci2, stats2] = ttest(Set1, Set2);

%% Question 2.b
X = oldX;
X = X(:,1:2);
X = [X, [eye(50); eye(50)]]; %, [zeros(50); eye(50)]];

% i
Checkrank3 = rank(X);

% ii 
%C222 = [0, 1, 0];
C222 = [0, 1, zeros(1,50)]; % , -ones(1, 50)];
%
beta222 = pinv(X'*X) * X' * Y;
Yhat222 = X*beta222;
error222 = Y - Yhat222;
sigmasq222 = (norm(error222))^2 / (length(Y) - 51);
covMat222 = sigmasq222 * pinv(X'*X);

% iii
t222 = (C222 * beta222) ./ sqrt(C222 * covMat222 * C222');

%% Part 2 
%% Question 1.a
clear all

% question a
Set1 = randn(6,1);
Set2 = randn(8,1) * 1.25;

[h, p, ci, stats] = ttest2(Set1,Set2);

%% Question 1.b
% i
D = [Set1; Set2];

% ii 
C1 = combnk(D, 6);
for i=1:length(C1)
    C2(i,:) = setdiff(D',C1(i,:));
end

for i=1:length(C1)
    [~,~,~,STATS] = ttest2(C1(i,:), C2(i,:));
    t(i) = STATS.tstat;
end

figure;
hist(t,100)
hold on
plot(stats.tstat,0, 'x', 'Color', 'r','LineWidth',5);
title('Question 1.b');
set(gca, 'FontSize',16);
hold off


PVALUE = nnz(stats.tstat < t) / length(t);


%% Question 1.c

tstat2 = mean(Set1) - mean(Set2);

for i=1:length(C1)
    t2(i) = mean(C1(i,:)) - mean(C2(i,:));
end

figure;
hist(t2,100);
hold on
plot(stats.tstat,0, 'x', 'Color', 'r','LineWidth',5);
set(gca, 'FontSize',16);
title('Question 1.c');
hold off
PVALUE2 = nnz(tstat2 < t2) / length(t2);


%% Question 1.d

for i=1:1000
    I1 = randperm(14);
    Dn = D(I1);
    C3 = Dn(1:6);
    C4 = Dn(7:end);
    BIGCHECK(i,:)= [C3', C4'];
    t4(i) = mean(C3) - mean(C4);
    [~,~,~,STATS2] = ttest2(C3,C4);
    t3(i) = STATS2.tstat;
end
check = size(unique(BIGCHECK,'rows'),1);
duplicates = 1000 - check;
figure;
hist(t3,100)
hold on
plot(stats.tstat,0, 'x', 'Color', 'r','LineWidth',5);
set(gca, 'FontSize',16);
title('Question 1.d');
hold off

figure;
hist(t4,100);
hold on
plot(stats.tstat,0, 'x', 'Color', 'r','LineWidth',5);
set(gca, 'FontSize',16);
title('Question 1.d');
hold off
PVALUE3 = nnz(stats.tstat < t3) / 1000;
PVALUE4 = nnz(tstat2 < t4) / 1000;

 
%% Question 2
addpath('glm');
CPA = (4:11);
PPA = [3, 6, 9, 10, 13, 14, 15, 16];
n = 0;

for i=CPA
    n = n + 1;
    fid = fopen(['CPA',num2str(i),'_diffeo_fa.img'], 'r', 'l');
    datas = fread(fid, 'float'); % 16-bit floating point
    dataCPA(:,:,:,n) = reshape(datas, [40 40 40]); % dimension 40x40x40
end
n = 0;
for j=PPA
    n = n+1;
    fid = fopen(['PPA',num2str(j),'_diffeo_fa.img'], 'r', 'l');
    datas = fread(fid, 'float'); % 16-bit floating point
    dataPPA(:,:,:,n) = reshape(datas, [40 40 40]); % dimension 40x40x40
end
fid2 = fopen('wm_mask.img', 'r', 'l');
data2 = fread(fid2, 'float'); % 16-bit floating point
mask = reshape(data2, [40 40 40]); % dimension 40x40x40





%[d1, d2, d3] = ind2sub( size(dataCPA) ,find(mask ~= 0) );

I = find(mask ~= 0);

X = [ones(8,1), zeros(8,1);...
    zeros(8,1), ones(8,1)];

INVXX = pinv(X'*X) * X';
INVX = pinv(X'* X);

DATA = cat(4,dataCPA,dataPPA);
for i=1:16
    data = DATA(:,:,:,i);
    Y(i,:) = data(I);
end

beta = INVXX * Y;
Yhat = X * beta;

% calulating the error:
error = bsxfun(@minus, Y, Yhat);
%    error = P * Y;
error = bsxfun(@power, error, 2);
error = sum(error);
error = sqrt(error);
error = bsxfun(@power, error, 2);
error = error./14;
denominator = sqrt(error * 0.1250 * 2);
numerator = bsxfun(@minus, beta(1,:), beta(2,:));
t_test_a = numerator ./ denominator;

TRUEmax = max(t_test_a);



%% find all permutations
clear  Y C covMat sigmasq Yhat beta voxCPA voxPPA

subjects = (1:16);

CPAcomb = combnk(subjects,8);

for i=1:length(CPAcomb)
    PPAcomb(i,:) = setdiff(subjects, CPAcomb(i,:));
end
COMBINATIONS = [CPAcomb, PPAcomb];

tic
display('start big stuffs')

% optimise calculations before loop


I = find(mask ~= 0);
%P = eye(16) - X*(X'*X)*X';
% need to be into the loop:
tic
for it = 1:length(PPAcomb)
%     if rem(it,1000) == 0;
%         fprintf('iteration number: %i\n',it);
%         toc
%     end
    newDATA = DATA(:,:,:,COMBINATIONS(it,:));
    for i=1:16
        data = newDATA(:,:,:,i);
        Y(i,:) = data(I);
    end
    
    beta = INVXX * Y;
    Yhat = X * beta;
    
    % calulating the error:
    error = bsxfun(@minus, Y, Yhat);
%    error = P * Y;
    error = bsxfun(@power, error, 2);
    error = sum(error);
    error = sqrt(error);
    error = bsxfun(@power, error, 2);
    error = error./14;
    denominator = sqrt(error * 0.1250 * 2);
    numerator = bsxfun(@minus, beta(1,:), beta(2,:));
    t_test_b = numerator ./ denominator;
    
    MAXit(it) = max(t_test_b);
    
    
end

MULTPVALUE = nnz(MAXit > TRUEmax) / it;

% get the p-value = 0.05

ORDER = sort(MAXit);

%PVALUE05 = prctile(ORDER, 5);
PVALUE95 = prctile(ORDER, 95);

figure;
hist(MAXit,100);
hold on
plot(TRUEmax,0, 'x', 'Color', 'r','LineWidth',5);
plot(PVALUE95,0, 'o', 'Color', 'g', 'LineWidth', 5);
set(gca, 'FontSize',16);
legend('Distributions', 'Max t-stat', 'P-05');
title('Question 2');
hold off

BINIMAGE = mask;

BINIMAGE(I) = ( (t_test_a' > PVALUE95));
TTESTIMAGE(I) = t_test_a;

CC = bwconncomp(BINIMAGE, 6);
S = regionprops(CC, 'Centroid');
for o=1:length(CC.PixelIdxList)
    Recap{o,1} = num2str(S(o).Centroid);
    Recap{o,2} = length(CC.PixelIdxList{o});
    Recap{o,3} = mean(TTESTIMAGE(CC.PixelIdxList{o}));
end
