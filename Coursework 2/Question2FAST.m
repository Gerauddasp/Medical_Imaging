clear all; close all;
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
    if rem(it,1000) == 0;
        fprintf('iteration number: %i\n',it);
        toc
    end
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




