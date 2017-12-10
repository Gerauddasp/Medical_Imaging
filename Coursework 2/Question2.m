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

[d1, d2, d3] = ind2sub( size(dataCPA) ,find(mask ~= 0) );

tic
for v=1:length(d1)
    
    voxCPA = squeeze( dataCPA( d1(v), d2(v), d3(v), : ) );
    voxPPA = squeeze( dataPPA( d1(v), d2(v), d3(v), : ) );
    % calculate the GLM
    Y = [voxCPA; voxPPA];
    X = [ones(length(voxCPA),1), zeros(length(voxCPA),1);...
         zeros(length(voxPPA),1), ones(length(voxPPA),1)];
    beta = pinv(X'*X) * X' * Y;
    Yhat = X * beta;
    sigmasq = (norm(Y - Yhat)^2) / (length(Y) - 2);
    covMat = sigmasq * pinv(X'*X);
    C = [1, -1];
    t_test_a(d1(v), d2(v), d3(v)) = (C * beta) / sqrt((C * covMat * C'));
    
end
TRUEmax = max(max(max(t_test_a)));

if true
% find all permutations


subjects = (1:16);

CPAcomb = combnk(subjects,8);

for i=1:length(CPAcomb)
    PPAcomb(i,:) = setdiff(subjects, CPAcomb(i,:));
end

DATA = cat(4,dataCPA,dataPPA);
tic
display('start big stuffs')
for it=1:length(CPAcomb)
    if rem(it,100) == 0
        fprintf('iteration: %i\n', it);
        toc
    end
    newdataCPA = DATA(:,:,:,CPAcomb(it,:));
    newdataPPA = DATA(:,:,:,PPAcomb(it,:));
    for v=1:length(d1)
        voxCPA = squeeze( newdataCPA( d1(v), d2(v), d3(v), : ) );
        voxPPA = squeeze( newdataPPA( d1(v), d2(v), d3(v), : ) );
        % calculate the GLM
        Y = [voxCPA; voxPPA];
        X = [ones(length(voxCPA),1), zeros(length(voxCPA),1);...
            zeros(length(voxPPA),1), ones(length(voxPPA),1)];
        beta = pinv(X'*X) * X' * Y;
        Yhat = X * beta;
        sigmasq = (norm(Y - Yhat)^2) / (length(Y) - 2);
        covMat = sigmasq * pinv(X'*X);
        C = [1, -1];
        t_test_b(d1(v), d2(v), d3(v)) = (C * beta) / sqrt((C * covMat * C'));
    end
    MAXit(it) = max(max(max(t_test_b)));
end
save('VoxStats.mat','MAXit');
else
    load('VoxStats.mat');

end
MULTPVALUE = nnz(MAXit > TRUEmax) / it;


display('total time')
toc