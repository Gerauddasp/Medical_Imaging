clear all; close all;
ws = load('wsQ114');
Parameters = ws.Parameters;
AllResnorm = ws.AllResnorm;

%%%%%% First we represent the data as they are.
figure;
subplot(2,2,1)
S0im = squeeze(Parameters(1,:,:));
immin1 = min(min(S0im));
immax1 = max(max(S0im));
image1 = mat2gray(S0im,[immin1, immax1]);
imshow(image1,[]);
title('S0');

subplot(2,2,2)
dim = squeeze(Parameters(2,:,:));
immin2 = min(min(dim));
immax2 = max(max(dim));
image2 = mat2gray(dim,[immin2, immax2]);
imshow(image2,[]);
title('diffusion');

subplot(2,2,3)
fim = squeeze(Parameters(3,:,:));
immin3 = min(min(fim));
immax3 = max(max(fim));
image3 = mat2gray(fim,[immin3, immax3]);
imshow(image3,[]);
title('f');

subplot(2,2,4)
immin4 = min(min(AllResnorm));
immax4 = max(max(AllResnorm));
image4 = mat2gray(AllResnorm, [immin4, immax4]);
imshow(image4,[]);
title('SSD');

%%%%%%%%%%%%%%%%%%%%%%
% Plot histograms

figure;
subplot(2,2,1);
hist(S0im(:),100);

subplot(2,2,2);
hist(dim(:),100);

subplot(2,2,3);
hist(fim(:),100);

subplot(2,2,4);
hist(AllResnorm(:),100);


%%%%%%%%%%%%%%%%%%%%%%


figure;
subplot(2,2,1)
S0im = squeeze(Parameters(1,:,:));
immin1 = prctile(S0im(:), 5);
immax1 = prctile(S0im(:), 95);
image1 = mat2gray(S0im,[immin1, immax1]);
imshow(image1,[]);
title('S0');

subplot(2,2,2)
dim = squeeze(Parameters(2,:,:));
immin2 = prctile(dim(:), 0);
immax2 = prctile(dim(:), 67);
image2 = mat2gray(dim,[immin2, immax2]);
imshow(image2,[]);
title('diffusion');

subplot(2,2,3)
fim = squeeze(Parameters(3,:,:));
immin3 = prctile(fim(:), 5);
immax3 = prctile(fim(:), 95);
image3 = mat2gray(fim,[immin3, immax3]);
imshow(image3,[]);
title('f');


subplot(2,2,4)
immin4 = prctile(AllResnorm(:), 5);
immax4 = prctile(AllResnorm(:), 95);
image4 = mat2gray(AllResnorm, [immin4, immax4]);
imshow(image4,[]);
title('SSD');

figure;
Phi = squeeze(Parameters(4,:,:));
Theta = squeeze(Parameters(5,:,:));
row = (1:112);
col = (1:112);
u = (cos(Phi).* sin(Theta)) .* squeeze(Parameters(3,:,:));
v = (sin(Phi).* sin(Theta)) .* squeeze(Parameters(3,:,:));
title('fibre direction');

quiver(row, col, u, v);

