% Q5.3
%%% 

clear all; 
close all;
subject = 5;
PLD = [50, 50, 50, 50, 300, 700, 1200];
TD = [500, 750, 1000, 3000, 3000, 3000, 3000];
echo = [29, 41, 53, 68.2, 80.2, 92.2, 107.4, 119.4, 131.4, 146.6, 158.6, 170.6]; 
% Import data from question Q1.2

%load('T2vsTE.mat');
load('data');

%%
% Here we decide wich subject we are going to choose and also the
% parameters for fminunc
subject = 1;
h=optimset('MaxFunEvals',200000,...
    'Algorithm','levenberg-marquardt',...
    'LargeScale','off',...
    'Display','off',...
    'TolX',1e-5,...
    'TolFun',1e-5);
% Now these are the starting points for S0 (startx(1)), and T2 (startx(2))
startx(1) = 1;
startx(2) = 120;

% We create two new figures to plot our results:
map1 = figure;
map2 = figure;

% We cut the parts outside of the brain, put the value to true to re-cut,
% we do it once for every subjects.
if false
    for i=1:9
        figure;
        mask(i) = struct('BW',[], 'xi', [],'yi', []);
        immax1 = max(max(data(i).new_images_off(:,:,1,1)));
        immin1 = min(min(data(i).new_images_off(:,:,1,1)));
        image1 = mat2gray(data(i).new_images_off(:,:,1,1),[immin1, immax1]);
        [mask(i).BW, mask(i).xi, mask(i).yi] = roipoly(image1);
    end
    save('braincut.mat','mask');
% Or loading previous cut
else
    load('braincut.mat');
end

% deleting data that is outside of the cut, we put the data to zero and
% will excluse it after from our fitting procedure
for subject = 1:9
    data(subject).new_images_off = data(subject).new_images_off .* repmat(mask(subject).BW,1,1,12,14);
    data(subject).new_images_on = data(subject).new_images_on .* repmat(mask(subject).BW,1,1,12,14);
    data(subject).perf_w_off = data(subject).perf_w_off .* repmat(mask(subject).BW,1,1,12,7);
    data(subject).perf_w_on = data(subject).perf_w_on .* repmat(mask(subject).BW,1,1,12,7);
end


% Now we are calclating the T2 value in each voxel for every configuration
% with gradient on and gradient off.
computevox = false;
if not(computevox)
    load('voxMap.mat');
end
pd = makedist('Normal');
t = truncate(pd,0.8,1.2);
for control=1:7
    name=strcat('tau: ', num2str(TD(control)),'   PLD: ', num2str(PLD(control)));
    if computevox
        fprintf('Control: %i/7\n', control)
        Maps(control) = struct('GradOFF', [], 'GradON', []);
        for lin=1:size(data(subject).perf_w_off,1)
            fprintf('line number: %i\n', lin);
            for col=1:size(data(subject).perf_w_off,2)
            %testing if we are in our outside the brain with the signal:
            if mean(squeeze(data(subject).perf_w_off(lin,col,:,control))) == 0
                Maps(control).GradOFF(lin,col) = 0;
                Maps(control).GradON(lin,col) = 0;
            else
                for i =1:10
                    r = random(t,1,2);
                    startxn = startx .* r;
                    [x1, fval, exitflag, output] = fminunc('OneExp', startxn, h,...
                        squeeze(data(subject).perf_w_off(lin,col,:,control)),echo');
                    res1(i,1:2) = abs(x1);
                    res1(i,3) = OneExp(abs(x1), squeeze(data(subject).perf_w_off(lin,col,:,...
                        control)),echo');
                    
                    [x2, fval, exitflag, output] = fminunc('OneExp', startxn, h,...
                        squeeze(data(subject).perf_w_on(lin,col,:,control)),echo');
                    res2(i,1:2) = abs(x2);
                    res2(i,3) = OneExp(abs(x2), squeeze(data(subject).perf_w_on(lin,col,:,control)),...
                        echo');
                end
                [zob, I1] = min(res1(:,3));
                [zob, I2] = min(res2(:,3));
                Maps(control).GradOFF(lin,col) = res1(I1,2);
                Maps(control).GradON(lin,col) = res2(I2,2);
            end
            end
        end
        
    else
        
        maximages = 200;
        
        figure(map1);
        
        subplot(3,3,control)
        image1 = Maps(control).GradOFF;
        imshow(image1,[0, maximages]);
        title(strcat(name))
        
        figure(map2);
        subplot(3,3, control)
        image2 = Maps(control).GradON;
        imshow(image2,[0, maximages]);
        title(strcat(name))
        
        
        
    end
end
save('voxMap.mat','Maps');
figure(map1);
suptitle('One Exponential fit of T2 Gradient off');

figure(map2);
suptitle('One Exponential fit of T2 Gradient on');

%% compute the mean of both Maps
control = 4;

for i=1:7
    MEANOFF(:,:,i) = Maps(i).GradOFF(:,:);
    MEANON(:,:,i) = Maps(i).GradON(:,:);
end
MEANOFF = MEANOFF(:,:,control);
MEANON = MEANON(:,:,control);
I = find(MEANOFF < prctile(MEANOFF(:),99));
I2 = find(MEANON < prctile(MEANON(:), 99));
MEANON = MEANON(I2);
MEANOFF = MEANOFF(I);
MEANON = MEANON(MEANON ~=0 );
MEANOFF = MEANOFF(MEANOFF ~=0 );
I = find(MEANOFF > prctile(MEANOFF(:),1));
I2 = find(MEANON > prctile(MEANON(:), 1));
display('grad on T2 mean')
mean(MEANON(I2))

display('grad off  T2 mean');
mean(MEANOFF(I))