% Q5.2

%%%%
% Import data


clear all; 
close all;
subject = 1;
PLD = [50, 50, 50, 50, 300, 700, 1200];
TD = [500, 750, 1000, 3000, 3000, 3000, 3000];
echo = [29, 41, 53, 68.2, 80.2, 92.2, 107.4, 119.4, 131.4, 146.6, 158.6, 170.6]; 

% imort the data and take the difference of control vs marquage
load('data');

% choose ROI for each subject
if false
    for i=1:9
        figure;
        mask(i) = struct('BW',[], 'xi', [],'yi', []);
        immax1 = max(max(data(i).perf_w_off(:,:,1,7)));
        immin1 = min(min(data(i).perf_w_off(:,:,1,7)));
        image1 = mat2gray(data(i).perf_w_off(:,:,1,7),[0, immax1]);
        [mask(i).BW, mask(i).xi, mask(i).yi] = roipoly(image1);
    end
    save('maskpoly.mat','mask');
else
    load('maskpoly.mat');
end

figMs1 = figure;
figMs2 = figure;
s = 0;



for param=1:7
    % add the string for the titles here:
    name=strcat('tau: ', num2str(TD(param)),'   PLD: ', num2str(PLD(param)));
    T2vsTE(param) = struct('ASLon', [], 'ASLoff', [], 'CONTROLoff', [], 'CONTROLon', []);
    for e=1:length(echo)
        
        
        for subject=1:9
            normEcho1 = sum(sum((mask(subject).BW .* data(subject).perf_w_off(:,:,1,param)) ./ nnz(mask(subject).BW)));
            normEcho2 = sum(sum((mask(subject).BW .* data(subject).perf_w_on(:,:,1,param)) ./ nnz(mask(subject).BW)));
            normEcho3 = sum(sum((mask(subject).BW .* data(subject).new_images_off(:,:,1,param*2)) ./ nnz(mask(subject).BW)));
            normEcho4 = sum(sum((mask(subject).BW .* data(subject).new_images_on(:,:,1,param*2)) ./ nnz(mask(subject).BW)));
            
            T2vsTE(param).ASLoff(e,subject) = sum(sum((mask(subject).BW .* data(subject).perf_w_off(:,:,e,param)) ./ nnz(mask(subject).BW)))/normEcho1;
            T2vsTE(param).ASLon(e, subject) = sum(sum((mask(subject).BW .* data(subject).perf_w_on(:,:,e,param)) ./ nnz(mask(subject).BW)))/normEcho2;
            T2vsTE(param).CONTROLoff(e, subject) = sum(sum((mask(subject).BW .* data(subject).new_images_off(:,:,e,param*2)) ./ nnz(mask(subject).BW)))/normEcho3;
            T2vsTE(param).CONTROLon(e, subject) = sum(sum((mask(subject).BW .* data(subject).new_images_on(:,:,e,param*2)) ./ nnz(mask(subject).BW)))/normEcho4;            
            
        end
        
    end

    meanTOplot1 = mean(T2vsTE(param).ASLoff,2)';
    stdTOplot1 = std(T2vsTE(param).ASLoff,0,2)';
    
    meanTOplot2 = mean(T2vsTE(param).ASLon,2)';
    stdTOplot2 = std(T2vsTE(param).ASLon,0,2)';
    
    meanTOplot3 = mean(T2vsTE(param).CONTROLoff,2)';
    stdTOplot3 = std(T2vsTE(param).CONTROLoff,0,2)';
    
    meanTOplot4 = mean(T2vsTE(param).CONTROLon,2)';
    stdTOplot4 = std(T2vsTE(param).CONTROLon,0,2)';
    
    
    
    
    figure(figMs1);
    subplot(2,4,param);
    errorbar(echo, meanTOplot1, stdTOplot1)
    hold on
    errorbar(echo, meanTOplot2, stdTOplot2,'Color','r')
    %set(gca, 'LineWidth', 1, 'FontSize',16);
    xlabel('Echos')
    ylabel('Mean signal')
    title(name);
    hold off
    
    figure(figMs2);
    subplot(2,4,param);
    errorbar(echo, meanTOplot3, stdTOplot3)
    hold on
    errorbar(echo, meanTOplot4, stdTOplot4,'Color','r')
    %set(gca, 'LineWidth', 1, 'FontSize',16);
    xlabel('Echos')
    ylabel('Mean signal')
    title(name);
    hold off    
    
end

figure(figMs1);
suptitle('The mean signal of ROI with gradient off');
legend('ASL','Control','Location','EastOutside')
figure(figMs2);
suptitle('The mean signal of ROI with gradient on');
legend('ASL','Control','Location','EastOutside')
hold off

save('T2vsTE.mat','T2vsTE');
