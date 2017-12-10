% Q5.1

clear all; 
close all;
subject = 1;
PLD = [50, 50, 50, 50, 300, 700, 1200];
TD = [500, 750, 1000, 3000, 3000, 3000, 3000];
echo = [29, 41, 53, 68.2, 80.2, 92.2, 107.4, 119.4, 131.4, 146.6, 158.6, 170.6]; 

% Part 5.1

% imort the data and take the difference of control vs marquage
for i=1:9
    name = sprintf('base_images_subject_%i',i);
    field = sprintf('subject%i',i);
    value = importData(name);
    data(i) = struct('new_images_off', value.new_images_off, 'new_images_on', value.new_images_on, 'perf_w_off',[], 'perf_w_on',[]); 
    data(i).new_images_off = squeeze(data(i).new_images_off);
    data(i).new_images_on = squeeze(data(i).new_images_on);
    
    % cut the unuseful data outside of the brain
    data(i).new_images_off = data(i).new_images_off(10:35,15:50,:,:);
    data(i).new_images_on = data(i).new_images_on(10:35,15:50,:,:);

    %%%%% 
    % clalculate the perfusion weighted images
    
    for e=1:12
        for tp=2:2:14
            data(i).perf_w_off(:,:,e,tp/2) = (data(i).new_images_off(:,:,e,tp) - data(i).new_images_off(:,:,e,tp-1));
            data(i).perf_w_on(:,:,e,tp/2) = (data(i).new_images_on(:,:,e,tp) - data(i).new_images_on(:,:,e,tp-1));
        end
    end
    
end


% plot the figures for different xx combination
fig1 = figure;
fig2 = figure;
immax1 = max(max(data(1).perf_w_off(:,:,subject,1)));
immin1 = min(min(data(1).perf_w_off(:,:,subject,1)));
immax2 = max(max(data(1).perf_w_on(:,:,subject,1)));
immin2 = min(min(data(1).perf_w_on(:,:,subject,1)));

for im=1:7
    pld = num2str(PLD(im));
    td = num2str(TD(im));
    
    
    figure(fig1);
    subplot(2,4,im);
    %image1 = mat2gray(data(1).perf_w_off(:,:,1,im),[immin1, immax1]);
    %imshow(image1,[]);
    imshow(data(1).perf_w_off(:,:,1,im),[0,Inf]);
    title(strcat('pld: ',pld,'   tau: ',td));
    
    
    figure(fig2);
    subplot(2,4,im);
%     image2 = mat2gray(data(1).perf_w_on(:,:,1,im),[immin2, immax2]);
%     imshow(image2,[]);
    imshow(data(1).perf_w_on(:,:,1,im),[0,Inf]);
    title(strcat('pld: ',pld,'   tau: ',td));
end

ECHO = num2str(echo(1));
SUBJECT = num2str(subject);
figure(fig1);
suptitle(['Perfusion weighted images for patient:  ',SUBJECT,'   at echo:  ', ECHO, '   with gradient off']);
figure(fig2);
suptitle(['Perfusion weighted images for patient:  ',SUBJECT,'   at echo:  ', ECHO, '   with gradient on']);
hold off

save('data.mat','data');