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
    image1 = mat2gray(data(1).perf_w_off(:,:,subject,im),[immin1, immax1]);
    imshow(image1,[]);
    title(strcat(pld,' / ',td));
    
    
    figure(fig2);
    subplot(2,4,im);
    image2 = mat2gray(data(1).perf_w_on(:,:,subject,im),[immin2, immax2]);
    imshow(image2,[]);
    title(strcat(pld,' / ',td));
end

ECHO = num2str(echo(1));
SUBJECT = num2str(subject);
figure(fig1);
ha = axes('Position',[0 0 1 1],'Xlim',[0 1],'Ylim',[0 1],'Box','off','Visible','off','Units','normalized', 'clipping' , 'off');
text(0.5, 1,strcat('Perfusion weighted images for patient: ',SUBJECT,' at echo: ', ECHO, ' with gradient off'),'HorizontalAlignment' ,'center','VerticalAlignment', 'top')
figure(fig2);
ha = axes('Position',[0 0 1 1],'Xlim',[0 1],'Ylim',[0 1],'Box','off','Visible','off','Units','normalized', 'clipping' , 'off');
text(0.5, 1,strcat('Perfusion weighted images for patient: ',SUBJECT,' at echo: ', ECHO, ' with gradient on'),'HorizontalAlignment' ,'center','VerticalAlignment', 'top')
hold off
% Part 5.2

% choose ROI for each subject
if false
    for i=1:9
        figure;
        mask(i) = struct('BW',[], 'xi', [],'yi', []);
        immax1 = max(max(data(i).perf_w_off(:,:,1,1)));
        immin1 = min(min(data(i).perf_w_off(:,:,1,1)));
        image1 = mat2gray(data(i).perf_w_off(:,:,1,1),[immin1, immax1]);
        [mask(i).BW, mask(i).xi, mask(i).yi] = roipoly(image1);
    end
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
    set(gca, 'LineWidth', 1, 'FontSize',16);
    xlabel('Echos')
    ylabel('Mean signal')
    title(name);
    hold off
    
    figure(figMs2);
    subplot(2,4,param);
    errorbar(echo, meanTOplot3, stdTOplot3)
    hold on
    errorbar(echo, meanTOplot4, stdTOplot4,'Color','r')
    set(gca, 'LineWidth', 1, 'FontSize',16);
    xlabel('Echos')
    ylabel('Mean signal')
    title(name);
    hold off    
    
end
            

figure(figMs1);
suptitle('The mean ASL signal of ROI with gradient off');
legend('ASL','Control','Location','EastOutside')
figure(figMs2);
suptitle('The mean ASL signal of ROI with gradient on');
legend('ASL','Control','Location','EastOutside')
hold off
% Part 5.3


if false
    % find the best T2 value with one exponential
    h=optimset('MaxFunEvals',200000,...
        'Algorithm','levenberg-marquardt',...
        'LargeScale','off',...
        'Display','off');
    
    Fit1 = struct('ASLon', [], 'ASLoff', [], 'CONTROLon' , [], 'CONTROLoff', []);
    Conf1 = struct('ASLon', [], 'ASLoff', [], 'CONTROLon' , [], 'CONTROLoff', []);
    
    for control=1:7
        for subject=1:9
            startx(1) = 1;
            fprintf('Control: %i and Subject: %i\n', control, subject);
            %mean(T2vsTE(control).ASLoff(:,subject));
            startx(2) = 120;
            [x1,fval,exitflag,output] = fminunc('OneExp',startx,h,T2vsTE(control).ASLoff(:,subject),echo');
            Fit1.ASLoff(subject,control) = abs(x1(2));
            [Conf1.ASLoff(subject,control,1), Conf1.ASLoff(subject,control,2), Conf1.ASLoff(subject,control,3), Conf1.ASLoff(subject,control,4)] = Bootstrap('OneExp', T2vsTE(control).ASLoff(:,subject),100,h , startx);
            
            [x2, fval, exitflag, output] = fminunc('OneExp', startx, h, T2vsTE(control).ASLon(:,subject),echo');
            Fit1.ASLon(subject,control) = abs(x1(2));
            [Conf1.ASLon(subject,control,1), Conf1.ASLon(subject,control,2), Conf1.ASLon(subject,control,3), Conf1.ASLon(subject,control,4)] = Bootstrap('OneExp', T2vsTE(control).ASLon(:,subject),100,h , startx);
            
            [x3, fval, exitflag, output] = fminunc('OneExp', startx, h, T2vsTE(control).CONTROLoff(:,subject),echo');
            Fit1.CONTROLoff(subject, control) = abs(x3(2));
            [Conf1.CONTROLoff(subject,control,1), Conf1.CONTROLoff(subject,control,2), Conf1.CONTROLoff(subject,control,3), Conf1.CONTROLoff(subject,control,4)] = Bootstrap('OneExp', T2vsTE(control).CONTROLoff(:,subject),100,h , startx);
            
            [x4, fval, exitflag, output] = fminunc('OneExp', startx, h, T2vsTE(control).CONTROLon(:,subject),echo');
            Fit1.CONTROLon(subject, control) = abs(x4(2));
            [Conf1.CONTROLon(subject,control,1), Conf1.CONTROLon(subject,control,2), Conf1.CONTROLon(subject,control,3), Conf1.CONTROLon(subject,control,4)] = Bootstrap('OneExp', T2vsTE(control).CONTROLon(:,subject),100,h , startx);
            
        end
    end
    
else
    load('Fit1');
    load('Conf1');
end

subject = 1;
h=optimset('MaxFunEvals',200000,...
    'Algorithm','levenberg-marquardt',...
    'LargeScale','off',...
    'Display','off');
startx(1) = 1;
startx(2) = 120;

map1 = figure;

for control=1:7
    if false
        name=strcat('tau: ', num2str(TD(control)),'   PLD: ', num2str(PLD(control)));
        fprintf('Control: %i/7\n', control)
        Maps(control) = struct('GradOFF', [], 'GradON', []);
        for lin=1:size(data(subject).perf_w_off,1)
            fprintf('Lin: %i/64\n', lin);
            for col=1:size(data(subject).perf_w_off,2)
                [x1, fval, exitflag, output] = fminunc('OneExp', startx, h, squeeze(data(subject).perf_w_off(lin,col,:,control)),echo');
                Maps(control).GradOFF(lin,col) = abs(x1(2));
                
                [x2, fval, exitflag, output] = fminunc('OneExp', startx, h, squeeze(data(subject).perf_w_off(lin,col,:,control)),echo');
                Maps(control).GradON(lin,col) = abs(x2(2));
            end
        end
        
        min1 = min(min(Maps(control).GradOFF));
        min2 = min(min(Maps(control).GradON));
        minall = min([min1, min2]);
        max1 = max(max(Maps(control).GradOFF));
        max2 = max(max(Maps(control).GradON));
        maxall = max([max1, max2]);
        
        figure(map1);
        
        subplot(2,7,control)
        image1 = mat2gray(Maps(control).GradOFF,[minall, maxall]);
        imshow(image1,[]);
        title(strcat(name,' ' ,'Grad Off'))
        
        subplot(2,7, control+7)
        image2 = mat2gray(Maps(control).GradON,[minall, maxall]);
        imshow(image2,[]);
        title(strcat(name,' ' ,'Grad On'))
    else
        load('Maps.mat');
        min1 = min(min(Maps(control).GradOFF));
        min2 = min(min(Maps(control).GradON));
        minall = min([min1, min2]);
        max1 = max(max(Maps(control).GradOFF));
        max2 = max(max(Maps(control).GradON));
        maxall = max([max1, max2]);
        
        figure(map1);
        
        subplot(2,7,control)
        image1 = mat2gray(Maps(control).GradOFF,[minall, maxall]);
        imshow(image1,[]);
        title(strcat(name,' ' ,'Grad Off'))
        
        subplot(2,7, control+7)
        image2 = mat2gray(Maps(control).GradON,[minall, maxall]);
        imshow(image2,[]);
        title(strcat(name,' ' ,'Grad On'))
    end
end


