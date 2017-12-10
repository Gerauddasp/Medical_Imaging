%% preparing the datas
clear all; 
close all;
subject = 1;
PLD = [50, 50, 50, 50, 300, 700, 1200];
TD = [500, 750, 1000, 3000, 3000, 3000, 3000];
echo = [29, 41, 53, 68.2, 80.2, 92.2, 107.4, 119.4, 131.4, 146.6, 158.6, 170.6];

for i=1:7
    Name{i} = strcat(num2str(TD(i)),(' / '),num2str(PLD(i)));
end

load('maskpoly.mat');
load('data');

ASL = struct('off',[],'on',[],'IV',[]);
Control = struct('off',[],'on',[],'IV',[]);
e = 1; % echo
for subject=1:9
        for control=1:7
            mat = mask(subject).BW .* data(subject).perf_w_on(:,:,e,control);
            ASL.on(subject, control) = mean(mat(mat~=0));
            
            mat2 = mask(subject).BW .* data(subject).perf_w_off(:,:,e,control);
            ASL.off(subject,control) = mean(mat2(mat2~=0));
            
            mat3 = mask(subject).BW .* data(subject).new_images_on(:,:,e,control * 2);
            Control.on(subject, control) = mean(mat3(mat3~=0));
            
            mat4 = mask(subject).BW .* data(subject).new_images_off(:,:,e, control * 2);
            Control.off(subject, control) = mean(mat4(mat4~=0));
        end
end

%% calulating phiIV

ASL.IV = 1 - (ASL.on ./ ASL.off);
Control.IV = 1 - (Control.on ./ Control.off);

save('Phi','ASL','Control');

%% Plotting results

figure;
X = (1:7);
bar(X, ASL.IV','b');
hold on
errorbar(X,mean(ASL.IV),std(ASL.IV),'r');
set(gca,'XTickLabel',Name)
xlabel('tau / PLD','fontsize',16);
ylabel('% IV','fontsize',16);
set(gca, 'LineWidth', 1, 'FontSize',16);

title('Percentage of Intra Vascular signal - ASL');

figure;
X = (1:7);
bar(X, Control.IV','b');
hold on
errorbar(X,mean(Control.IV),std(Control.IV),'r');

set(gca,'XTickLabel',Name)
xlabel('tau / PLD','fontsize',16);
ylabel('% IV','fontsize',16);
set(gca, 'LineWidth', 1, 'FontSize',16);

title('Percentage of Intra Vascular signal - Control');
