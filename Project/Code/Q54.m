% Q5.4
% clear all; 
close all;
subject = 1;
PLD = [50, 50, 50, 50, 300, 700, 1200];
TD = [500, 750, 1000, 3000, 3000, 3000, 3000];
echo = [29, 41, 53, 68.2, 80.2, 92.2, 107.4, 119.4, 131.4, 146.6, 158.6, 170.6];

for i=1:7
    Name{i} = strcat(num2str(TD(i)),(' / '),num2str(PLD(i)));
end

%load('T2vsTE')

% Retrain parameters
NotFix = false;
Fix = false;
%% preparing the datas:
load('maskpoly.mat');
load('data');

% taking the ASL signal of the ROI defined in Q52 with diffusion gradients
% for vascular suppression on
ASL = struct('signal',[]);
for subject=1:9
    for e=1:12
        for control=1:7
            mat = mask(subject).BW .* data(subject).perf_w_on(:,:,e,control);
            ASL(subject).signal(e,control) = mean(mat(mat~=0));
        end
    end
end


% taking the Control signal of the ROI defined in Q52 with diffusion gradients
% for vascular suppression on
Control = struct('signal',[]);
for subject=1:9
    for e=1:12
        for control=2:2:14
            mat = mask(subject).BW .* data(subject).new_images_on(:,:,e,control);
            Control(subject).signal(e,control/2) = mean(mat(mat~=0));
        end
    end
end


%% Here we are calculating T2IC and T2EC for every parameters combinations time with ASL and Control

if NotFix

% The parameter that are fix during our loops:
pd = makedist('Normal');
t = truncate(pd,0.8,1.2);

h=optimset('MaxFunEvals',200000,...
    'Algorithm','levenberg-marquardt',...
    'LargeScale','off',...
    'Display','off',...
    'TolX',1e-10,...
    'TolFun',1e-1);

ASLS0 = struct('IC',[], 'EC',[]);
ASLT2 = struct('IC',[],'EC',[],'error',[]);
ASLBootStrap1 = struct('s2lbS0IC', [], 's2ubS0IC', [], 's2lbT2IC', [], 's2ubT2IC',...
    [], 's2lbS0EC', [], 's2ubS0EC', [], 's2lbT2EC',[] , 's2ubT2EC', []);

ControlS0 = struct('IC',[], 'EC',[]);
ControlT2 = struct('IC',[],'EC',[],'error',[]);
ControlBootStrap1 = struct('s2lbS0IC', [], 's2ubS0IC', [], 's2lbT2IC', [], 's2ubT2IC',...
    [], 's2lbS0EC', [], 's2ubS0EC', [], 's2lbT2EC',[] , 's2ubT2EC', []);
for control=1:7

    fprintf('Control Settings: %i\n', control);
    
    for subject=1:9
        % ASL starting points
        startx(1) = mean(mean(ASL(subject).signal(:,control))); % S0 IC
        startx(2) = 53; % T2 IC
        startx(3) = mean(mean(ASL(subject).signal(:,control))); % S0 EC
        startx(4) = 133; % T2 EC

        % control starting points
        startx2(1) = mean(mean(Control(subject).signal(:,control))); % S0 IC
        startx2(2) = 53; % T2 IC
        startx2(3) = mean(mean(Control(subject).signal(:,control))); % S0 EC
        startx2(4) = 133; % T2 EC
        for i=1:10
            % for the ASL
            r1 = random(t,1,1);
            r2 = random(t,1,1);
            r3 = random(t,1,1);
            r4 = random(t,1,1);
            startxn(1) = startx(1)*r1;
            startxn(2) = startx(2)*r2;
            startxn(3) = startx(3)*r3;
            startxn(4) = startx(4)*r4;
            [x1, fval, exitflag, output] = fminunc('TwoExp',startxn,h, ASL(subject).signal(:,control));
            if x1(2) > x1(4)
                temp = x1(1);
                x1(1) = x1(4);
                x1(4) = temp;
            end
            res(i,1) = abs(x1(1));
            res(i,2) = abs(x1(2));
            res(i,3) = abs(x1(3));
            res(i,4) = abs(x1(4));
            res(i,5) = TwoExp(x1,ASL(subject).signal(:,control));
            % for control
            r1 = random(t,1,1);
            r2 = random(t,1,1);
            r3 = random(t,1,1);
            r4 = random(t,1,1);
            startx2n(1) = startx(1)*r1;
            startx2n(2) = startx(2)*r2;
            startx2n(3) = startx(3)*r3;
            startx2n(4) = startx(4)*r4;
            [x2, fval, exitflag, output] = fminunc('TwoExp',startx2n,h, Control(subject).signal(:,control));
            if x2(2) > x2(4)
                temp = x2(1);
                x2(1) = x2(4);
                x2(4) = temp;
            end
            res2(i,1) = abs(x2(1));
            res2(i,2) = abs(x2(2));
            res2(i,3) = abs(x2(3));
            res2(i,4) = abs(x2(4));
            res2(i,5) = TwoExp(x2,Control(subject).signal(:,control));
        end
        
        % ASL results
        [zob, I] = min(res(:,5));
        ASLS0.IC(control,subject) = res(I,1);
        ASLT2.IC(control, subject) = res(I,2);
        ASLS0.EC(control,subject) = res(I,3);
        ASLT2.EC(control, subject) = res(I,4);
        ASLT2.error(control, subject) = res(I,5);
        % save data for the plot Mic/Mec
        MicMec(control,subject) = res(I,1)/ (res(I,3) + res(I,1));
        
        % Do Bootstrap Here:
        [ASLBootStrap1.s2lbS0IC(control,subject), ASLBootStrap1.s2ubS0IC(control,subject),...
            ASLBootStrap1.s2lbT2IC(control,subject), ASLBootStrap1.s2ubT2IC(control,subject),...
            ASLBootStrap1.s2lbS0EC(control,subject), ASLBootStrap1.s2ubS0EC(control,subject),...
            ASLBootStrap1.s2lbT2EC(control,subject), ASLBootStrap1.s2ubT2EC(control,subject)] = ...
            BootstrapTwo('TwoExp',ASL(subject).signal(:,control),100,h,startx); 
       
        % Control results
        [zob, I2] = min(res(:,5));
        ControlS0.IC(control,subject) = res2(I2,1);
        ControlT2.IC(control, subject) = res2(I2,2);
        ControlS0.EC(control,subject) = res2(I2,3);
        ControlT2.EC(control, subject) = res2(I2,4);
        ControlT2.error(control, subject) = res2(I2,5);
        % save data for the plot Mic/Mec
        MicMec2(control,subject) = res2(I,1)/ (res2(I,3) + res2(I,1));
        
        % Do Bootstrap Here:
        [ControlBootStrap1.s2lbS0IC(control,subject), ControlBootStrap1.s2ubS0IC(control,subject),...
            ControlBootStrap1.s2lbT2IC(control,subject), ControlBootStrap1.s2ubT2IC(control,subject),...
            ControlBootStrap1.s2lbS0EC(control,subject), ControlBootStrap1.s2ubS0EC(control,subject),...
            ControlBootStrap1.s2lbT2EC(control,subject), ControlBootStrap1.s2ubT2EC(control,subject)] = ...
            BootstrapTwo('TwoExp',Control(subject).signal(:,control),100,h,startx2);
    end
end
save('TwoExpNotFix.mat','ASLS0','ASLT2','ASLBootStrap1','ControlS0','ControlT2','ControlBootStrap1','MicMec','MicMec2');
else
load('TwoExpNotFix');
end
%% Plotting the Confidence intervals for one subject - ASL signal
subject = 1;

Intervals = figure;
S0ICcurve = ASLS0.IC(:,subject);
S0IClbcurve = ASLBootStrap1.s2lbS0IC(:,subject);
S0ICubcurve = ASLBootStrap1.s2ubS0IC(:,subject);

T2ICcuvre = ASLT2.IC(:,subject);
T2IClbcurve = ASLBootStrap1.s2lbT2IC(:,subject);
T2ICubcurve = ASLBootStrap1.s2ubT2IC(:,subject);

S0ECcuvre = ASLS0.EC(:,subject);
S0EClbcurve = ASLBootStrap1.s2lbS0EC(:,subject);
S0ECubcurve = ASLBootStrap1.s2ubS0EC(:,subject);

T2ECcurve = ASLT2.EC(:,subject);
T2EClbcurve = ASLBootStrap1.s2lbT2EC(:,subject);
T2ECubcurve = ASLBootStrap1.s2ubT2EC(:,subject);

X = (1:7);
subplot(4,1,1)
plot(X,S0ICcurve,'b',X,S0IClbcurve,'r',...
    X,S0ICubcurve,'g');
title('Intra-Cellular Signal');
set(gca, 'LineWidth', 1, 'FontSize',16);
xlabel('                                                                                                                                                                  PLD / tau')
ylabel('Signals')


subplot(4,1,2)
plot(X,T2ICcuvre,'b',X,T2IClbcurve,'r',...
    X,T2ICubcurve,'g');
title('Intra-Cellular T2 value');
set(gca, 'LineWidth', 1, 'FontSize',16);
xlabel('                                                                                                                                                                  PLD / tau')
ylabel('T2 (in ms)')

subplot(4,1,3)
plot(X,S0ECcuvre,'b',X,S0EClbcurve,'r',...
    X,S0ECubcurve,'g');
title('Extra-Cellular Signal');
set(gca, 'LineWidth', 1, 'FontSize',16);
xlabel('                                                                                                                                                                  PLD / tau')
ylabel('Signal')

subplot(4,1,4)
plot(X,T2ECcurve,'b',X,T2EClbcurve,'r',...
    X,T2ECubcurve,'g');
title('Extra-Cellular T2 value');
set(gca, 'LineWidth', 1, 'FontSize',16);
xlabel('                                                                                                                                                                  PLD / tau')
ylabel('T2 (in ms)')
legend('parameter','lower bound','upper bound')
suptitle('2 sigma confidence intervals for bi-exponential model with 4 parameters - ASL signal')


%% Plotting the Confidence intervals for one subject - Control Signal
Intervals2 = figure;

S0ICcurve2 = ControlS0.IC(:,subject);
S0IClbcurve2 = ControlBootStrap1.s2lbS0IC(:,subject);
S0ICubcurve2 = ControlBootStrap1.s2ubS0IC(:,subject);

T2ICcuvre2 = ControlT2.IC(:,subject);
T2IClbcurve2 = ControlBootStrap1.s2lbT2IC(:,subject);
T2ICubcurve2 = ControlBootStrap1.s2ubT2IC(:,subject);

S0ECcuvre2 = ControlS0.EC(:,subject);
S0EClbcurve2 = ControlBootStrap1.s2lbS0EC(:,subject);
S0ECubcurve2 = ControlBootStrap1.s2ubS0EC(:,subject);

T2ECcurve2 = ControlT2.EC(:,subject);
T2EClbcurve2 = ControlBootStrap1.s2lbT2EC(:,subject);
T2ECubcurve2 = ControlBootStrap1.s2ubT2EC(:,subject);

X = (1:7);
subplot(4,1,1)
plot(X,S0ICcurve2,'b',X,S0IClbcurve2,'r',...
    X,S0ICubcurve2,'g');
title('Intra-Cellular Signal');
set(gca, 'LineWidth', 1, 'FontSize',16);
xlabel('                                                                                                                                                                  PLD / tau')
ylabel('Signals')


subplot(4,1,2)
plot(X,T2ICcuvre2,'b',X,T2IClbcurve2,'r',...
    X,T2ICubcurve2,'g');
title('Intra-Cellular T2 value');
set(gca, 'LineWidth', 1, 'FontSize',16);
xlabel('                                                                                                                                                                  PLD / tau')
ylabel('T2 (in ms)')

subplot(4,1,3)
plot(X,S0ECcuvre2,'b',X,S0EClbcurve2,'r',...
    X,S0ECubcurve2,'g');
title('Extra-Cellular Signal');
set(gca, 'LineWidth', 1, 'FontSize',16);
xlabel('                                                                                                                                                                  PLD / tau')
ylabel('Signals')

subplot(4,1,4)
plot(X,T2ECcurve2,'b',X,T2EClbcurve2,'r',...
    X,T2ECubcurve2,'g');
title('Extra-Cellular T2 value');
set(gca, 'LineWidth', 1, 'FontSize',16);
xlabel('                                                                                                                                                                  PLD / tau')
ylabel('T2 (in ms)')
legend('parameter','lower bound','upper bound')
suptitle('2 sigma confidence intervals for bi-exponential model with 4 parameters - Control signal')

%Here we plot Mic/Mec against all control parameters:
%% for ASL
MicMecFig = figure;
figure(MicMecFig);
% Colors = {'b','r','y','g','c','m','k',[165/255, 42/255, 42/255],[255/255, 127/255, 80/255]};
% for i=1:9
%     semilogy([1:7],MicMec(:,i),'LineWidth',2,'Color',Colors{i});
%     hold on
% end
% legend('Rat1','Rat2','Rat3','Rat4','Rat5','Rat6','Rat7','Rat8','Rat9','Location','best');
plot([1:7],mean(MicMec,2), '-o','LineWidth',2);
%plot([1:7],mean(MicMec,2));
hold on
title('The ratio of the MRI signal deriving from the intra to the extra-cellular space (4 parameters) - ASL');
set(gca,'XTickLabel',Name)
xlabel('tau / PLD','fontsize',16);
ylabel('log(Mic / Mec)','fontsize',16);
set(gca, 'LineWidth', 1, 'FontSize',16);
% add line in y = 0
%graph2d.constantline(0,'Color',[0.7, 0.7, 0.7]);
%plot([1:7],(zeros(1,7)),'Color',[0.7,0.7,0.7],'LineWidth',4);

%% for control
% MicMecFig2 = figure;
% figure(MicMecFig2);
% Colors = {'b','r','y','g','c','m','k',[165/255, 42/255, 42/255],[255/255, 127/255, 80/255]};
% for i=1:9
%     semilogy([1:7],MicMec2(:,i),'LineWidth',2,'Color',Colors{i});
%     hold on
% end
% legend('Rat1','Rat2','Rat3','Rat4','Rat5','Rat6','Rat7','Rat8','Rat9','Location','best');
plot([1:7],mean(MicMec2,2),'-o','LineWidth',2,'Color','r');
%plot([1:7],mean(MicMec2,2),'Color','r');
legend('ASL','Control');
%hold on
% title('The ratio of the MRI signal deriving from the intra to the extra-cellular space (4 parameters) - Control');
% set(gca,'XTickLabel',Name)
% xlabel('tau / PLD','fontsize',16);
% ylabel('log(Mic / Mec)','fontsize',16);
% set(gca, 'LineWidth', 1, 'FontSize',16);
% add line in y = 0
%graph2d.constantline(0,'Color',[0.7, 0.7, 0.7]);
%plot([1:7],(zeros(1,7)),'Color',[0.7,0.7,0.7],'LineWidth',4);
%% Calculating the AIC
AICtwoExp = 8 + 12 * log( (1/12) * ASLT2.error );


%% Now we do the same but fixing the values of T2
if Fix
    clear startx startx2

% The parameter that are fix during our loops:
pd = makedist('Normal');
t = truncate(pd,0.8,1.2);

h=optimset('MaxFunEvals',200000,...
    'Algorithm','levenberg-marquardt',...
    'LargeScale','off',...
    'Display','off',...
    'TolX',1e-10,...
    'TolFun',1e-1);

ASLS0FIX = struct('IC',[], 'EC',[],'error',[]);
ASLBootStrap1FIX = struct('s2lbS0IC', [], 's2ubS0IC', [], 's2lbT2IC', [], 's2ubT2IC',...
    [], 's2lbS0EC', [], 's2ubS0EC', [], 's2lbT2EC',[] , 's2ubT2EC', []);

ControlS0FIX = struct('IC',[], 'EC',[],'error',[]);
ControlBootStrap1FIX = struct('s2lbS0IC', [], 's2ubS0IC', [], 's2lbT2IC', [], 's2ubT2IC',...
    [], 's2lbS0EC', [], 's2ubS0EC', [], 's2lbT2EC',[] , 's2ubT2EC', []);
for control=1:7

    fprintf('Control Settings: %i\n', control);
    
    for subject=1:9
        % ASL starting points
        startx(1) = mean(mean(ASL(subject).signal(:,control))); % S0 IC
        startx(2) = mean(mean(ASL(subject).signal(:,control))); % S0 EC

        % control starting points
        startx2(1) = mean(mean(Control(subject).signal(:,control))); % S0 IC
        startx2(2) = mean(mean(Control(subject).signal(:,control))); % S0 EC
        for i=1:10
            % for the ASL
            r1 = random(t,1,1);
            r2 = random(t,1,1);
            startx(1) = startx(1)*r1;
            startx(2) = startx(2)*r2;
            [x1, fval, exitflag, output] = fminunc('TwoExpFix',startx,h, ASL(subject).signal(:,control));
            res(i,1) = abs(x1(1));
            res(i,2) = abs(x1(2));
            res(i,3) = TwoExpFix(x1,ASL(subject).signal(:,control));
            % for control
            r1 = random(t,1,1);
            r2 = random(t,1,1);
            startx2(1) = startx(1)*r1;
            startx2(2) = startx(2)*r2;
            [x2, fval, exitflag, output] = fminunc('TwoExpFix',startx2,h, Control(subject).signal(:,control));
            res2(i,1) = abs(x2(1));
            res2(i,2) = abs(x2(2));
            res2(i,3) = TwoExpFix(x2,ASL(subject).signal(:,control));
        end
        
        % ASL results
        [zob, I] = min(res(:,3));
        ASLS0FIX.IC(control,subject) = res(I,1);
        ASLS0FIX.EC(control,subject) = res(I,2);
        ASLS0FIX.error(control, subject) =res(I,3);
        % save data for the plot Mic/Mec
        MicMecFIX(control,subject) = res(I,1)/ (res(I,2)+ res(I,1));
        
        % Do Bootstrap Here:
        [ASLBootStrap1FIX.s2lbS0IC(control,subject), ASLBootStrap1FIX.s2ubS0IC(control,subject),...
            ASLBootStrap1FIX.s2lbS0EC(control,subject), ASLBootStrap1FIX.s2ubS0EC(control,subject)] = ...
            BootstrapTwo('TwoExpFix',ASL(subject).signal(:,control),100,h,startx); 
       
        % Control results
        [zob, I2] = min(res(:,3));
        ControlS0FIX.IC(control,subject) = res2(I2,1);
        ControlS0FIX.EC(control,subject) = res2(I2,2);
        ControlS0FIX.error(control, subject) = res2(I2,3);
        % save data for the plot Mic/Mec
        MicMec2FIX(control,subject) = res2(I,1)/(res2(I,2) + res2(I,1));
        
        % Do Bootstrap Here:
        [ControlBootStrap1FIX.s2lbS0IC(control,subject), ControlBootStrap1FIX.s2ubS0IC(control,subject),...
            ControlBootStrap1FIX.s2lbS0EC(control,subject), ControlBootStrap1FIX.s2ubS0EC(control,subject)] =...
            BootstrapTwo('TwoExpFix',Control(subject).signal(:,control),100,h,startx2);
    end
end
save('TwoExpFix.mat','ASLS0FIX','ASLBootStrap1FIX','ControlS0FIX','ControlBootStrap1FIX','MicMecFIX','MicMec2FIX');
else
load('TwoExpFix');
end

%% Plotting the Confidence intervals for one subject - ASL signal
subject = 1;
IntervalsFIX = figure;
S0ICcurveFIX = ASLS0FIX.IC(:,subject);
S0IClbcurveFIX = ASLBootStrap1FIX.s2lbS0IC(:,subject);
S0ICubcurveFIX = ASLBootStrap1FIX.s2ubS0IC(:,subject);


S0ECcuvreFIX = ASLS0FIX.EC(:,subject);
S0EClbcurveFIX = ASLBootStrap1FIX.s2lbS0EC(:,subject);
S0ECubcurveFIX = ASLBootStrap1FIX.s2ubS0EC(:,subject);


X = (1:7);
subplot(2,1,1)
plot(X,S0ICcurveFIX,'b',X,S0IClbcurveFIX,'r',...
    X,S0ICubcurveFIX,'g');
title('Intra-Cellular Signal');
set(gca, 'LineWidth', 1, 'FontSize',16);
xlabel('                                                                                                                                                                  PLD / tau')
ylabel('Signals')

subplot(2,1,2)
plot(X,S0ECcuvreFIX,'b',X,S0EClbcurveFIX,'r',...
    X,S0ECubcurveFIX,'g');
title('Extra-Cellular Signal');
set(gca, 'LineWidth', 1, 'FontSize',16);
xlabel('                                                                                                                                                                  PLD / tau')
ylabel('Signals')

legend('parameter','lower bound','upper bound')
suptitle('2 sigma confidence intervals for bi-exponential model with 2 parameters - ASL signal')


%% Plotting the Confidence intervals for one subject - Control Signal
Intervals2FIX = figure;
S0ICcurve2FIX = ControlS0FIX.IC(:,subject);
S0IClbcurve2FIX = ControlBootStrap1FIX.s2lbS0IC(:,subject);
S0ICubcurve2FIX = ControlBootStrap1FIX.s2ubS0IC(:,subject);


S0ECcuvre2FIX = ControlS0FIX.EC(:,subject);
S0EClbcurve2FIX = ControlBootStrap1FIX.s2lbS0EC(:,subject);
S0ECubcurve2FIX = ControlBootStrap1FIX.s2ubS0EC(:,subject);

X = (1:7);
subplot(2,1,1)
plot(X,S0ICcurve2FIX,'b',X,S0IClbcurve2FIX,'r',...
    X,S0ICubcurve2FIX,'g');
title('Intra-Cellular Signal');
set(gca, 'LineWidth', 1, 'FontSize',16);
xlabel('                                                                                                                                                                  PLD / tau')
ylabel('Signals')


subplot(2,1,2)
plot(X,S0ECcuvre2FIX,'b',X,S0EClbcurve2FIX,'r',...
    X,S0ECubcurve2FIX,'g');
title('Extra-Cellular Signal');
set(gca, 'LineWidth', 1, 'FontSize',16);
xlabel('                                                                                                                                                                  PLD / tau')
ylabel('Signals')

legend('parameter','lower bound','upper bound')
suptitle('2 sigma confidence intervals for bi-exponential model with 2 parameters - Control signal')
%% plotting Mic/Mec
MicMecFixFig = figure;
figure(MicMecFixFig);
%Colors = {'b','r','y','g','c','m','k',[165/255, 42/255, 42/255],[255/255, 127/255, 80/255]};
%for i=1:9
    %semilogy([1:7],MicMecFIX(:,i),'LineWidth',2,'Color',Colors{i});
    %hold on
%end
plot([1:7],mean(MicMecFIX,2),'-o','LineWidth',2)
hold on
%legend('Rat1','Rat2','Rat3','Rat4','Rat5','Rat6','Rat7','Rat8','Rat9','Location','best');
title('The ratio of the MRI signal deriving from the intra to the extra-cellumar space (2 parameters) - ASL');
set(gca,'XTickLabel',Name)
xlabel('tau / PLD','fontsize',16);
ylabel('log(Mic / Mec)','fontsize',16);
set(gca, 'LineWidth', 1, 'FontSize',16);
% add line in y=0
%plot([1:7],(zeros(1,7)),'Color',[0.7,0.7,0.7],'LineWidth',4);
%graph2d.constantline(0,'Color',[0.7, 0.7, 0.7]);

%% for control
% MicMecFig2FIX = figure;
% figure(MicMecFig2FIX);
% Colors = {'b','r','y','g','c','m','k',[165/255, 42/255, 42/255],[255/255, 127/255, 80/255]};
% for i=1:9
%     semilogy([1:7],MicMec2FIX(:,i),'LineWidth',2,'Color',Colors{i});
%     hold on
% end
% legend('Rat1','Rat2','Rat3','Rat4','Rat5','Rat6','Rat7','Rat8','Rat9','Location','best');
plot([1:7],mean(MicMec2FIX,2),'-o','LineWidth',2,'Color','r');
legend('ASL','Control')
% hold on
% title('The ratio of the MRI signal deriving from the intra to the extra-cellular space (2 parameters) - Control');
% set(gca,'XTickLabel',Name)
% xlabel('tau / PLD','fontsize',16);
% ylabel('log(Mic / Mec)','fontsize',16);
% set(gca, 'LineWidth', 1, 'FontSize',16);
% add line in y = 0
%graph2d.constantline(0,'Color',[0.7, 0.7, 0.7]);
%plot([1:7],(zeros(1,7)),'Color',[0.7,0.7,0.7],'LineWidth',4);
%% And now we calculate the AIC

AICtwoEXPfix = 4 + 12 * log( (1/12) * ASLS0FIX.error);

%% computing percentage of the signal coming from the IC and EC
