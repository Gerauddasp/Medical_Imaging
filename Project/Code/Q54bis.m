% Q5.4
clear all; 
close all;
subject = 1;
PLD = [50, 50, 50, 50, 300, 700, 1200];
TD = [500, 750, 1000, 3000, 3000, 3000, 3000];
echo = [29, 41, 53, 68.2, 80.2, 92.2, 107.4, 119.4, 131.4, 146.6, 158.6, 170.6];

for i=1:7
    Name{i} = strcat(num2str(TD(i)),(' / '),num2str(PLD(i)));
end


NotFix = true;
Fix = true;
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

%% Here we are calculating T2IC and T2EC with 2 parameters for ASL and Control

if NotFix

% The parameter that are fix during our loops:
pd = makedist('Normal');
t = truncate(pd,0.8,1.2);

h=optimset('MaxFunEvals',200000,...
    'Algorithm','levenberg-marquardt',...
    'LargeScale','off',...
    'Display','off',...
    'TolX',1,...
    'TolFun',0.001);

for c=1:7
    fprintf('Control Settings: %i\n', c);
    for s=1:9
        % starting points for ASL
        startx(1) = mean(ASL(s).signal(:,c));
        startx(2) = 53;
        startx(3) = mean(ASL(s).signal(:,c));
        startx(4) = 153;
        
        % starting points for control
        startx2(1) = mean(Control(s).signal(:,c));
        startx2(2) = 53;
        startx2(3) = mean(Control(s).signal(:,c));
        startx2(4) = 153;
        
        for it=1:10
            % fit the ASL
            r = random(t,1,4);
            
            startxn = r .* startx;
            
            [x1,~,~,~] = fminunc('TwoExp',startxn,h, ASL(s).signal(:,c));
            x1 = abs(x1);
            if x1(2) > x1(4)
                temp = x1(1);
                x1(1) = x1(4);
                x1(4) = temp;
            end
            res1(it,1:4) = x1;
            res1(it,5) = TwoExp(x1, ASL(s).signal(:,c));
            
            % fit the control
            
            startx2n = r .* startx2;
            
            [x2, ~,~,~] = fminunc('TwoExp', startx2n, h, Control(s).signal(:,c));
            x2 = abs(x2);
            if x2(2) > x2(4)
                temp = x2(1);
                x2(1) = x2(4);
                x2(4) = temp;
            end
            res2(it,1:4) = x2;
            res2(it,5) = TwoExp(x2, Control(s).signal(:,c));
        end
        
        % save variables
        [~,I1] = min(res1(:,5));
        ASLNFIX(s,c,:) = res1(I1,:); % S0IC, T2IC, S0EC, T2EC, error
        
        [~,I2] = min(res2(:,5));
        ControlNFIX(s,c,:) = res2(I2,:); % S0IC, T2IC, S0EC, T2EC, error
        
        ErrorNCasl(s,c) = res1(I1,5);
        ErrorNCcontrol(s,c) = res2(I2,5);
        
    end
end

%% plotting S0 IC/EC

ASLICECNFIX = ASLNFIX(:,:,1) ./ (ASLNFIX(:,:,3) + ASLNFIX(:,:,1));
ControlICECNFIX = ControlNFIX(:,:,1) ./ (ControlNFIX(:,:,3) + ControlNFIX(:,:,1));

X = [1:7];

plot(X,mean(ASLICECNFIX),'-o','Color','b','LineWidth',2);
hold on
plot(X,mean(ControlICECNFIX),'-o','Color','r','LineWidth',2);
set(gca,'XTickLabel',Name)
xlabel('tau / PLD','fontsize',16);
ylabel('log(Mic / Mec)','fontsize',16);
set(gca, 'LineWidth', 1, 'FontSize',16);

legend('ASL','Control','Location','Best');
title('bi-exponential non constrain');




end

%% Here we are calculating T2IC and T2EC with 4 parameters for ASL and Control

clear pd t h c s startx startx2 x1 x2 res1 res2 ASLICEC ControlICEC

if Fix

% The parameter that are fix during our loops:
pd = makedist('Normal');
t = truncate(pd,0.8,1.2);

h=optimset('MaxFunEvals',200000,...
    'Algorithm','levenberg-marquardt',...
    'LargeScale','off',...
    'Display','off',...
    'TolX',1,...
    'TolFun',0.001);

for c=1:7
    fprintf('Control Settings: %i\n', c);
    for s=1:9
        % starting points for ASL
        startx(1) = mean(ASL(s).signal(:,c));
        startx(2) = 53;
        startx(3) = mean(ASL(s).signal(:,c));
        startx(4) = 153;
        
        % starting points for control
        startx2(1) = mean(Control(s).signal(:,c));
        startx2(2) = 53;
        startx2(3) = mean(Control(s).signal(:,c));
        startx2(4) = 153;
        
        for it=1:10
            % fit the ASL
            r = random(t,1,4);
            
            startxn = r .* startx;
            
            [x1,~,~,~] = fminunc('TwoExpFix',startxn,h, ASL(s).signal(:,c));
            x1 = abs(x1);
            if x1(2) > x1(4)
                temp = x1(1);
                x1(1) = x1(4);
                x1(4) = temp;
            end
            res1(it,1:4) = x1;
            res1(it,5) = TwoExp(x1, ASL(s).signal(:,c));
            
            % fit the control
            
            startx2n = r .* startx2;
            
            [x2, ~,~,~] = fminunc('TwoExpFix', startx2n, h, Control(s).signal(:,c));
            x2 = abs(x2);
            if x2(2) > x2(4)
                temp = x2(1);
                x2(1) = x2(4);
                x2(4) = temp;
            end
            res2(it,1:4) = x2;
            res2(it,5) = TwoExp(x2, Control(s).signal(:,c));
        end
        
        % save variables
        [~,I1] = min(res1(:,5));
        ASLFIX(s,c,:) = res1(I1,:); % S0IC, T2IC, S0EC, T2EC, error
        
        [~,I2] = min(res2(:,5));
        ControlFIX(s,c,:) = res2(I2,:); % S0IC, T2IC, S0EC, T2EC, error
        
        ErrorCasl(s,c) = res1(I1,5);
        ErrorCcontrol(s,c) = res2(I2,5);
        
    end
end

%% plotting S0 IC/EC
figure;
ASLICECFIX = ASLFIX(:,:,1) ./ (ASLFIX(:,:,3) + ASLFIX(:,:,1));
ControlICECFIX = ControlFIX(:,:,1) ./ (ControlFIX(:,:,3) + ControlFIX(:,:,1));

X = [1:7];

plot(X,mean(ASLICECFIX),'-o','Color','b','LineWidth',2);
hold on
plot(X,mean(ControlICECFIX),'-o','Color','r','LineWidth',2);
set(gca,'XTickLabel',Name)
xlabel('tau / PLD','fontsize',16);
ylabel('log(Mic / Mec)','fontsize',16);
set(gca, 'LineWidth', 1, 'FontSize',16);

legend('ASL','Control','Location','Best');
title('bi-exponential constrain');


save('ICEC','ASLICECFIX','ControlICECFIX', 'ASLICECNFIX', 'ControlICECFIX');
save('ErrorBiexp.mat', 'ErrorCasl','ErrorNCasl','ErrorCcontrol', 'ErrorNCcontrol');

end