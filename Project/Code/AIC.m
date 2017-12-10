clear all;
close all;
% COMPUTING AIC VALUES

load('ErrorsOneExpGradon.mat');
load('ErrorBiexp.mat');


PLD = [50, 50, 50, 50, 300, 700, 1200];
TD = [500, 750, 1000, 3000, 3000, 3000, 3000];
echo = [29, 41, 53, 68.2, 80.2, 92.2, 107.4, 119.4, 131.4, 146.6, 158.6, 170.6];

for i=1:7
    Name{i} = strcat(num2str(TD(i)),(' / '),num2str(PLD(i)));
end

%AICtwoEXPfix = 4 + 12 * log( (1/12) * ASLS0FIX.error);
%% for ASL

AICmono = aicbic(sum(ErrorASLon), repmat(2,1,7),12);
AICbiNC = aicbic(sum(ErrorNCasl), repmat(4,1,7),12);
AICbiC = aicbic(sum(ErrorCasl), repmat(2,1,7),12);

bar([AICmono', AICbiNC', AICbiC']);

set(gca,'XTickLabel',Name)
xlabel('tau / PLD','fontsize',16);
ylabel('AIC','fontsize',16);
set(gca, 'LineWidth', 1, 'FontSize',16);

legend('mono','bi-exp NC','bi-exp C' ,'Location','Best');
title('AIC for ASL');

%% for control
figure;
AICmono = aicbic(sum(ErrorCONTROLon), repmat(2,1,7),12);
AICbiNC = aicbic(sum(ErrorNCcontrol), repmat(4,1,7),12);
AICbiC = aicbic(sum(ErrorCcontrol), repmat(2,1,7),12);

bar([AICmono', AICbiNC', AICbiC']);

set(gca,'XTickLabel',Name)
xlabel('tau / PLD','fontsize',16);
ylabel('AIC','fontsize',16);
set(gca, 'LineWidth', 1, 'FontSize',16);

legend('mono','bi-exp NC','bi-exp C' ,'Location','Best');
title('AIC for control');
