clear all;
close all;

load('Phi');
load('ICEC');

PLD = [50, 50, 50, 50, 300, 700, 1200];
TD = [500, 750, 1000, 3000, 3000, 3000, 3000];
echo = [29, 41, 53, 68.2, 80.2, 92.2, 107.4, 119.4, 131.4, 146.6, 158.6, 170.6];

for i=1:7
    Name{i} = strcat(num2str(TD(i)),(' / '),num2str(PLD(i)));
end

%% ASL
IVasl = ASL.IV;
ICasl = (ones(size(ASL.IV)) - ASL.IV) .* ASLICECFIX;
ECasl = (ones(size(ASL.IV)) - ASL.IV) .* (ones(size(ASLICECFIX)) - ASLICECFIX);

PLOTasl = [mean(ECasl) ; mean(ICasl); mean(IVasl)];

bar(PLOTasl','stacked');
hold on
errorbar(mean(ECasl),std(ECasl),'Color','y');
errorbar(mean(ECasl) + mean(ICasl),std(ICasl),'Color','m');
errorbar(mean(ECasl) + mean(ICasl) + mean(IVasl), std(IVasl),'Color','c');


set(gca,'XTickLabel',Name)
xlabel('tau / PLD','fontsize',16);
ylabel('Contribution from each compartment','fontsize',16);
set(gca, 'LineWidth', 1, 'FontSize',16);

legend('EC','IC','IV' ,'Location','Best');
title('Compartmental distribution ASL');


%% Control
figure;
IVasl = Control.IV;
ICasl = (ones(size(Control.IV)) - Control.IV) .* ControlICECFIX;
ECasl = (ones(size(Control.IV)) - Control.IV) .* (ones(size(ControlICECFIX)) - ControlICECFIX);

PLOTasl = [mean(ECasl) ; mean(ICasl); mean(IVasl)];

bar(PLOTasl','stacked');
hold on
errorbar(mean(ECasl),std(ECasl),'Color','y');
errorbar(mean(ECasl) + mean(ICasl),std(ICasl),'Color','m');
errorbar(mean(ECasl) + mean(ICasl) + mean(IVasl), std(IVasl),'Color','c');


set(gca,'XTickLabel',Name)
xlabel('tau / PLD','fontsize',16);
ylabel('Contribution from each compartment','fontsize',16);
set(gca, 'LineWidth', 1, 'FontSize',16);

legend('EC','IC','IV' ,'Location','Best');
title('Compartmental distribution Control');
