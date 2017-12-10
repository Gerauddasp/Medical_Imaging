clear all; close all;
ws = load('wsQ121');
Parameters_boot = ws.Parameters_boot;
x = ws.x;

figure;
S0 = Parameters_boot(:,1);
hist(S0,100);
hold on
sigmaS0 = std(S0);
meanS0 = mean(S0);
conf = prctile(S0,[2.5, 97.5]);
lb1 = meanS0 - 2 * sigmaS0;
ub1 = meanS0 + 2* sigmaS0;
set(gca, 'FontSize', 14);
xlabel('S0');
ylabel('# of solutions');

hSLines = plot([lb1, ub1], [200, 200],'r', [lb1, lb1], [170, 230], 'r',[ub1,ub1],[170,230],'r');
hCLines = plot([conf(1), conf(2)], [300,300], 'g', [conf(1), conf(1)], [270,330],'g', [conf(2), conf(2)],[270,330],'g');
hSGroup = hggroup;
hCGroup = hggroup;
set(hSLines,'Parent',hSGroup)
set(hCLines,'Parent',hCGroup)
% Include these hggroups in the legend:
set(get(get(hSGroup,'Annotation'),'LegendInformation'),...
    'IconDisplayStyle','on'); 
set(get(get(hCGroup,'Annotation'),'LegendInformation'),...
    'IconDisplayStyle','on'); 
legend('S0','2 Sigma','95%')
plot(x(1),0,'X','Linewidth',5,'Color','r');

figure;
d = Parameters_boot(:,2);
hist(d,100);
hold on
sigmad = std(d);
meand = mean(d);
conf = prctile(d,[2.5, 97.5]);
lb1 = meand - 2 * sigmad;
ub1 = meand + 2* sigmad;
set(gca, 'FontSize', 14);
xlabel('diffusion');
ylabel('# of solutions');

hSLines = plot([lb1, ub1], [200, 200],'r', [lb1, lb1], [170, 230], 'r',[ub1,ub1],[170,230],'r');
hCLines = plot([conf(1), conf(2)], [300,300], 'g', [conf(1), conf(1)], [270,330],'g', [conf(2), conf(2)],[270,330],'g');
hSGroup = hggroup;
hCGroup = hggroup;
set(hSLines,'Parent',hSGroup)
set(hCLines,'Parent',hCGroup)
% Include these hggroups in the legend:
set(get(get(hSGroup,'Annotation'),'LegendInformation'),...
    'IconDisplayStyle','on'); 
set(get(get(hCGroup,'Annotation'),'LegendInformation'),...
    'IconDisplayStyle','on'); 
legend('d','2 Sigma','95%')
plot(x(2),0,'X','Linewidth',5,'Color','r');

figure;
f = Parameters_boot(:,3);
hist(f,100);
hold on
sigmaf = std(f);
meanf = mean(f);
conf = prctile(f,[2.5, 97.5]);
lb1 = meanf - 2 * sigmaf;
ub1 = meanf + 2* sigmaf;
set(gca, 'FontSize', 14);
xlabel('f');
ylabel('# of solutions');

hSLines = plot([lb1, ub1], [200, 200],'r', [lb1, lb1], [170, 230], 'r',[ub1,ub1],[170,230],'r');
hCLines = plot([conf(1), conf(2)], [300,300], 'g', [conf(1), conf(1)], [270,330],'g', [conf(2), conf(2)],[270,330],'g');
hSGroup = hggroup;
hCGroup = hggroup;
set(hSLines,'Parent',hSGroup)
set(hCLines,'Parent',hCGroup)
% Include these hggroups in the legend:
set(get(get(hSGroup,'Annotation'),'LegendInformation'),...
    'IconDisplayStyle','on'); 
set(get(get(hCGroup,'Annotation'),'LegendInformation'),...
    'IconDisplayStyle','on'); 
legend('f','2 Sigma','95%')
plot(x(3),0,'X','Linewidth',5,'Color','r');