% MCMC

close all; clear all;

fid = fopen('dwi.Bfloat', 'r', 'b');
dwis = fread(fid, 'float');
fclose(fid);


dwis = reshape(dwis, 33, 112, 112, 50);

qhat = load('grad_dirs.txt')';

bvals = 1000*sum(qhat.*qhat);

Avox = dwis(:,52,62,25);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% We take random initials points:

pd = makedist('Normal','mu',1,'sigma',2);
t1 = truncate(pd,-1,inf);
r = random(t1,2,1);
x(1) = 2.5E5 + (2.5E5) * r(1);
x(2) = 1E-3 + (1E-3) * r(2);

t2 = truncate(pd,-0.99,1);
r2 = random(t2,1,1);
x(3) = 0.5 + 0.5 * r2(1);

t3 = truncate(pd,-pi/2,pi/2);
r3 = random(t3,2,1);
x(4) = r3(1);
x(5) = r3(2);


i = 0;
% Loop from here
% Make transformation for parameters:
burn = 1000;

for t=1:10000;
    if rem(t,1000)==0
        fprintf('iteration number: %i\n',t);
    end
    % Now we calculate prob of A if x
    
    [err1, Sx] = BallStickSSDMCMC(x, Avox, bvals, qhat);
    x(1) = abs(x(1));
    x(2) = abs(x(2));
    x(3) = sin(x(3))^2;
    
    if t == burn
        Etrack(1) = err1;
    end
    
    sigma = mean(Avox-Sx);
    pAx = prod(exp(-((Avox-Sx).^2/(2*sigma)^2)));
    
    % We define a new point from the previous points:
    
    pd = makedist('Normal','mu',1,'sigma',1);
    t1 = truncate(pd,-0.01,0.01);
    r = random(t1,2,1);
    y(1) = x(1) + x(1) * r(1);
    
    
    t2 = truncate(pd,-0.001,0.001);
    r2 = random(t2,1,1);
    y(2) = x(2) + x(2) * r2;
    
    tf = truncate(pd,0,1);
    rf = random(tf,1,1);
    y(3) = rf ;

    
    
    t3 = truncate(pd,-pi/4,pi/4);
    r3 = random(t3,2,1);
    y(4) = x(4) + r3(1) * x(4);
    y(5) = x(4) + r3(2) * x(4);
    
    
    [err2, Sy] = BallStickSSDMCMC(y, Avox, bvals, qhat);
    y(1) = abs(y(1));
    y(2) = abs(y(2));
    y(3) = sin(y(3))^2;
    
    if t >= burn
        Etrack((t+1)-burn) = err2;
    end
    
    sigma = mean(Avox - Sy);
    pAy = prod(exp(-((Avox-Sy).^2/(2*sigma)^2)));
    
    ALPHA = pAy / pAx;
    lucky = rand;
    alpha(t) = ALPHA;
    
    if ALPHA > 1
        x = y;
        if t >= burn
            i=i+1;
            parameter_mcmc(i,:) = y;
        end
    elseif ALPHA > lucky
        x = y;
        if t >= burn
            i=i+1;
            parameter_mcmc(i,:) = y;
        end
    end
    
    
end

figure;
plot(log(Etrack));
title('log(SSD)');
xlabel('iteration');
ylabel('SSD');

figure;
S0 = parameter_mcmc(:,1);
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

figure;
d = parameter_mcmc(:,2);
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

figure;
f = parameter_mcmc(:,3);
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



