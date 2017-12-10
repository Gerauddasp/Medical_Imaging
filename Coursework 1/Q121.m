
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Question 1.2.1
clear all; close all;

fid=fopen('dwi.Bfloat', 'r', 'b');
dwis = fread(fid, 'float');
fclose(fid);

dwis = reshape(dwis, 33, 112, 112, 50);

qhat = load('grad_dirs.txt')';

bvals = 1000*sum(qhat.*qhat);

Avox = dwis(:,70,70,20);


i=0;

while i<100
    i=i+1;
    
    if rem(i,10)==0
        fprintf('i = %i\n',i);
    end    
    % Define a starting point for the non-linear fit
    pd = makedist('Normal','mu',1,'sigma',2);
    t1 = truncate(pd,-1,inf);
    r = random(t1,2,1);
    noise1 = 2.5E5 + (2.5E5) * r(1);
    noise2 = 1E-3 + (1E-3) * r(2);
    
    t2 = truncate(pd,-1,1);
    r2 = random(t2,1,1);
    noise3 = 0.5 + 0.5 * r2(1);
    
    t3 = truncate(pd,-pi/2,pi/2);
    r3 = random(t3,2,1);   
    noise4 = r3(1);
    noise5 = r3(2);    
    startx = [noise1, noise2, noise3, noise4, noise5];
    
    h=optimset('MaxFunEvals',1000,...
        'Algorithm','levenberg-marquardt',...
        'LargeScale','off',...
        'MaxIter', 1000,...
        'Display','off',...
        'TolX',1e-10,...
        'TolFun',1e-10);
    % Now run the fitting
    [parameter_hat,RESNORM,EXITFLAG,OUTPUT]=fminunc('BallStickSSD',startx,h,Avox,bvals,qhat);
    parameter_hat(1) = abs(parameter_hat(1));
    parameter_hat(2) = abs(parameter_hat(2));
    parameter_hat(3) = sin(parameter_hat(3))^2;
    
    Resnorm(i)= RESNORM;
    parameterszob(i,:) = parameter_hat;
    
    
    
end
[C, I] = min(Resnorm);

% Here we try to compare the results we obtain:
x = parameterszob(I,:);
S0 = x(1);
diff = x(2);
f = x(3);
theta = x(4);
phi = x(5);

% Synthesize the signals
fibdir = [cos(phi)*sin(theta) sin(phi)*sin(theta) cos(theta)];
fibdotgrad = sum(qhat.*repmat(fibdir, [length(qhat), 1])');
model = S0*(f*exp(-bvals*diff.*(fibdotgrad.^2)) + (1-f)*exp(-bvals*diff));


% calculate the A_hat:
K = 10000;
S = model';
E = std(S);
pd = makedist('Normal','mu',1,'sigma',1);
b = truncate(pd,-0.5,0.5);
n = random(b, size(S,1),K);
A_hat = repmat(S,1,K) + n .* E;

% Now estimates best fit for all A_hat:
fprintf('Now estimates best fit for all A_hat\n')

for k=1:K
    
    if rem(k,100) == 0
        fprintf('We are in bootstrap number: %f\n',k);
    end
    
    Avox = A_hat(:,k);
    
    i=0;
    
    while i<10
        i = i+1;
        
        pd = makedist('Normal','mu',1,'sigma',1);
        t1 = truncate(pd,-1,Inf);
        r = random(t1,2,1);
        noise1 = x(1) + x(1) * r(1);
        noise2 = x(2) + x(2) * r(2);
        
        t2 = truncate(pd,-1,Inf);
        r2 = random(t2,1,1);
        noise3 = x(3) + x(3) * r2(1) + 0.1;
        
        t3 = truncate(pd,-pi/2,pi/2);
        r3 = random(t3,2,1);
        noise4 = x(4) + r3(1);
        noise5 = x(5) + r3(2);


        startx = [noise1, noise2, noise3, noise4, noise5];
        
        h=optimset('MaxFunEvals',200000,...
            'Algorithm','levenberg-marquardt',...
            'LargeScale','off',...
            'MaxIter', 200000,...
            'Display','off',...
            'TolX',1e-10,...
            'TolFun',1e-10);
        % Now run the fitting
        [parameter_hat,RESNORM,EXITFLAG,OUTPUT]=fminunc('BallStickSSD',startx,h,Avox,bvals,qhat);
        parameter_hat(1) = abs(parameter_hat(1));
        parameter_hat(2) = abs(parameter_hat(2));
        parameter_hat(3) = sin(parameter_hat(3))^2;
        
        Resnorm2(i)= RESNORM;
        parameterszob2(i,:) = parameter_hat;
        

        
    end
    [C, I] = min(Resnorm2);
    Parameters_boot(k,:) = parameterszob2(I,:);
    R(k) = Resnorm2(I);
    
end

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