clear all; close all;

fid=fopen('dwi.Bfloat', 'r', 'b');
dwis = fread(fid, 'float');
fclose(fid);


dwis = reshape(dwis, 33, 112, 112, 50);

qhat = load('grad_dirs.txt')';

bvals = 1000*sum(qhat.*qhat);

Avox = dwis(:,52,62,25);

i=0;
n = 0;
totn = 1;
TotIterations = 0;

tic;
while (n/totn < 0.4 || i < 50) && i < 1000
    i=i+1;
    
    if rem(i,10)==0
        fprintf('i = %i\n',i);
    end
    
    pd = makedist('Normal');
    t = truncate(pd,0,inf);
    r = random(t,2,1);
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
    
    h=optimset('MaxFunEvals',20000,...
        'Algorithm','levenberg-marquardt',...
        'LargeScale','off',...
        'MaxIter', 20000,...
        'Display','off',...
        'TolX',1e-10,...
        'TolFun',1e-10);
    % Now run the fitting
    [parameter_hat,RESNORM,EXITFLAG,OUTPUT]=fminunc('BallStickSSD',startx,h,Avox,bvals,qhat);
    parameter_hat(1) = abs(parameter_hat(1));
    parameter_hat(2) = abs(parameter_hat(2));
    parameter_hat(3) = sin(parameter_hat(3))^2;
    
    Resnorm(i)= RESNORM;
    n = sum(Resnorm==mode(Resnorm));
    totn = length(Resnorm);
    if Resnorm(i) == min(Resnorm)
        TrueParameters = parameter_hat;
    end
    
    TotIterations = TotIterations + OUTPUT.iterations;
    
end
AverageTime =  toc / i;
AverageIteration = TotIterations / i;

hist(Resnorm,100);
set(gca, 'FontSize', 14);
xlabel('SSD');
ylabel('# of solutions');