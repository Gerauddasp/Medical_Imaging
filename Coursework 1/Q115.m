% Q115

clear all; close all;

fid=fopen('dwi.Bfloat', 'r', 'b');
dwis = fread(fid, 'float');
fclose(fid);


dwis = reshape(dwis, 33, 112, 112, 50);

qhat = load('grad_dirs.txt')';

bvals = 1000*sum(qhat.*qhat);

TotIterations = 0;

Avox = dwis(:,52,62,25);

q = qhat';
b = bvals';

y = [ones(size(b)), - b .* q(:,1).^2, -2 * b .* q(:,1) .* q(:,2), - 2 * b .* ...
    q(:,1) .* q(:,3), - b .* q(:,2).^2, - 2 * b .* q(:,2) .* q(:,3), - b .* q(:,3).^2]; 

x = y\log(Avox);

S0init = exp(x(1));
Dinit = mean([x(2), x(5), x(7)]);
startx(1) = S0init;
startx(2) = Dinit;
D = [ x(2), x(3), x(4);...
    x(3), x(5), x(6);...
    x(4), x(6), x(7)];

lambda = eig(D);
FA = (3/2)*sqrt((sum((lambda - startx(2)).^2)/sum(lambda)));

tic;
for i=1:1000
    
    pd = makedist('Normal');
    t2 = truncate(pd,-1,1);
    r2 = random(t2,1,1);
    startx(3) = 0.5 + 0.5 * r2(1);
    
    t3 = truncate(pd,-pi/2,pi/2);
    r3 = random(t3,2,1);
    startx(4) = r3(1);
    startx(5) = r3(2);
    
    
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
    
    Parameterit(i,:) = parameter_hat;
    Resnormit(i) = RESNORM;
    
    t = truncate(pd,0,inf);
    r = random(t,2,1);
    noise1 =  S0init + S0init * r(1);
    noise2 = Dinit + Dinit * r(2);
    
    TotIterations = TotIterations + OUTPUT.iterations;
    
    
end
AverageTime =  toc / i;
AverageIteration = TotIterations / i;
[C, I] = min(Resnormit);
BestParameters = Parameterit(I,:); 



