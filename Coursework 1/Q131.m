% Q131
close all; clear all;

datafile = 'challengeOPEN.txt'; fid = fopen(datafile, 'r', 'b');
% Read in the header
A = fgetl(fid);
% Read in the data
A = fscanf(fid, '%f', [8, inf]); fclose(fid);
% Create the protocol
meas = A(1,:)';
grad_dirs = A(2:4,:);
G = A(5,:);
delta = A(6,:);
smalldel = A(7,:);
TE = A(8,:);
GAMMA = 2.675987E8;
bvals = ((GAMMA*smalldel.*G).^2).*(delta-smalldel/3);

BestStart = [0.9791, 1.3491e-09, 0.7709, 1.5968, -0.0287];
zob = 0;
for T=1:10
    
    fprintf('%i\n',T);
    
    for i=1:1000
        
        pd = makedist('Normal','mu',1,'sigma',1);
        t1 = truncate(pd,-0.1,0.1);
        r = random(t1,2,1);
        startx(1) = BestStart(1) + BestStart(1) * r(1);
        startx(2) = BestStart(2) + BestStart(2) * r(2);
        
        t2 = truncate(pd,-0.1,0.1);
        r2 = random(t2,1,1);
        startx(3) = BestStart(3) + BestStart(3) * r2(1);
        
        t3 = truncate(pd,-0.1,0.1);
        r3 = random(t3,2,1);
        startx(4) = BestStart(4) + BestStart(4) * r3(1);
        startx(5) = BestStart(5) + BestStart(5) * r3(2);
        %starthist(i,:) = startx;
        
        
        h=optimset('MaxFunEvals',20000,...
            'Algorithm','levenberg-marquardt',...
            'LargeScale','off',...
            'Display','off',...
            'MaxFunEvals',20000,...
            'TolX',1E-10,...
            'TolFun',1E-10);
        
        % Now run the fitting
        [parameter_hat,RESNORM,EXITFLAG,OUTPUT]=fminunc('BallStickSSD',startx,h,meas,bvals,grad_dirs);
        parameter_hat(1) = abs(parameter_hat(1));
        parameter_hat(2) = abs(parameter_hat(2));
        parameter_hat(3) = sin(parameter_hat(3))^2;
        
        resnorm(i) = RESNORM;
        Parameters(i,:) = parameter_hat;
        
        
    end
    
    [C, I] = min(resnorm);
    RESNORMhist(T) = resnorm(I);
    
    if RESNORMhist(T) == min(RESNORMhist)
        BestStart = Parameters(I,:);
        BestStart(3) = asin(sqrt(BestStart(3)));
        fprintf('SSD: %d\n',RESNORMhist(T));
        zob = zob+1;
        BestStarthist(zob,:) = BestStart;
    end
end

% Here we try to compare the results we obtain:
qhat = grad_dirs;
Avox = meas;

S0 = abs(x(1));
diff = abs(x(2));
f = sin(x(3))^2;
theta = x(4);
phi = x(5);

% Synthesize the signals
fibdir = [cos(phi)*sin(theta) sin(phi)*sin(theta) cos(theta)];
fibdotgrad = sum(qhat.*repmat(fibdir, [length(qhat), 1])');
model = S0*(f*exp(-bvals*diff.*(fibdotgrad.^2)) + (1-f)*exp(-bvals*diff));

figure;
% Plot the actual data points!
plot(Avox, ' bs', 'MarkerSize', 1, 'LineWidth', 2);
hold on;
% Add the predictions to the plot!
plot(model, ' rx', 'MarkerSize', 1, 'LineWidth', 2);
% Add labels and legend.!
set(gca, 'FontName', 'Times');
set(gca, 'FontSize', 26);
xlabel('\bf{q} index');
ylabel('S');
legend('Data', 'Model');

SSD = sum((Avox - model').^2);

