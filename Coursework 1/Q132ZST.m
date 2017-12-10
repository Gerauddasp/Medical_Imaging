% Q132

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

BestStart = [1.0087, 1.1692e-09, 0.9080, -1.5425, 0.0105, 486.9675];

for T=1:10
    fprintf('%i\n',T)
    
    for i=1:2000
        if rem(i,100) == 0
            fprintf('Iteration number: %i\n',i);
        end
        
        pd = makedist('Normal','mu',1,'sigma',2);
        t1 = truncate(pd,-0.1,0.1);
        r = random(t1,4,1);
        startx(1) = BestStart(1) + BestStart(1) * r(1);
        startx(2) = BestStart(2) + BestStart(2) * r(2);
        startx(6) = BestStart(6) + BestStart(6) * r(3);
        
        t2 = truncate(pd,-0.1,0.1);
        r2 = random(t2,1,1);
        startx(3) = BestStart(3) + BestStart(3) * r2(1);
        
        t3 = truncate(pd,-0.1,0.1);
        r3 = random(t3,2,1);
        startx(4) = BestStart(4) + BestStart(4) * r3(1);
        startx(5) = BestStart(5) + BestStart(5) *r3(2);
        %starthist(i,:) = startx;
        
        
        
        
        h=optimset('MaxFunEvals',20000,...
            'Algorithm','levenberg-marquardt',...
            'LargeScale','off',...
            'Display','off',...
            'MaxFunEvals',20000,...
            'TolX',1e-10,...
            'TolFun',1e-10);
        % Now run the fitting
        [parameter_hat,RESNORM,EXITFLAG,OUTPUT]=fminunc('ZeppelinAndStickTortion',startx,h,meas,bvals,grad_dirs);
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
        fprintf('SSD = %d\n', RESNORMhist(T))
    end
end

% Here we try to compare the results we obtain:
% Extract the parameters
%x = real(x);
x = BestStart;
S0 = abs(x(1));
diff = abs(x(2));
f = sin(x(3))^2;
theta = x(4);
phi = x(5);
lambda1 = abs(x(6));
lambda2 = (1-f) * lambda1;

qhat = grad_dirs;
Avox = meas;
% Synthesize the signals
fibdir = [cos(phi)*sin(theta) sin(phi)*sin(theta) cos(theta)];
fibdotgrad = sum(qhat.*repmat(fibdir, [length(qhat), 1])');

S = S0*(f*exp(-bvals*diff.*(fibdotgrad.^2)) + (1-f)*exp(-bvals.*(lambda2 + (lambda1 - lambda2)*(fibdotgrad.^2))));
% Compute the sum of square differences
sumRes = sum((Avox - S').^2);

figure;
% Plot the actual data points!
plot(Avox, ' bs', 'MarkerSize', 1, 'LineWidth', 1);
hold on;
% Add the predictions to the plot!
plot(S, ' rx', 'MarkerSize', 1, 'LineWidth', 1);
% Add labels and legend.!
set(gca, 'FontName', 'Times');
set(gca, 'FontSize', 26);
xlabel('\bf{q} index');
ylabel('S');
legend('Data', 'Model');