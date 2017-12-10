%% Q111
close all; clear all;

fid = fopen('dwi.Bfloat', 'r', 'b');
dwis = fread(fid, 'float');
fclose(fid);


dwis = reshape(dwis, 33, 112, 112, 50);

qhat = load('grad_dirs.txt')';

bvals = 1000*sum(qhat.*qhat);

Avox = dwis(:,52,62,25);

% Middle slice of the first image volume, which has b=0 
%imshow(squeeze(dwis(1,:,:,25)), []);

% Middle slice of the second image volume with b=1000 
imshow(squeeze(dwis(1,:,:,25)), []);

% Define a starting point for the non-linear fit

startx = [250000, 1E-3, 0.5, 0, 0];


h=optimset('MaxFunEvals',20000,...
   'Algorithm','levenberg-marquardt',...
   'LargeScale','off',...
   'Display','off',...
   'TolX',1e-1000,...
   'TolFun',1e-1000);
% Now run the fitting
[parameter_hat,RESNORM,EXITFLAG,OUTPUT]=fminunc('BallStickSSDone',startx,h,Avox,bvals,qhat);

% Here we try to compare the results we obtain:
x = parameter_hat;
S0 = x(1);
diff = x(2);
f = x(3);
theta = x(4);
phi = x(5);

% Synthesize the signals
fibdir = [cos(phi)*sin(theta) sin(phi)*sin(theta) cos(theta)];
fibdotgrad = sum(qhat.*repmat(fibdir, [length(qhat), 1])');
model = S0*(f*exp(-bvals*diff.*(fibdotgrad.^2)) + (1-f)*exp(-bvals*diff));

figure;
% Plot the actual data points!
plot(Avox, ' bs', 'MarkerSize', 16, 'LineWidth', 4);
hold on;
% Add the predictions to the plot!
plot(model, ' rx', 'MarkerSize', 16, 'LineWidth', 4);
% Add labels and legend.!
set(gca, 'FontName', 'Times');
set(gca, 'FontSize', 26);
xlabel('\bf{q} index');
ylabel('S');
legend('Data', 'Model');
%% Q112

clear all; close all;

fid=fopen('dwi.Bfloat', 'r', 'b');
dwis = fread(fid, 'float');
fclose(fid);


dwis = reshape(dwis, 33, 112, 112, 50);

qhat = load('grad_dirs.txt')';

bvals = 1000*sum(qhat.*qhat);

Avox = dwis(:,52,62,25);

% Define a starting point for the non-linear fit

startx = [250000, 1E-3, 0.5, 0, 0];
startx(3) = asin(sqrt(startx(3)));

h=optimset('MaxFunEvals',20000,...
   'Algorithm','levenberg-marquardt',...
   'LargeScale','off',...
   'Display','off',...
   'TolX',1e-1000,...
   'TolFun',1e-1000);

% Now run the fitting
[parameter_hat,RESNORM,EXITFLAG,OUTPUT]=fminunc('BallStickSSD',startx,h,Avox,bvals,qhat);
parameter_hat(1) = abs(parameter_hat(1));
parameter_hat(2) = abs(parameter_hat(2));
parameter_hat(3) = sin(parameter_hat(3))^2;

% Here we try to compare the results we obtain:
x = parameter_hat;
S0 = x(1);
diff = x(2);
f = x(3);
theta = x(4);
phi = x(5);

% Synthesize the signals
fibdir = [cos(phi)*sin(theta) sin(phi)*sin(theta) cos(theta)];
fibdotgrad = sum(qhat.*repmat(fibdir, [length(qhat), 1])');
model = S0*(f*exp(-bvals*diff.*(fibdotgrad.^2)) + (1-f)*exp(-bvals*diff));

figure;
% Plot the actual data points!
plot(Avox, ' bs', 'MarkerSize', 16, 'LineWidth', 4);
hold on;
% Add the predictions to the plot!
plot(model, ' rx', 'MarkerSize', 16, 'LineWidth', 4);
% Add labels and legend.!
set(gca, 'FontName', 'Times');
set(gca, 'FontSize', 26);
xlabel('\bf{q} index');
ylabel('S');
legend('Data', 'Model');

%% Q113

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
    
end


hist(Resnorm,100);
set(gca, 'FontSize', 14);
xlabel('SSD');
ylabel('# of solutions');

%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Question 4: Train part
clear all; close all;

for k=1:112
    for l=1:112
        
        i=0;
        n = 0;
        totn = 1;
        
        while n/totn < 0.4 && i<10
            i=i+1;
            
            fid=fopen('dwi.Bfloat', 'r', 'b');
            dwis = fread(fid, 'float');
            fclose(fid);
            
            
            dwis = reshape(dwis, 33, 112, 112, 50);
            
            qhat = load('grad_dirs.txt')';
            
            bvals = 1000*sum(qhat.*qhat);
            
            Avox = dwis(:,k,l,25);
            
            % Define a starting point for the non-linear fit
            %             noise1 = 250000 + 250000 * randn;
            %             noise2 = 1E-5 + (0 - (1E-5)).* rand;
            %             noise3 = rand;
            %             noise4 = -pi/2 + (pi/2 - (-pi/2)).* rand;
            %             noise5 = -pi/2 + (pi/2 - (-pi/2)).* rand;
            pd = makedist('Normal');
            t = truncate(pd,0,inf);
            r = random(t,2,1);
            noise1 = 2.5E5 + (2.5E5) * r(1);
            noise2 = 1E-3 + (1E-3) * r(2);
            
            t2 = truncate(pd,-1,1);
            r2 = random(t,1,1);
            noise3 = 0.5 + 0.5 * r2(1);
            
            t3 = truncate(pd,-pi/2,pi/2);
            r3 = random(t,2,1);
            noise4 = r3(1);
            noise5 = r3(2);
            
            startx = [noise1, noise2, noise3, noise4, noise5];
            
            h=optimset('MaxFunEvals',20000,...
                'Algorithm','levenberg-marquardt',...
                'LargeScale','off',...
                'MaxIter', 20000,...
                'Display','off',...
                'TolX',1e-100,...
                'TolFun',1e-100);
            % Now run the fitting
            [parameter_hat,RESNORM,EXITFLAG,OUTPUT]=fminunc('BallStickSSD',startx,h,Avox,bvals,qhat);
            parameter_hat(1) = abs(parameter_hat(1));
            parameter_hat(2) = abs(parameter_hat(2));
            parameter_hat(3) = exp(-(parameter_hat(3)^2));
            
            Resnorm(i)= RESNORM;
            n = sum(Resnorm==mode(Resnorm));
            totn = length(Resnorm);
            parameterszob(i,:) = parameter_hat;
            
            [C, I] = min(Resnorm);
            
        end
        
        Parameters(:,k,l) = parameter_hat(I,:);

    end
end

save('ParametersVox25',Parameters);



%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Question 4: Plot part.
Parameters = load('ParametersVox25.mat','-mat');
Parameters = Parameters.Parameters;



%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Question 1.2.1
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

while n/totn < 0.4 && i<100
    i=i+1;
    
    if rem(i,10)==0
        fprintf('i = %i\n',i);
    end    
    % Define a starting point for the non-linear fit
%     noise1 = (1E6) * rand;
%     noise2 = (1E-1)* rand;
%     noise3 = rand;
%     noise4 = -pi/2 + (pi/2 - (-pi/2)).* rand;
%     noise5 = -pi/2 + (pi/2 - (-pi/2)).* rand;
    pd = makedist('Normal');
    t = truncate(pd,0,inf);
    r = random(t,2,1);
    noise1 = 2.5E5 + (2.5E5) * r(1);
    noise2 = 1E-3 + (1E-3) * r(2);
    
    t2 = truncate(pd,-1,1);
    r2 = random(t,1,1);
    noise3 = 0.5 + 0.5 * r2(1);
    
    t3 = truncate(pd,-pi/2,pi/2);
    r3 = random(t,2,1);   
    noise4 = r3(1);
    noise5 = r3(2);    
    startx = [noise1, noise2, noise3, noise4, noise5];
    
    h=optimset('MaxFunEvals',20000,...
        'Algorithm','levenberg-marquardt',...
        'LargeScale','off',...
        'MaxIter', 20000,...
        'Display','off',...
        'TolX',1e-100,...
        'TolFun',1e-100);
    % Now run the fitting
    [parameter_hat,RESNORM,EXITFLAG,OUTPUT]=fminunc('BallStickSSD',startx,h,Avox,bvals,qhat);
    parameter_hat(1) = abs(parameter_hat(1));
    parameter_hat(2) = abs(parameter_hat(2));
    parameter_hat(3) = exp(-(parameter_hat(3)^2));
    
    Resnorm(i)= RESNORM;
    n = sum(Resnorm==mode(Resnorm));
    totn = length(Resnorm);
    parameterszob(i,:) = parameter_hat;
    
    [C, I] = min(Resnorm);
    
end


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
K = 100;
newA = model';
A_hat = repmat(newA,1,K) + randn(size(newA,1),K) .* 1;

% Now estimates best fit for all A_hat:
fprintf('Now estimates best fit for all A_hat\n')

for k=1:K
    
    if rem(k,10) == 0
        fprintf('We are in bootstrap number: %f\n',k);
    end
    
    Avox = A_hat(:,k);
    
    i=0;
    n = 0;
    totn = 1;
    
    while n/totn < 0.4 && i<20
        i=i+1;
        
        
        % Define a starting point for the non-linear fit
        %         noise1 = (1E9) * rand;
        %         noise2 = (1E-1)* rand;
        %         noise3 = rand;
        %         noise4 = -pi/2 + (pi/2 - (-pi/2)).* rand;
        %         noise5 = -pi/2 + (pi/2 - (-pi/2)).* rand;
        pd = makedist('Normal');
        t = truncate(pd,0,inf);
        r = random(t,2,1);
        noise1 = x(1) + x(1) * r(1);
        noise2 = x(2) + x(2) * r(2);
        
        t2 = truncate(pd,-1,1);
        r2 = random(t,1,1);
        noise3 = x(3) + x(3) * r2(1);
        
        t3 = truncate(pd,-pi/2,pi/2);
        r3 = random(t,2,1);
        noise4 = x(3) + x(3) * r3(1);
        noise5 = x(3) + x(3) * r3(2);
        startx = [noise1, noise2, noise3, noise4, noise5];
        
        h=optimset('MaxFunEvals',20000,...
            'Algorithm','levenberg-marquardt',...
            'LargeScale','off',...
            'MaxIter', 20000,...
            'Display','off',...
            'TolX',1e-100,...
            'TolFun',1e-100);
        % Now run the fitting
        [parameter_hat,RESNORM,EXITFLAG,OUTPUT]=fminunc('BallStickSSD',startx,h,Avox,bvals,qhat);
        parameter_hat(1) = abs(parameter_hat(1));
        parameter_hat(2) = abs(parameter_hat(2));
        parameter_hat(3) = exp(-(parameter_hat(3)^2));
        
        Resnorm(i)= RESNORM;
        n = sum(Resnorm==mode(Resnorm));
        totn = length(Resnorm);
        parameterszob(i,:) = parameter_hat;
        
        [C, I] = min(Resnorm);
        
    end
    Parameters_boot(k,:) = parameterszob(I,:);
    R(k) = Resnorm(I);
    
end

%%
%%%%%%%%%%%%%%%%%%%%%%
% Question 1.2.2


