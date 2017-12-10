
clear all; close all;

fid=fopen('dwi.Bfloat', 'r', 'b');
dwis = fread(fid, 'float');
fclose(fid);


dwis = reshape(dwis, 33, 112, 112, 50);

qhat = load('grad_dirs.txt')';

bvals = 1000*sum(qhat.*qhat);



for k=1:112
    fprintf('Line number: %i\n',k);
    for l=1:112

        
        Avox = dwis(:,k,l,25);
        if all(Avox == 0);
            AllResnorm(k,l) = 0;
            Parameters(:,k,l) = zeros(5,1);
            
        else
            
            for i=1:10
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
                parameterszob(i,:) = parameter_hat;
                
                [C, I] = min(Resnorm);
                
            end
            
            AllResnorm(k,l) = Resnorm(I);
            Parameters(:,k,l) = parameterszob(I,:);
            
        end

    end
end


figure;
subplot(2,2,1)
imshow(squeeze(Parameters(1,:,:)));
title('S0');

subplot(2,2,2)
imshow(squeeze(Parameters(2,:,:)));
title('diffusion');

subplot(2,2,3)
imshow(squeeze(Parameters(3,:,:)));
title('f');

subplot(2,2,4)
imshow(AllResnorm);
title('SSD');

figure;
Phi = squeeze(Parameters(4,:,:));
Theta = squeeze(Parameters(5,:,:));
row = (1:112);
col = (1:112);
u = (cos(Phi).* sin(Theta)) .* squeeze(Parameters(3,:,:));
v = sin(Phi).* sin(Theta) .* squeeze(Parameters(3,:,:));
title('fibre direction');

quiver(row, col, u, v);

