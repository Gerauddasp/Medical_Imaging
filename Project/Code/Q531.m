clear all; 
close all;
RUN = false;
subject = 2;
PLD = [50, 50, 50, 50, 300, 700, 1200];
TD = [500, 750, 1000, 3000, 3000, 3000, 3000];
echo = [29, 41, 53, 68.2, 80.2, 92.2, 107.4, 119.4, 131.4, 146.6, 158.6, 170.6]; 
% Import data from question Q1.2

%load('T2vsTE.mat');

%% treat data
load('maskpoly.mat');
load('data');

% taking the ASL signal of the ROI defined in Q52 with diffusion gradients
% for vascular suppression on
T2vsTE(1:7) = struct('ASLon', [], 'ASLoff', [], 'CONTROLoff', [], 'CONTROLon', []);
for subject=1:9
    for e=1:12
        for control=1:7
            mat = mask(subject).BW .* data(subject).perf_w_on(:,:,e,control);
            T2vsTE(control).ASLon(e,subject) = mean(mat(mat~=0));
        end
    end
end

for subject=1:9
    for e=1:12
        for control=1:7
            mat = mask(subject).BW .* data(subject).perf_w_off(:,:,e,control);
            T2vsTE(control).ASLoff(e,subject) = mean(mat(mat~=0));
        end
    end
end

% taking the Control signal of the ROI defined in Q52 with diffusion gradients
% for vascular suppression on

for subject=1:9
    for e=1:12
        for control=2:2:14
            mat = mask(subject).BW .* data(subject).new_images_on(:,:,e,control);
            T2vsTE(control/2).CONTROLon(e,subject) = mean(mat(mat~=0));
        end
    end
end


for subject=1:9
    for e=1:12
        for control=2:2:14
            mat = mask(subject).BW .* data(subject).new_images_off(:,:,e,control);
            T2vsTE(control/2).CONTROLoff(e,subject) = mean(mat(mat~=0));
        end
    end
end
if RUN
%% find the best T2 value with one exponential
h=optimset('MaxFunEvals',200000,...
    'Algorithm','levenberg-marquardt',...
    'LargeScale','off',...
    'Display','off',...
    'TolX',1e-10,...
    'TolFun',1e-1);
h2=optimset('MaxFunEvals',200000,...
    'Algorithm','levenberg-marquardt',...
    'LargeScale','off',...
    'Display','off',...
    'TolX',1e-10,...
    'TolFun',1e-10);

Fit1 = struct('ASLon', [], 'ASLoff', [], 'CONTROLon' , [], 'CONTROLoff', []);
Conf1 = struct('ASLon', [], 'ASLoff', [], 'CONTROLon' , [], 'CONTROLoff', []);
% Here we find the T2 value for each subject and for each configuration
% with it's confidence interval and this for gradient off and on.

pd = makedist('Normal');
t = truncate(pd,0,Inf);
for control=1:7
    for subject=1:9
        startx(1) = mean(T2vsTE(control).ASLoff(:,subject));
        fprintf('Control: %i and Subject: %i\n', control, subject);
        startx(2) = 120;
        
        for i=1:10
            r = random(t,1,2);
            startxn = startx .* r;
            [x1,fval,exitflag,output] = fminunc('OneExp',...
                startxn,h,T2vsTE(control).ASLoff(:,subject),echo');
            res(i,1) = abs(x1(1));
            res(i,2) = abs(x1(2));
            res(i,3) = OneExp(x1,T2vsTE(control).ASLoff(:,subject),echo');
        end
        [zob, I] = min(res(i,3));
        Fit1.ASLoff(subject,control) = abs(res(I,2));
        % saving results for AIC
        ErrorASLoff(subject, control) = res(I,3);
        startBS = res(I,1:2);
        [Conf1.ASLoff(subject,control,1), Conf1.ASLoff(subject,control,2),...
            Conf1.ASLoff(subject,control,3), Conf1.ASLoff(subject,control,4)] =...
            Bootstrap('OneExp', T2vsTE(control).ASLoff(:,subject),100,h ,startx);
        
        startx(1) = mean(T2vsTE(control).ASLon(:,subject));
        startx(2) = 120;
        for i=1:10
            r = random(t,1,2);
            startxn = startx .* r;
            [x2, fval, exitflag, output] = fminunc('OneExp',...
                startxn, h, T2vsTE(control).ASLon(:,subject),echo');
            res(i,1) = abs(x2(1));
            res(i,2) = abs(x2(2));
            res(i,3) = OneExp(x2,T2vsTE(control).ASLon(:,subject),echo');
        end
        [zob, I] = min(res(i,3));
        Fit1.ASLon(subject,control) = abs(res(I,2));
        % saving results for AIC
        ErrorASLon(subject, control) = res(I,3);
        %AICerror(subject, control) = res(I,3);
        %startBS = res(I,1:2);
        [Conf1.ASLon(subject,control,1), Conf1.ASLon(subject,control,2),...
            Conf1.ASLon(subject,control,3), Conf1.ASLon(subject,control,4)] =...
            Bootstrap('OneExp', T2vsTE(control).ASLon(:,subject),100,h , startx);
        
        
        startx(1) = mean(T2vsTE(control).CONTROLoff(:,subject));
        startx(2) = 120;
        for i=1:10
            r = random(t,1,2);
            startxn = startx .* r;
            [x3, fval, exitflag, output] = fminunc('OneExp', startxn, h2,...
                T2vsTE(control).CONTROLoff(:,subject),echo');
            res(i,1) = abs(x2(1));
            res(i,2) = abs(x2(2));
            res(i,3) = OneExp(x2,T2vsTE(control).CONTROLoff(:,subject),echo');
        end
        [zob, I] = min(res(i,3));
        Fit1.CONTROLoff(subject, control) = abs(res(I,2));
        % saving results for AIC
        ErrorCONTROLoff(subject, control) = res(I,3);
        
        startBS = res(I,1:2);
        [Conf1.CONTROLoff(subject,control,1), Conf1.CONTROLoff(subject,control,2),...
            Conf1.CONTROLoff(subject,control,3), Conf1.CONTROLoff(subject,control,4)] =...
            Bootstrap('OneExp', T2vsTE(control).CONTROLoff(:,subject),100,h2 , startBS);
        
        startx(1) = mean(T2vsTE(control).CONTROLon(:,subject));
        startx(2) = 120;
        for i=1:10
            r = random(t,1,2);
            startxn = startx .* r;
            [x4, fval, exitflag, output] = fminunc('OneExp', startxn, h2,...
                T2vsTE(control).CONTROLon(:,subject),echo');
            res(i,1) = abs(x2(1));
            res(i,2) = abs(x2(2));
            res(i,3) = OneExp(x2,T2vsTE(control).CONTROLon(:,subject),echo');
        end
        [zob, I] = min(res(i,3));
        Fit1.CONTROLon(subject, control) = abs(res(I,2));
        % saving results for AIC
        ErrorCONTROLon(subject, control) = res(I,3);
        %startBS = res(I,1:2);
        [Conf1.CONTROLon(subject,control,1), Conf1.CONTROLon(subject,control,2),...
            Conf1.CONTROLon(subject,control,3), Conf1.CONTROLon(subject,control,4)] =...
            Bootstrap('OneExp', T2vsTE(control).CONTROLon(:,subject),100,h2 , startx);
        
    end
end
save('Q531.mat','Fit1', 'Conf1', 'ErrorASLon', 'ErrorCONTROLoff', 'ErrorCONTROLon', 'ErrorASLoff');
else
load('Q531.mat')
end
%save('AICerror.mat','AICerror');
%% Now we need to write the Table into a csv file
indice = (1:3:7*3);
for control=1:7
    for subject=1:9
        ASLoffmat(indice(control) + 1,subject) = Fit1.ASLoff(subject,control);
        ASLonmat(indice(control) + 1,subject) = Fit1.ASLon(subject, control);
        CONTROLoffmat(indice(control) + 1, subject) = Fit1.CONTROLoff(subject,control);
        CONTROLonmat(indice(control) + 1,subject) = Fit1.CONTROLon(subject, control);
        
%         ASLoffmat(indice,subject) = Conf1.ASLoff(subject, control,1);
%         ASLonmat(indice, subject) = Conf1.ASLon(subject, control,1);
%         CONTROLoffmat(indice, subject) = Conf1.CONTROLoff(subject, control,1);
%         CONTROLonmat(indice, subject) = Conf1.CONTROLon(subject, control,1);
        
        ASLoffmat(indice(control),subject) = Conf1.ASLoff(subject, control,3);
        ASLonmat(indice(control), subject) = Conf1.ASLon(subject, control,3);
        CONTROLoffmat(indice(control), subject) = Conf1.CONTROLoff(subject, control,3);
        CONTROLonmat(indice(control), subject) = Conf1.CONTROLon(subject, control,3);
        
        ASLoffmat(indice(control) + 2,subject) = Conf1.ASLoff(subject, control,4);
        ASLonmat(indice(control) + 2, subject) = Conf1.ASLon(subject, control,4);
        CONTROLoffmat(indice(control) + 2, subject) = Conf1.CONTROLoff(subject, control,4);
        CONTROLonmat(indice(control) + 2, subject) = Conf1.CONTROLon(subject, control,4);
        
%         ASLoffmat(indice + 4,subject) = Conf1.ASLoff(subject, control,2);
%         ASLonmat(indice + 4, subject) = Conf1.ASLon(subject, control,2);
%         CONTROLoffmat(indice + 4, subject) = Conf1.CONTROLoff(subject, control,2);
%         CONTROLonmat(indice + 4, subject) = Conf1.CONTROLon(subject, control,2);
    end
end
csvwrite('Q53ASLoffTable.csv',ASLoffmat);
csvwrite('Q53ASLonTable.csv',ASLonmat);
csvwrite('Q53CONTROLoffTable.csv', CONTROLoffmat);
csvwrite('Q53CONTROLonTable.csv', CONTROLonmat);
      

% Now we calculate the AIC for the oneExp here above in our ROI
%AIConeExp = 4 + 12 * log((1/12) * AICerror);

%fprintf('\nThe AIC for the One-exponential model is: %i\n',AIConeExp);