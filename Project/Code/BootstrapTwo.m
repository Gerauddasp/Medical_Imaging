function [s2lbS0IC, s2ubS0IC, s2lbT2IC, s2ubT2IC, s2lbS0EC, s2ubS0EC, s2lbT2EC, s2ubT2EC] = BootstrapTwo(func, data, iterations, h, startx)
% input: the function we are using, the data as a column vector, and
% iterations the number of bootstrap we are doing
% output: s2 is the 2*std interval and perc95 is the 95% confidence
% interval




A = repmat(data,1,iterations);

pd = makedist('Normal');
t = truncate(pd,0.5,1.5);
r = random(t,length(data),iterations);

A = A .* r;

pd = makedist('Normal');
t = truncate(pd,0.8,1.2);



for i=1:iterations
    
    startx(1) = 1; % S0 IC
    startx(2) = 53; % T2 IC
    startx(3) = 1; % S0 EC
    startx(4) = 133; % T2 EC
    for j=1:10
        r1 = random(t,1,1);
        r2 = random(t,1,1);
        r3 = random(t,1,1);
        r4 = random(t,1,1);
        startxb(1) = startx(1)*r1;
        startxb(2) = startx(2)*r2;
        startxb(3) = startx(3)*r3;
        startxb(4) = startx(4)*r4;
        [x1,fval,exitflag,output] = fminunc(func,startxb,h,A(:,i));
        if x1(2) > x1(4)
            temp = x1(1);
            x1(1) = x1(4);
            x1(4) = temp;
        end
        res(j,1) = abs(x1(1));
        res(j,2) = abs(x1(2));
        res(j,3) = abs(x1(3));
        res(j,4) = abs(x1(4));
        res(j,5) = TwoExp(x1,A(:,i));
    end
    [zob, I] = min(res(:,5));
    S0IC(i) = abs(res(I,1));
    T2IC(i) = abs(res(I,2));
    S0EC(i) = abs(res(I,3));
    T2EC(i) = abs(res(I,4));
    
    
    
end

s2lbS0IC = mean(S0IC) - 2*std(S0IC);
s2ubS0IC = mean(S0IC) + 2*std(S0IC);

s2lbT2IC = mean(T2IC) - 2*std(T2IC);
s2ubT2IC = mean(T2IC) + 2*std(T2IC);

s2lbS0EC = mean(S0EC) - 2*std(S0EC);
s2ubS0EC = mean(S0EC) + 2*std(S0EC);

s2lbT2EC = mean(T2EC) - 2*std(T2EC);
s2ubT2EC = mean(T2EC) + 2*std(T2EC);



end