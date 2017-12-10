function [s2lbS0IC, s2ubS0IC, s2lbS0EC, s2ubS0EC] = BootstrapTwoFix(func, data, iterations, h, startx)
% input: the function we are using, the data as a column vector, and
% iterations the number of bootstrap we are doing
% output: s2 is the 2*std interval and perc95 is the 95% confidence
% interval


echo = [29, 41, 53, 68.2, 80.2, 92.2, 107.4, 119.4, 131.4, 146.6, 158.6, 170.6]; 


A = repmat(data,1,iterations);

pd = makedist('Normal');
t = truncate(pd,0.5,1.5);
r = random(t,length(data),iterations);

A = A .* r;

pd = makedist('Normal');
t = truncate(pd,0.8,1.2);



for i=1:iterations
    
    for j=1:10
        r1 = random(t,1,1);
        r2 = random(t,1,1);
        startxb(1) = startx(1)*r1;
        startxb(2) = startx(2)*r2;
        [x1,fval,exitflag,output] = fminunc(func,startxb,h,A(:,i));
        res(j,1) = abs(x1(1));
        res(j,2) = abs(x1(2));
        res(j,5) = TwoExpFix(x1,A(:,i),echo');
    end
    [zob, I] = min(res(:,5));
    S0IC(i) = abs(res(I,1));
    S0EC(i) = abs(res(I,2));
    
    
    
end

s2lbS0IC = mean(S0IC) - 2*std(S0IC);
s2ubS0IC = mean(S0IC) + 2*std(S0IC);


s2lbS0EC = mean(S0EC) - 2*std(S0EC);
s2ubS0EC = mean(S0EC) + 2*std(S0EC);




end