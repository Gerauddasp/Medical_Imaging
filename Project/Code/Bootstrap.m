function [s2lb, s2ub, perc95lb, perc95ub] = Bootstrap(func, data, iterations, h, startx)
% input: the function we are using, the data as a column vector, and
% iterations the number of bootstrap we are doing
% output: s2 is the 2*std interval and perc95 is the 95% confidence
% interval


echo = [29, 41, 53, 68.2, 80.2, 92.2, 107.4, 119.4, 131.4, 146.6, 158.6, 170.6]; 


A = repmat(data,1,iterations);

pd = makedist('Normal');
t = truncate(pd,0,Inf);
r = random(t,length(data),iterations);
g = truncate(pd,0.8,1.2);

A = A .* r;


for i=1:iterations
    startx(1) = mean(A(:,i));
    startx(2) = 120;
    for it=1:10
        r = random(g,1,2);
        startxb = r .* startx;
        [x1,fval,exitflag,output] = fminunc(func,startxb,h,A(:,i),echo');
        res(it,1:2) = abs(x1);
        res(it,3) = OneExp(abs(x1), A(:,i), echo');
    end
    [~, I] = min(res(:,3));
    T2(i) = res(I,2);
end

s2lb = mean(T2) - 2*std(T2);
s2ub = mean(T2) + 2*std(T2);

perc95lb = prctile(T2, 2.5);
perc95ub = prctile(T2, 97.5);

end

