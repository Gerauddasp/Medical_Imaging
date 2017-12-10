function [res, zob] = TwoExp(startx, Measurements)
S0IC = abs(startx(1));
T2IC = abs(startx(2));
S0EC = abs(startx(3));
T2EC = abs(startx(4));

ET = [29, 41, 53, 68.2, 80.2, 92.2, 107.4, 119.4, 131.4, 146.6, 158.6, 170.6];
ET = ET';
if T2IC > T2EC
    temp = T2IC;
    T2IC = T2EC;
    T2EC = temp;
end

if size(Measurements,2) > 1
    Measurements = Measurements';
end

M = S0IC .* exp(-ET/T2IC) + S0EC .* exp(-ET/T2EC) ;

res = sum((Measurements - M).^2);

end