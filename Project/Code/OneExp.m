function [res, zob] = OneExp(startx, Measurements, ET)
S0 = abs(startx(1));
T2 = abs(startx(2));

if size(Measurements,2) > 1
    Measurements = Measurements';
end

M = S0 * exp(-ET/T2);

res = sum((Measurements - M).^2);

end

