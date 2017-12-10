function [res, zob] = TwoExpFix(startx, Measurements, ET)
S0IC = abs(startx(1));
S0EC = abs(startx(2));

ET = [29, 41, 53, 68.2, 80.2, 92.2, 107.4, 119.4, 131.4, 146.6, 158.6, 170.6];
ET = ET';


if size(Measurements,2) > 1
    Measurements = Measurements';
end

M = S0IC * exp(-ET/53) + S0EC * exp(-ET/133) ;

res = sum((Measurements - M).^2);

end