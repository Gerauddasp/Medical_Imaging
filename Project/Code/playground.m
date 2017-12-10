clear all; close all;

echo = [29, 41, 53, 68.2, 80.2, 92.2, 107.4, 119.4, 131.4, 146.6, 158.6, 170.6];

T2a = 40;
T2b = 200;
S0 = 3;

Signal = S0*exp(-echo/T2a) + S0*exp(-echo/T2b);

plot(echo, Signal);
hold on

T2a = 60;
T2b = 60;

Signal = S0*exp(-echo/T2a) + S0*exp(-echo/T2b);

plot(echo, Signal,'Color','r');