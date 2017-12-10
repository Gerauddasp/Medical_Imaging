% test direction:
% fibdir = [cos(phi)*sin(theta) sin(phi)*sin(theta) cos(theta)];

phi = [pi/2, -pi/2];
theta = [pi/2, -pi/2];


%u = cos(phi) .* sin(theta);
%v = sin(phi) .* sin(theta);
%v = cos(theta);
u = [pi, 0];
v = [1, 1];
x = [1, 2];
y = [1, 2];

quiver(x,y,u,v);