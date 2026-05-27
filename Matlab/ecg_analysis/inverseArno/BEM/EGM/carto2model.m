% function [newcartopoints] = carto2model(cartopoints, referencepoint, rotx, roty, rotz)
%
%
% Version 1.0: 5-aug-2015

function [newcartopoints] = carto2model(cartopoints, referencepoint, rotx, roty, rotz)

%% Rotation:

% Transform angle in degrees to angle in radians:
angle_z = rotz/180*pi;
angle_y = roty/180*pi;
angle_x = rotx/180*pi;

% Construct rotation matrices:
Rz = [cos(angle_z)   -sin(angle_z)  0
      sin(angle_z) 	 cos(angle_z)	0
      0          	 0              1];

Ry =  [cos(angle_y)  0            	sin(angle_y)
       0             1             	0
       -sin(angle_y) 0              cos(angle_y)];

Rx = [1             0                       0
      0             cos(angle_x) -sin(angle_x)
      0             sin(angle_x) cos(angle_x)];

Rot =  Rz * Ry * Rx;

% Rotate cartopoints:
rot_cartopoints = zeros(size(cartopoints,1),size(cartopoints,2));

for k = 1:size(cartopoints,1), rot_cartopoints(k,:) = (Rot*cartopoints(k,:)')'; end

%% Translation:
newcartopoints = zeros(size(cartopoints,1),size(cartopoints,2));

for k = 1:size(cartopoints,1), newcartopoints(k,:) = rot_cartopoints(k,:) + referencepoint; end
