%% check graphdist:
% This script is used to verify the calculations of graphdist.c in
% different modes. There are several parts

clear all
close all
clc

%% Part 1; create simple mesh with two triangles:
%

TRIANGLE.VER = [    0 0  0;
                    0 1  0;
                    0 1  1;
                    1 1  -1];

TRIANGLE.ITRI = [1 2 3; 4 3 2];

% run graphdist mode 4 [over the surface]:
[TRIANGLE.ADJ,TRIANGLE.DIST,TRIANGLE.PATH] = graphdist(TRIANGLE.ITRI,TRIANGLE.VER,4);

% save results:
save('/Users/arnojanssen/Documents/stw/BEM/check_graphdist/triangle.mat', 'TRIANGLE');

% visualize results:
qtriplot(TRIANGLE.VER,TRIANGLE.ITRI);
qtriplot('edge on');
qtriplot('back both');
qtriplot(TRIANGLE.ADJ(1,:)');
qtriplot('bgdcolor white');

pause

% save figure:
qtriplot('png /Users/arnojanssen/Documents/stw/BEM/check_graphdist/triangle.png 1000 1000'); pause(0.2);

qtriplot('exit');

% The calculation is done over the surface: correct!

%% Part 2; create mesh with multiple second neighbours:
% This is a check to verify if all second neighbour nodes are included:

PLANE.VER = [ 2 0 0;    % 1
            4 0 0;      % 2
            6 0 0;      % 3
            1 2 0;      % 4
            3 2 0;      % 5
            5 2 0;      % 6
            7 2 0;      % 7
            0 4 0;      % 8
            2 4 0;      % 9
            4 4 0;      % 10
            6 4 0;      % 11
            8 4 0;      % 12
            1 6 0;      % 13
            3 6 0;      % 14
            5 6 0;      % 15
            7 6 0;      % 16
            2 7 0;      % 17
            4 8 0;      % 18
            6 8 0;];    % 19

PLANE.ITRI = [1 4 5;
            1 5 2;
            2 5 6;
            2 6 3;
            3 6 7;
            4 8 9;
            4 9 5;
            5 9 10;
            5 10 6;
            6 10 11;
            6 11 7;
            7 11 12;
            8 13 9;
            9 13 14;
            9 14 10;
            10 14 15;
            10 15 11;
            11 15 16;
            11 16 12;
            13 17 14;
            14 17 18;
            14 18 15;
            15 18 19;
            15 19 16];

% run graphdist mode 4 [over the surface]:
[PLANE.ADJ,PLANE.DIST,PLANE.PATH] = graphdist(PLANE.ITRI,PLANE.VER,4);

% save results:
save('/Users/arnojanssen/Documents/stw/BEM/check_graphdist/plane.mat', 'PLANE');

% visualize results:
qtriplot(PLANE.VER,PLANE.ITRI);
qtriplot('edge on');
qtriplot('back both');
qtriplot(PLANE.ADJ(10,:)');
qtriplot('bgdcolor white');

pause

% save figure:
qtriplot('png /Users/arnojanssen/Documents/stw/BEM/check_graphdist/plane_2nd.png 1000 1000'); pause(0.2);

qtriplot('exit');

% Why are not all second order neighbours taken into account: incorrect?

%% Part 3; create mesh with multiple second neighbours:
% This is a check to verify if all second neighbour nodes are included in 3D:

BEAM.VER = [0 0 0;      % 1
            1 0 0;      % 2
            1 1 0;      % 3
            0 1 0;      % 4
            0 0 1;      % 5
            1 0 1;      % 6
            1 1 1;      % 7
            0 1 1;      % 8
            0 0 2;      % 9
            1 0 2;      % 10
            1 1 2;      % 11
            0 1 2];     % 12

BEAM.ITRI = [1 6 2;
             1 5 6;
             1 3 2;
             1 4 3;
             1 5 8;
             1 8 4;
             2 6 3;
             6 7 3;
             3 7 8;
             3 8 4;
             5 9 6;
             9 10 6;
             5 12 8;
             5 9 12;
             6 11 7;
             10 11 6;
             11 8 7;
             11 12 8;
             9 11 10;
             9 12 11;];

% run graphdist mode 3:
[BEAM.ADJ,BEAM.DIST,BEAM.PATH] = graphdist(BEAM.ITRI,BEAM.VER,3);

% save results:
save('/Users/arnojanssen/Documents/stw/BEM/check_graphdist/beam.mat', 'BEAM');

% visualize results:
qtriplot(BEAM.VER,BEAM.ITRI);
qtriplot('edge on');
qtriplot('back both');
qtriplot(BEAM.ADJ(1,:)');
qtriplot('bgdcolor white');

pause

% save figure:
qtriplot('png /Users/arnojanssen/Documents/stw/BEM/check_graphdist/BEAM.png 1000 1000'); pause(0.2);

qtriplot('exit');


%% Part 4; load in a complicated [heart] mesh to calculate the distance matrices in 2D and in 3D
% 

load('/Users/arnojanssen/Documents/stw/BEM/check_graphdist/heart.mat');

% run graphdist mode 1 [over the surface]:
[heart.ADJ2D,heart.DIST2D,heart.PATH2D] = graphdist(heart.ITRI,heart.VER,4);

% visualize results:
qtriplot(heart.VER,heart.ITRI);
qtriplot('edge on');
qtriplot('back both');
qtriplot(heart.ADJ2D(100,:)');
qtriplot('bgdcolor white');
qtriplot('scale 50');

pause

% save figure:
qtriplot('png /Users/arnojanssen/Documents/stw/BEM/check_graphdist/heart_2D.png 1000 1000'); pause(0.2);

qtriplot('exit');

% run graphdist mode 2 [also through walls]:
[heart.ADJ3D,heart.DIST3D,heart.PATH3D] = graphdist(heart.ITRI,heart.VER,3);

% visualize results:
qtriplot(heart.VER,heart.ITRI);
qtriplot('edge on');
qtriplot('back both');
qtriplot(heart.ADJ3D(100,:)');
qtriplot('bgdcolor white');
qtriplot('scale 50');

pause

% save figure:
qtriplot('png /Users/arnojanssen/Documents/stw/BEM/check_graphdist/heart_3D.png 1000 1000'); pause(0.2);

qtriplot('exit');

save('/Users/arnojanssen/Documents/stw/BEM/check_graphdist/heart.mat', 'heart');
