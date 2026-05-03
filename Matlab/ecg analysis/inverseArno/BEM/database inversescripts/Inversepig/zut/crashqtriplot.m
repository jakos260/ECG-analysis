clearvars
[VER,ITRI]=loadtri('/Users/peteroosterhoff/Documents/Werk/STW/Data/geometries/Pig09/Pig09_ventricles.tri');
qtriplot(VER,ITRI,'Ventricles');
load('crashfunVentrPig09.mat','A');
% A(isnan(A))=0; % fix
qtriplot(A);