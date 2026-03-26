

% init model and heart distances
MODEL = readGeomPeacsModel_prof('./data/BSMUMC048Extended','BSMUMC048Extended');
HD = MODEL.VENTR.HEARTDIST;

% ECG sim vertexes
foci_1 = 275;      % left septal wall     
foci_2 = 1720;     % left free wall, at the base of one of the pupillary muscle 
foci_3 = 271;      % left free wall, at the base of one of the pupillary muscle            
foci_4 = 2329;     % right free wall at the base of the trab               
foci_5 = 668;      % right septal apical wall

% foci setup
foci     = [foci_1,foci_2,foci_3,foci_4, foci_5 ];  % foci vertexes ids
foci_pos = MODEL.VENTR.geom.VER(foci,:);            % foci vertexes
fociact  = [0, 15, 15, 15, 10];                     % foci activation times [ms]

% calculate depolarization and repolarization
dep=[];
for i=1:length(foci)
    dep(i,:) = bsxfun(@plus, HD(foci(i),:), fociact(i));
end
dep = min(dep,[],1)';
rep = 350 + mean(dep) - dep*0.6;

% time vector setup
start = 0;
ending = 1000;
maxt = ending-start;
T = ones(length(dep),1)*(0:maxt);

% old function S
% S_old = gets_v2(T,dep,rep, [0.5, 0, -0.05, -0.05, 0],4);

% fixed base parameters
p = [0.2, 0, -0.04, 0, 0];
values = [0, 0.1, 0.25, 0.5, 0.75, 0.9, 1];

% MODE
% mode = 2; % ST linear
% mode = 3; % ST cosinus
% mode = 4; % ST constatn dep and rep time, mod st angle
mode = getsMode.LinST_rept_ang_dept;


figure;
ax1 = subplot(1, 2, 1);
hold(ax1, 'on'); 

new_p = p;
for x = values(1:end)
    new_p(5) = x*0.2;
    S = gets(T,dep,rep,new_p,mode);
    plot(ax1, S(1,1:600), 'DisplayName', sprintf('p(5) = %.2f', new_p(5)));
end
% title(ax1, "Angle of ST inclination"); % mode == 4
title(ax1, "ST drop");
grid(ax1, 'on');
legend(ax1, 'show');

ax2 = subplot(1, 2, 2);
hold(ax2, 'on');
new_p = p;
for x = values(1:end)
    new_p(5) = x*0.2;
    S = gets(T,dep,rep,new_p,mode);
    plot(ax2, diff(S(1,1:600)), 'DisplayName', sprintf('p(5) = %.2f', new_p(5)));
end
grid(ax2, 'on');
legend(ax2, 'show');


sgtitle(sprintf('gets v2 function, %s', mode));
