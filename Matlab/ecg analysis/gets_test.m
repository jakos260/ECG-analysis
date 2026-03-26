% -------------------------------------------------------------------------
% ----------------------- Testing gets_v2() -------------------------------
% -------------------------------------------------------------------------

close all; clear all;

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
start = 700;
ending = 1800;
maxt = ending-start;
T = ones(length(dep),1)*(0:maxt);

% old function S
% S_old = gets_v2(T,dep,rep, [0.5, 0, -0.05, -0.05, 0],4);

figure;

% fixed base parameters
p = [0.5, 0.2, 0.1, 0.0, 0.0];

params_to_iterate = {
    3, 'ST drop',   [0, 0.2, 0.6]; % ampl = [1.0 : 0.8] -> p(3) = [0: 1]
    4, 'ST dv',     [-0.5, 0, 0.6]; % dv = [-1, 1] x1000 [-1/s, 1/s] -> p(4) = [-1, 1]
    5, 'rep dv',    [0.1, 0.5, 0.9]
};
[n_rows, ~] = size(params_to_iterate);

mode = getsMode.Two_dy_Spline;

for k = 1:n_rows
    ax1 = subplot(2, n_rows, k);
    hold(ax1, 'on'); 
    
    param           = params_to_iterate{k,1};
    param_name      = params_to_iterate{k,2};
    param_values    = params_to_iterate{k,3};
    
    new_p = p;
    for x = 1:length(param_values)
        new_p(param) = param_values(x);
        S = gets(T,dep,rep,new_p,mode);

        plot(ax1, S(1,1:600), 'DisplayName', sprintf('p(%d) = %.2f', param, new_p(param)));
    end
    title(ax1, param_name);
    grid(ax1, 'on');
    legend(ax1, 'show', 'Location', 'southwest');

    ax2 = subplot(2, n_rows, k+n_rows);
    hold(ax2, 'on');

    new_p = p;
    for x = 1:length(param_values)
        new_p(param) = param_values(x);
        S = gets(T,dep,rep,new_p,mode);

        plot(ax2, diff(S(1,1:600)), 'DisplayName', sprintf('diff p(%d) = %.2f', param, new_p(param)));
    end
    grid(ax2, 'on');
    legend(ax2, 'show', 'Location', 'northeast');

end

sgtitle(sprintf('gets v2 function, mode = %s', mode));
