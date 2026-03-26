MAT = load('data\BSMUMC048_CARTOMAP_Allinfo\BSMUMC048_CARTOMAP_Allinfo.mat');
MODEL = readGeomPeacsModel_prof('data\PPD1_ECGSIM_NEW','PPD1_ECGSIM_NEW');
M = MAT.Map2Copy2;


%%  Signal computation

% mode = 2; % ST section linear
% mode = 3; % ST section cosinus
% mode = 4; % ST section linear, downslope position correction
mode = 5;

refPos = MODEL.GEOM.thorax.VER(38,:) - MODEL.GEOM.thorax.NORMV(38,:) * 20;
HD = MODEL.VENTR.HEARTDIST;

%ECG sim vertexes
foci_1 = 263;      % left septal wall     
foci_2 = 201;      % left septal wall
foci_3 = 214;      % right septal wall              
foci_4 = 95;       % right free wall                


% FOCI AND FOCI ACT
foci        = [foci_1,foci_2,foci_3,foci_4 ];   
fociact     = [6, 8, 13, 22];   
foci_pos    = MODEL.VENTR.geom.VER(foci,:);


dep=[];
for i=1:length(foci)
    dep(i,:) = bsxfun(@plus, HD(foci(i),:), fociact(i));
end
dep = min(dep,[],1)';
rep = 350 + mean(dep) - dep*0.6;

% compute the transmembrane potentials
start = 700;
ending = 1800;
maxt = ending-start;
T = ones(length(dep),1)*(0:maxt);
p(1) = 2;
p(2) = 0;
p(3) = -0.025;
p(4) = 0;
p(5) = 0.2;
S_nodes= gets_v2(T,dep,rep,p,5); % old mode
lp=5;
doplot = 0;

%% measured signal  

nV = size(MODEL.VENTR.VENTRICLES,1);
chan_idx = [1:8,19:38];
chan_names = M.channels(chan_idx);
A = M.ECG_all;
n_CH = length(chan_idx);

BestR = -inf(1,n_CH);
Corr = NaN(nV, n_CH); 
BestVX = NaN(1, n_CH);

posM = NaN(n_CH,1);   %positive peack of the measured signal
negM = NaN(n_CH,1);   %negative peack of the measured signal
AMP_M = NaN(n_CH,1);  %amplitude measured signal
posS = NaN(n_CH,1);   %positive peack of the simulated signal
negS = NaN(n_CH,1);   %negative peack of the simulated signal
AMP_S = NaN(n_CH,1);  %amplitude simulated signal
diff_AMP = NaN(n_CH,1);

wind = start:ending;
for c = 1:n_CH
    sig = A{1}(chan_idx(c),:);
    yM = sig(wind);    % measured
    yM = yM(:);
    Lm = length(yM);
      
    % Steepest downslope (minimum of the derivative) 
    [~, iMd] = min(diff(yM));  
    idxM = iMd + 1; 
    
    for vx = 1:nV
        E1 = (MODEL.VENTR.VENTRICLES(vx,:)*S_nodes);    %simulated
        yS = E1(:);
    
        [~, iSd] = min(diff(yS));  
        idxS = iSd + 1; 
    
        % Temporal Shift 
        lag = idxM - idxS;          
        if lag > 0
            yS_shift = [zeros(lag,1); yS(1:end-lag)];  %if >0 the simulated is brought forward
        elseif lag < 0
            k = abs(lag);
            yS_shift = [yS(1+k:end); zeros(k,1)];      %if <0 the simulated is brought backward
        else
            yS_shift = yS;
        end
        
        % Padding to have the same length of the signal
        if length(yS_shift) < Lm
            yS_shift = [yS_shift; zeros(Lm - length(yS_shift),1)];
        else
            yS_shift = yS_shift(1:Lm);
        end
    
        % correlation computation
        C = corrcoef(yM(:), yS_shift(:));
        r = C(1,2);
        Corr(vx,c) = r;   %table
        if r > BestR(c)
            BestR(c)  = r;
            BestVX(c) = vx;
            bestYS    = yS_shift; 
        end
        
    end
    posM(c) = max(yM);                  %positive peack of the measured signal
    negM(c) = min(yM);                  %negative peack of the measured signal
    AMP_M(c) = posM(c) - negM(c);       %amplitude measured signal
    posS(c) = max(bestYS);              %positive peack of the simulated signal
    negS(c) = min(bestYS);              %negative peack of the simulated signal
    AMP_S(c) = posS(c) - negS(c);       %amplitude simulated signal
    diff_AMP(c) = AMP_S(c) - AMP_M(c);

end


%% Parameter sweep p1, p3, p5 for the best vertex of each channel

p1 = 0.1:0.2:2;          % depolarization slope
p2 = 0;                  % depolarization baseline
p3 = -0.050:0.010:0.010; % repolarization slope
p4 = 0;                  % repolarization baseline
p5 = 0:0.20:1;           % plateau slope (ST section)

% create values set and get number of combinations
valueSets = {p1, p2, p3, p4, p5};
totalComb = prod(cellfun(@numel, valueSets));

BestP5 = NaN(n_CH,1);
BestR_tuned = -inf(n_CH,1);

RMSE   = NaN(n_CH,1);
MAE    = NaN(n_CH,1);
BIAS   = NaN(n_CH,1);
dAUC   = NaN(n_CH,1);
RMSE_QRS = NaN(n_CH,1);

for c = 1:n_CH
    vx_best = BestVX(c);

    sig = A{1}(chan_idx(c),:);
    yM  = sig(wind);
    yM  = yM(:);
    Lm  = length(yM);

    [~, iMd] = min(diff(yM));
    idxM     = iMd + 1;

    best_r  = -Inf;
    best_p1 = NaN;
    best_p3 = NaN;
    best_p5 = NaN;
    yS_best_shift = [];  %we store the best alligned signal 

    % Evaluation of different parameters
    for i = 1:totalComb
        % get values set for this iteration
        p = getParams(valueSets, i, totalComb, sprintf('ChanName %s', string(chan_names(c))));

        S_nodes_sweep = gets_v2(T, dep, rep, p, mode);

        % EGM for the best vertex
        yS = (MODEL.VENTR.VENTRICLES(vx_best,:)*S_nodes_sweep).';
    
        % steepest downslope allignment
        [~, iSd] = min(diff(yS));
        idxS     = iSd + 1;
        lag      = idxM - idxS;

        if lag > 0
            yS_shift = [zeros(lag,1); yS(1:end-lag)];
        elseif lag < 0
            k = abs(lag);
            yS_shift = [yS(1+k:end); zeros(k,1)];
        else
            yS_shift = yS;
        end

        % padding
        if length(yS_shift) < Lm
            yS_shift = [yS_shift; zeros(Lm - length(yS_shift),1)];
        else
            yS_shift = yS_shift(1:Lm);
        end

        % new evaluation for correlation
        C = corrcoef(yM, yS_shift);
        r = C(1,2);

        % update best values
        if  r > best_r
            best_r  = r;
            best_p1 = p(1);
            best_p3 = p(3);
            best_p5 = p(5);
            yS_best_shift = yS_shift;   
        end
            
    end

    BestR_tuned(c) = best_r;
    BestP1(c) = best_p1;
    BestP3(c) = best_p3;
    BestP5(c) = best_p5;

    % metrics for a common tract 
    Lc = min(length(yM), length(yS_best_shift)); 
    a = yM(1:Lc); 
    b = yS_best_shift(1:Lc);

    RMSE(c) = sqrt(mean((a-b).^2));
    MAE(c)  = mean(abs(a-b));
    BIAS(c) = mean(b-a);
    dAUC(c) = trapz(a) - trapz(b);

    W = 60;
    [~,iMd_loc] = min(diff(a)); iM = iMd_loc+1;
    i1 = max(1,iM-W); i2 = min(Lc,iM+W);
    RMSE_QRS(c) = sqrt(mean((a(i1:i2)-b(i1:i2)).^2));

end

ParamTable = table( ...
    chan_idx(:), ...
    string(chan_names(:)), ...
    BestVX(:), ...
    [BestP1(:) BestP3(:) BestP5(:)], ...
    BestR(:), ...
    BestR_tuned(:), ...
    diff_AMP(:), ...
    RMSE(:), ...
    dAUC(:), ...
    RMSE_QRS(:), ...
    'VariableNames',{ ...
    'ChanIdx', ...
    'ChanName', ...
    'BestVertex', ...
    'BestP1, BestP3, BestP5', ...
    'BestR_old', ...
    'BestR_new', ...
    'Amplitude difference', ...
    'rmse(:)', ...
    'dAUC(:)', ...
    'rmse_QRS(:)' ...
    } );
disp(ParamTable);


