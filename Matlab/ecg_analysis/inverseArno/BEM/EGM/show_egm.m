%% Script to show the EGM data from CARTO and select the latencys, based on the EGM's.
%
% Version 1.0: 10-aug-2015

clear all
close all
clc

% Options for visualization:
fig_opt                 = 2; % [1] show all 18 channels, [2] show only M1, Lead II and RMS curve
check_previous_result   = 0;
min_amp                 = 2;
max_dis                 = 10;
markersize_1            = 6;
markersize_2            = 8;

% If directory name is not defined:
subject = 'NICE001'; fpath = ['/Users/arnojanssen/Documents/STW/PVCs/' subject];

% load CARTO data:
load([fpath '/CATHDATA/CARTO.mat'], 'EGM', 'CARTO');

% load mesh of ventricles & translated/rotated carto points:
load([fpath '/CATHDATA/cartothorax_scirun.mat'], 'ventricles', 'newcartopoints', 'proj_cartopoints');

% make vector with EGM numbers:
egm_nr = 1:size(EGM,2);

% Find 'bad' points:
delsel = find(CARTO.points(:,8) == -1);

% file latencies:
lat     = CARTO.map_anno(:,8);
ed_lat  = CARTO.map_anno(:,7);

% select map index:
MAP = CARTO.points(:,2);

% Remove points from EGM & latencies:
EGM(delsel)     = [];
lat(delsel)     = [];
ed_lat(delsel)  = [];
egm_nr(delsel)  = [];
MAP(delsel)     = [];

% Signal names from SQL: POINT_TABLE & CONFIG_CHANNEL_TABLE:
ch_name = {'UNI-MAP-1'; 'UNI-MAP-2'; 'UNI-MAP-3'; 'UNI-MAP-4'; ...
    'V1'; 'V2'; 'V3'; 'V4'; 'V5'; 'V6'; ...
    'BI-M1-M2'; 'BI-M3-M4'; ...
    'I'; 'II'; 'III'; 'aVL'; 'aVR'; 'aVF'; ...
    'RMS'};

% window of interest according to SQL database:
woifrom = -150;
woito   = 300;

% referencetime according to SQL database:
reftime = 2000 + 1; % The extra [1] is to correct for teh selection of the data from the raw .mpt files. One sample is added at the beginning...

% x and y axes limitations:
xmin    = reftime + woifrom;
xmax    = reftime + woito;
tvec    = xmin:xmax;

% Add rms of v1 till v6, aVL, aVR & aVF
for z = 1:size(EGM,2),
    rmssigs = [];
    for selsigs = [5:10 16:18], % V1 till V6, aVL, aVR & aVF
        rmssigs = [rmssigs; EGM(z).sigs(selsigs,:)];
    end
    
    EGM(z).sigs = [EGM(z).sigs; rms(zeromean(rmssigs),2)]; 
end

% Check if earlier analyses exist:
if exist([fpath '/CATHDATA/EGM_selected.mat'],'file') == 2 && check_previous_result == 1,
    load([fpath '/CATHDATA/EGM_selected.mat']);
    pnt         = last_pnt+1;
else
    pnt         = 1;
    remove_egm  = [];
    EGM_markers = zeros(size(EGM,2),3);
end

% start qtriplot
qtriplot
qtriplot('delete *');
qtriplot('bgdcolor white');
qtriplot(ventricles.node,ventricles.face)
qtriplot('trans 0.2')


%% show overview for selected CARTO point [pnt]:
figure(fig_opt)
stopflag    = false;
fprintf('\nPress [space] to proceed with next EGM, [r] to mark EGM as bad point, [return] to jump 10 forward, [esc] to exit') ;

if fig_opt == 1, tot_subfigures = 1:18;         xsub = 3; ysub = 6; txtmarker = 18; end
if fig_opt == 2, tot_subfigures = [1 14 19];    xsub = 1; ysub = 3; txtmarker = 19; end

while ~stopflag && pnt < size(EGM,2),
    cartox  = reftime + ed_lat(pnt);
    carto8  = reftime + lat(pnt);
    
    % show location in qtriplot:
    qtriplot(newcartopoints.node(pnt,:),[]);
    qtriplot('color white')
    qtriplot(proj_cartopoints.node(pnt,:),[]);
    if MAP(pnt) == 1, mp = 'red'; elseif MAP(pnt) == 2, mp = 'blue'; elseif MAP(pnt) == 3, mp = 'green'; elseif MAP(pnt) == 4, mp = 'black'; end
    qtriplot(['color ' mp])
    
    handle_subplot = [];
    
    for k = tot_subfigures,
        
        if fig_opt == 1,
            if k > 4 && k < 11, n = k + 2; elseif k == 11 || k == 12, n = k - 6; else n = k; end
        else
            if k == 19, n = 3; elseif k == 14, n = 2; else n = 1; end
        end
        
        % calculate amplitude:
        sfig    = abs(max(EGM(pnt).sigs(k,tvec)) - min(EGM(pnt).sigs(k,tvec)));
        amp     = (max(EGM(pnt).sigs(k,xmin:reftime)) - min(EGM(pnt).sigs(k,xmin:reftime)))/1000;
        if amp < min_amp && n == 1, amp_color = [1,0,0]; sigcolor = 'r'; else amp_color = [0,0,0]; sigcolor = 'b'; end
        
        % check distance:
        if proj_cartopoints.field(:,pnt) > max_dis && n == 1, dis_color = [1,0,0]; sigcolor = 'r'; else dis_color = [0,0,0]; end
        
        handle_subplot = [handle_subplot subplot(xsub,ysub,n)];
        plot(tvec,EGM(pnt).sigs(k,tvec), sigcolor)
        hold on
        
        % show time derivative of signal:
        diffsig = [diff(EGM(pnt).sigs(k,:)) 0] + mean(EGM(pnt).sigs(k,tvec));
        plot(tvec, diffsig(tvec), 'r'); hold on;
        
        plot(reftime, EGM(pnt).sigs(k,reftime),'or','MarkerSize', markersize_1); hold on
        plot(reftime, EGM(pnt).sigs(k,reftime),'xr','MarkerSize', markersize_1); hold on
        
        % find minimum value in derivative:
        tsel    = find(diffsig(tvec) == min(diffsig(tvec)));
        seldif  = tvec(tsel);
                
        % show markers:
        plot(carto8,EGM(pnt).sigs(k,carto8),'om','MarkerSize', markersize_2); hold on
        plot(carto8,EGM(pnt).sigs(k,carto8),'xm','MarkerSize', markersize_2); hold on
        
        plot(seldif,EGM(pnt).sigs(k,seldif),'ok','MarkerSize', markersize_1); hold on
        plot(seldif,EGM(pnt).sigs(k,seldif),'xk','MarkerSize', markersize_1); hold on
                
        % define axes and figure labels:
        xlim([xmin xmax])
        ylim([min(EGM(pnt).sigs(k,tvec))-10 max(EGM(pnt).sigs(k,tvec))+10])
        title([ch_name{k,:}])
        xlabel('ms')
        ylabel('microV')
        
        % visualize amplitude of EGM:
        text(xmin+10, min(EGM(pnt).sigs(k,tvec)) + 0.95*sfig, ['amp = ' num2str(amp) ' mV'], 'Color', amp_color)
        
        if k == txtmarker,
            text(xmin+10, min(EGM(pnt).sigs(k,tvec)) + 0.85*sfig, ['dis = ' num2str(proj_cartopoints.field(:,pnt)) ' mm'], 'Color', dis_color);
            text(xmin+10, min(EGM(pnt).sigs(k,tvec)) + 0.75*sfig, ['lat = ' num2str(ed_lat(pnt)) ' ms']);
            text(xmin+10, min(EGM(pnt).sigs(k,tvec)) + 0.65*sfig, ['pnt = ' num2str(egm_nr(pnt)) ]);
        end
        
    end
    
    if fig_opt == 2,
        % get markers:
        tmarkers                        = zeros(2,3);
        [tmarkers(1,1), y, buttonpress] = ginput(1); %
        tmarkers(2,1)                   = gca;
        
        if buttonpress == 2, % middle mouse press
            stopflag = true;
            close all;
            clc;
            
        elseif buttonpress == 3, % right mouse press
            remove_egm      = [remove_egm pnt];
            pnt             = pnt + 1;
            clf;
            
        elseif buttonpress == 1, % left mouse press
            
            plot(tmarkers(1,1),y,'ok','MarkerSize', markersize_2, 'MarkerFaceColor','g'); hold on
            
            for d = 2:3,
                [tmarkers(1,d), y] = ginput(1);
                plot(tmarkers(1,d),y,'ok','MarkerSize', markersize_2, 'MarkerFaceColor','g'); hold on
                
                tmarkers(2,d)       = gca;
            end
            
            [dumy,sort_index_1]     = sort(tmarkers(2,:));
            t_sort_egm              = tmarkers(1,sort_index_1);
            [egm_rms,sort_index_2]  = sort(t_sort_egm(1,2:3));
            EGM_markers(pnt,:)      = round([t_sort_egm(1) egm_rms]);
            
            pnt = pnt + 1;
            clf;
            
        end
        
    elseif fig_opt == 1,
        
        [xdum, ydum, buttonpress] = ginput(1);
        
        if buttonpress == 2, stopflag = true; close all; clc;
        elseif buttonpress == 3, pnt = pnt + 10; clf;
        elseif buttonpress == 1, pnt = pnt + 1; clf;
        end
        
    end
    
    % remove points:
    qtriplot('delete')
    qtriplot('delete')
end

if fig_opt == 2,
    % show & save results:
    last_pnt = pnt;
    disp(remove_egm)
    disp(EGM_markers(1:pnt,:))
    save([fpath '/CATHDATA/EGM_selected.mat'], 'EGM', 'EGM_markers', 'remove_egm', 'last_pnt');
end