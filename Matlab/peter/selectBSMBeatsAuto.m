% Wrapper function for selectBSMBeat. Input shoudl be the BSM data file and
% the used layout (the file shoudl have been read) Additionally the sample
% frequency can be added. Default it is 1000 Hz. The window needs to be
% closed before the output parameters are given
% 1 parameter is the baseline corrected BSM (Spline)
% optional 2 are the used baseline points
% optional 3 are teh signals identified as being bad


function varargout = selectBSMBeatsAuto(varargin)

global DATA
sampT = 1/1000;
VER =[];
ITRI=[];
elecs=[];
LAY = [1 9 0;[ ones(9,1) (1:9)' (1:9)']];
filename=[ ];
usemean = 0;
zoomfull = 1;
pp=2;
BSM = varargin{1};
while pp<=length(varargin)
    if ischar(varargin{pp})
        key=lower(varargin{pp});
        switch key
            case 'sampt'
                sampT = varargin{pp+1};pp=pp+2;
            case 'ver'    
                VER = varargin{pp+1};pp=pp+2;
            case 'itri'
                ITRI = varargin{pp+1};pp=pp+2;
            case 'elecs'
                elecs = varargin{pp+1};pp=pp+2;
            case 'filename'
                filename = varargin{pp+1};pp=pp+2;
            case 'usemean'
                usemean = varargin{pp+1};pp=pp+2;
            case 'zoomfull'
                zoomfull=varargin{pp+1};pp=pp+2;
            case 'lay'
                LAY = varargin{pp+1};pp=pp+2;
            otherwise
                error('unknown parameter');
        end
    end
end

% if std(std(BSM)) > 0.1
%     BSM = BSM/1000;
% end
if sampT ~= 1/1000
    BSM = resample(BSM,1,size(BSM,2),size(BSM,2)*sampT*1000);   
end

waitfor( selectBSMbeatAuto(BSM,'lay',LAY,'sampt',1/1000,'zoomfull',zoomfull ) );


%interpolate the badsigs with the laplacian and save the output
if ~isempty(VER) && ~isempty(ITRI)
%     fn=filename;
%     index=strfind(fn,'.');
%     if ~isempty(index)
%         fn = fn(1:index(end)-1);
%     end
    nodes = 1:length(VER);
    if  ~isempty(elecs)
        nodes=elecs;
        nodes(elecs(DATA.remove))=[];
    else
        nodes(nodes > size(DATA.BSMOUT,1)) =[];
        nodes(DATA.remove)=[]; % assume the first n vertices represent the elctrodes        
    end
    
    T = intripol(VER,ITRI,nodes);
    if usemean && sum(DATA.SELBEATS) == 0
        BSMOUT{1} = T * DATA.meanBeat(nodes,:);
    else
        BSMOUT=[];
        TIMES=[];
        k=1;
        for i=1:length(DATA.SELBEATS)
            if DATA.SELBEATS(i)
                BSMcorrect = T * DATA.BSMOUT(nodes,DATA.BEATS(i):DATA.BEATS(i+1));
                BSMOUT{k} = BSMcorrect;
                TIMES{k}.t0 = DATA.BEATS(i);
                TIMES{k}.t1 = DATA.BEATS(i+1);
                if ~isempty(filename)
                    savemat([filename '_beat' num2str(i) '.selecg'], BSMcorrect);
                end
                k = k+1;
            end
        end    
    end
else
    BSMOUT=[];
    TIMES=[];
    k=1;figure(1091);clf;
    for i=1:length(DATA.SELBEATS)
        if DATA.SELBEATS(i)
            BSMcorrect = DATA.BSMOUT(:,DATA.BEATS(i):DATA.BEATS(i+1));
            figure(1091);sigplot(BSMcorrect,'lay',LAY);
            BSMOUT{k} = BSMcorrect;
            TIMES{k}.t0 = DATA.BEATS(i);
            TIMES{k}.t1 = DATA.BEATS(i+1);
            if ~isempty(filename)
                fn = [filename '_beat' num2str(i) '_t0_' num2str(TIMES{k}.t0)];
                savemat([fn '.selecg'], BSMcorrect);
                saveasci([fn '.remove'], DATA.remove);
                
%                 GEOM.SPECS = prepareECGAuto(BSMcorrect,LAY,'documsum',1,'filename',fn,'dosave',2);
            end
            k = k+1;
        end
    end  
    disp('median')
    
end
if nargout >= 1
    varargout{1} = BSMOUT;
end
if nargout > 1
    varargout{2}=DATA.remove;
end

if nargout > 2
    varargout{3}= TIMES;
end

