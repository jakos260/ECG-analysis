% Wrapper function for selectBSMBeat. 
% Input 
% -1- the bdf data filename and
% -2- start time (in seconds)
% -3- end time (in seconds)
% Additionally the sample frequency can be added. Default it is 1000 Hz. The window needs to be
% closed before the output parameters are given


function varargout = selectHagaBSMBeats(varargin)

global DATA
sampT = 1/1000;
VER =[];
ITRI=[];
elecs=[];
usemean = 0;
pp=4;
filename = varargin{1};
t0 = varargin{2};
t1 = varargin{3};
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
            case 'usemean'
                usemean = varargin{pp+1};pp=pp+2;
            otherwise
                error('unknown parameter');
        end
    end
end
if ~isempty(strfind(filename,'bdf'))
    [H,BSM] = readHagaBDF(filename,t0,t1);
else
    BSM = loadmat(filename);
end



if	sampT ~= 1/1000
    BSM = resample(BSM,1,size(BSM,2),size(BSM,2)*sampT * 1000);   
end
if ispc
    dirout= '.\beats\' ;
else
    dirout= './beats/' ;
end
if ~exist(dirout,'dir')
    mkdir(dirout)
end
filename = [strrep(filename,'.','_') '_' num2str(t0) '-' num2str(t1) 's'];

dosave = strcmp(questdlg('Do you want to save the output files'), 'Yes');

if dosave
    savemat([filename  '.mat'],BSM);
end


waitfor( selectBSMbeat(BSM,'lay',loadmat('nim65.mla'),'sampt',1/1000) );

if dosave
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
                        savemat([dirout filename '_beat' num2str(i) '.selecg'], BSMcorrect);
                    end
                    k = k+1;
                end
            end    
        end
    else
        BSMOUT=[];
        TIMES=[];
        k=1;
        for i=1:length(DATA.SELBEATS)
            if DATA.SELBEATS(i)
                BSMcorrect = DATA.BSMOUT(:,DATA.BEATS(i):DATA.BEATS(i+1));
                BSMOUT{k} = BSMcorrect;
                TIMES{k}.t0 = DATA.BEATS(i);
                TIMES{k}.t1 = DATA.BEATS(i+1);
                if ~isempty(filename)
                    savemat([dirout filename '_beat' num2str(i) '.selecg'], BSMcorrect);
                    saveasci([dirout filename '_beat' num2str(i) '.remove'], DATA.remove);
                end
                k = k+1;
            end
        end   

    end
end

varargout{1} = BSMOUT;
if nargout > 1
    varargout{2}=DATA.remove;
end

if nargout > 2
    varargout{3}= TIMES;
end

