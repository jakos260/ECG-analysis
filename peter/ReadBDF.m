% function D = Readbdf(Filename,[channels, starttime, endtime]['noui'])
% Read bdf file (filename) into D structure (fields: SelectChannels, chanlist,
% signals, startpos, endpos, startsample)
% Call with 1,2, 4 or 5 arguments.
% Starttime and endtime in seconds, starting at 0. So 0 to 0 is 1 second of
% data (i.e. the first second).
% Channels, starttime and endtime can be strings (pop-up default) or
% array/numbers ( pop-up suppressed ). Last input argument can force ('ui')
% or suppress pop-up ('noui').
function D = ReadBDF(Filename,varargin)
global bsg_options
numargflag=false; % Check wether channels, starttime or stoptime is numeric
if nargin==1 || nargin==2
    DefAns{1}='';
    DefAns{2}='';
    DefAns{3}='';
elseif nargin>=4
%         DefAns{1}=varargin{1};
%         DefAns{2}=varargin{2};
%         DefAns{3}=varargin{3};
    for i=1:3
        if ischar(varargin{i})
            DefAns{i}=varargin{i};
        elseif isempty(varargin{i})
            DefAns{i}='';
        else
            if i==1
                [~ , DefAns{i}]=index2blocks(varargin{i},0); % no preprocessing to preserve order
            else
                DefAns{i}=num2str(varargin{i});
            end
            numargflag=true;
        end
    end
end

nouiflag=(nargin==1) || strcmpi(varargin{end},'noui')|| (numargflag && ~strcmpi(varargin{end},'ui'));

if nouiflag
    inputstr=DefAns;
else
    dlg_title='Bdf import dialog.';
    prompt{1}=sprintf('Select channels (format: 1:5 14). [%s]',DefAns{1});
    prompt{2}=sprintf('Start time in seconds. (Empty is from start. Starts at 0)[%s]',DefAns{2});
    
    prompt{3}=sprintf('End time in seconds. (or +length. Empty is until end. 0 to 0 is 1s)[%s]',DefAns{3});
    
    inputstr=inputdlg(prompt,dlg_title,1,DefAns);
end
SelectChannels=eval(['[' inputstr{1} ']']);
D.SelectChannels=SelectChannels;

% XML2Struct(which('bsg_common.xml'));
% bsg_options.common.maxsize=500;
disp(sprintf('loading BDF file: %s',Filename));
file_info=dir(Filename);
fid=fopen(Filename);
header=fread(fid,256);
if header(1)~=255; return;end
ver=char(header(2:8)');
loc_pat_inf=char(header(9:88)');
loc_rec_inf=char(header(89:168)');
start_date=char(header(169:176)');
start_time=char(header(177:184)');
nhead=char(header(185:192)');
reserv=char(header(193:236)');
nrec=char(header(237:244)');
trec=char(header(245:252)');
ns=char(header(253:256)');

ns=eval(ns);
nhead=eval(nhead);
trec=eval(trec);
nrec=eval(nrec);

label=char(fread(fid,[16,ns])');
ttype=char(fread(fid,[80,ns])');
phdim=char(fread(fid,[8,ns])');
phmin=char(fread(fid,[8,ns])');
phmax=char(fread(fid,[8,ns])');
digmin=char(fread(fid,[8,ns])');
digmax=char(fread(fid,[8,ns])');
prefilt=char(fread(fid,[80,ns])');
fs=char(fread(fid,[8,ns])');
reserv2=char(fread(fid,[32,ns])');

fs=sscanf(fs','%d');
if ~all(fs==fs(1))
    bs_error('all channels should have the same sample rate');
end
fs=fs(1);
phmin=sscanf(phmin','%d');
phmax=sscanf(phmax','%d');
digmin=sscanf(digmin','%d');
digmax=sscanf(digmax','%d');
bdfblock=ns*fs;

if isempty(SelectChannels)
    SelectChannels=1:ns;
end
SelectChannelsOld=SelectChannels; % store original order of channels (or duplicates, heaven forbid)
[cblocks estr SelectChannels]=index2blocks(SelectChannels,1); % split into blocks, preproc SelectChannels
cblocks.gap(end)=ns-cblocks.orgend(end)+cblocks.orgstart(1)-1;

for i=1:length(SelectChannels)
    D.chanlist{i}=sprintf('s%02d',SelectChannels(i));
end




%SelectChannels(SelectChannels<0)=SelectChannels(SelectChannels<0)+ns+1;
%oostep1: negative numbers relative to last channel. Removed.

% bsg_bdf.LastBlock=(file_info.bytes-nhead)/(bdfblock*3);
% rload_total   =   bsg_bdf.LastBlock*fs;


if isempty(inputstr{2})
    startpos=0;
else
    startpos=sscanf(inputstr{2},'%d');
end

if isempty(inputstr{3})
    endpos=ceil((file_info.bytes-nhead)/(bdfblock*3))-1; %positon identifies whole block starting at 0
elseif inputstr{3}(1)=='+';
    endpos=startpos+sscanf(inputstr{3}(2:end),'%d')-1;
else
    endpos=sscanf(inputstr{3},'%d');
end







numblocks=endpos-startpos+1; % 0-0 is 1 block
status=fseek(fid,startpos*bdfblock*3,'cof');
if status~=0
    error('Error setting start position');
end
D.signals=zeros(fs*numblocks,length(SelectChannels));

status=fseek(fid,fs*(cblocks.orgstart(1)-1)*3,'cof'); % start loop at beginning of first selected channel
if status~=0
    error('Error setting start position for first channel');
end
for k=1:numblocks
    for ib=1:length(cblocks.orgstart)
        D.signals((k-1)*fs+1:k*fs,cblocks.deststart(ib):cblocks.destend(ib))=fread(fid,[fs,cblocks.length(ib)],'bit24');
        status=fseek(fid,fs*cblocks.gap(ib)*3,'cof');
        if status~=0
            if k~=numblocks
                warning('Unexpected end of file');
            end
            endpos=round(startpos+size(D.signals,1)/fs-1);
            break
        end
    end
end


fclose(fid);



D.signals=D.signals*0.001*phmax(1)/digmax(1);
% display('gain correction  from bdf header applied')

% D.signals=D.signals/8000; %oostep1 1: fixed factor because real factor was not allways configured
% display('fixed gain correction  applied!!!_')



%D.signals(:,257)=(bitand(D.signals(:,257),2)>0)-(bitand(D.signals(:,257),8)>0);
D.Ts=1000/fs;
D.Fs=fs;
D.startpos=startpos;
D.endpos=endpos;
D.startsample=D.startpos*fs+1; %samplenumbers start at 1

% Restore original channels order and duplicates, when required 
if ~isequal(SelectChannelsOld,SelectChannels)
    [Lia,Locb]=ismember(SelectChannelsOld,SelectChannels);
    if ~all(Lia)
        error('Panic, Channel(s) missing!');
    else
        D.signals=D.signals(:,Locb);
        D.chanlist=D.chanlist(Locb);
    end
end
    

