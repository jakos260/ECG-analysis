function writeSIG(varargin)
% status = writeSIG(sigfile,struct)
% 
% Examples
% Overwrite sigfile with content of SIG structure
% status = writeSIG('c:\test.sig',SIG) or writeSIG('c:\test.sig',SIG,'BOTH')
%
% Overwrite sigfile only with header information from SIG structure
% status = writeSIG('c:\test.sig',SIG,'HDR')
%
% Overwrite sigfile only with data information from SIG structure
% status = writeSIG('c:\test.sig',SIG,'DATA')
%
% Add content data from SIG structure to the excisting sigfile AND
% overwrite headerinfo 
% status = writeSIG('c:\test.sig',SIG,'BOTH','APPEND')
%
% Add only current data selection from SIG structure to the excisting
% sigfile AND overwrite headerinfo
% status = writeSIG('c:\test.sig',SIG,'BOTH','APPEND',[20 40])
%
% Add only current data selection from SIG structure to the excisting
% sigfile AND overwrite headerinfo. Store the selection between the current
% data at location '100'
% status = writeSIG('c:\test.sig',SIG,'BOTH','APPEND',[20 40],100)
%
% Add only current data selection from SIG structure to the excisting
% sigfile AND DO NOT overwrite headerinfo. Store the selection between the
% current data at location '100'.
% status = writeSIG('c:\test.sig',SIG,'DATA','APPEND',[20 40],100)
%
% Not valid ('HDR')
% writeSIG('c:\test.sig',SIG,'HDR','APPEND',[20 40],100)

% Version 1
% Patrick.scholten@medtronic.com
% February 27, 2007
% First version: write SIG to file
    
%     HDR.version         = uint8(1);
%     HDR.size            = typecast(uint16(424),'uint8');
%     HDR.file_created    = typecast(now,'uint8');
%     HDR.origin          = uint8(char('a'*ones(1,32)));
%     HDR.comment         = uint8(char('b'*ones(1,256)));                % Comment field (256 bytes)
%     HDR.Fs              = typecast(double(1000),'uint8');
%     HDR.startTime       = typecast(now,'uint8');
%     HDR.startTimeOffset = typecast(now,'uint8');
%     HDR.numOfCol        = uint8(2);
%     HDR.timeFormat      = uint8(0);
%     HDR.timeUnit        = uint8(1);
%     HDR.samplesFormat{1}= uint8(4);
%     HDR.gain{1}         = typecast(double(1),'uint8');
%     HDR.offset{1}       = typecast(double(0),'uint8');
%     HDR.samplesUnit{1}  = uint8(char('c'*ones(1,32)));
%     HDR.samplesFormat{2}= uint8(4);
%     HDR.gain{2}         = typecast(double(1),'uint8');
%     HDR.offset{2}       = typecast(double(0),'uint8');
%     HDR.samplesUnit{2}  = uint8(char('d'*ones(1,32)));

% TBD: APPEND data to file.
% case 14,FormatSamplesType{i} = 'charXXX'; % tbd: define own length of string

%% Check input parameters

    % check num of input arguments
    error(nargchk(1, 2, nargin));
    
    % check first argument (filename)
    if ~ischar(varargin{1}),
        error(['writeSIG: Filename variable is not a string: ' num2str(varargin{1})]);
    end
    sigfile = varargin{1};
    [pathstr, name, ext, versn]=fileparts(sigfile);
    if isempty(ext)
        error(['writeSIG: Filename variable is not a string: ' sigfile]);
    end
    
    % Check version
    VALID_FIELD_NAMES_SIG_v1 = {'HDR','DATA'};
    VALID_FIELD_NAMES_SIG_HDR_v1 = {'version','size','file_created',...
        'origin','comment','Fs','startTime','timecolpresent','numOfCol',...
        'label','delay','format','unit','gain','offset'};
    
    % Check second argument
    SIG = varargin{2};
    
    fnames_SIG = fieldnames(SIG);
    if sum(ismember(fnames_SIG,VALID_FIELD_NAMES_SIG_v1))~=length(VALID_FIELD_NAMES_SIG_v1),
        error('writeSIG: Structure SIG not valid, fieldnames are missing.');
    end
    fnames_SIG_HDR = fieldnames(SIG.HDR);
    if sum(ismember(fnames_SIG_HDR,VALID_FIELD_NAMES_SIG_HDR_v1))~=length(VALID_FIELD_NAMES_SIG_HDR_v1),
        error('writeSIG: Structure SIG.HDR not valid, fieldnames are missing.');
    end
    
    if isempty(SIG.DATA),
        display('writeSIG: SIG.DATA is empty.');
    end
    
    % Check third argument
    % option = varargin{3};
    
%% Create header

    if ~isnumeric(SIG.HDR.version),
        error('writeSIG: SIG.HDR.version is not numeric.');
    end
    totalhdr = uint8([]);
    totalhdr = [totalhdr uint8(SIG.HDR.version)];
    totalhdr = [totalhdr typecast(uint16(0),'uint8')]; % placeholder size of header

    if ~char(SIG.HDR.file_created),
        error('writeSIG: SIG.HDR.file_created is not a string.');
    end
    try
        if isempty(SIG.HDR.file_created)
            totalhdr = [totalhdr typecast(double(0),'uint8')];
        else
            totalhdr = [totalhdr typecast(datenum(SIG.HDR.file_created),'uint8')];
        end
    catch
        error('writeSIG: SIG.HDR.file_created can not be converted.');
    end
    
    if ~ischar(SIG.HDR.origin),
        error('writeSIG: SIG.HDR.origin is not a string.');
    end
    if length(SIG.HDR.origin) > 32,
        display('writeSIG: SIG.HDR.origin will be truncated.');
        SIG.HDR.origin = SIG.HDR.origin(1:32);
    end
    totalhdr = [totalhdr uint8([SIG.HDR.origin char(32*ones(1,32-length(SIG.HDR.origin)))])];
    
    if ~ischar(SIG.HDR.comment),
        error('writeSIG: SIG.HDR.comment is not a string.');
    end
    if length(SIG.HDR.origin) > 256,
        display('writeSIG: SIG.HDR.comment will be truncated.');
        SIG.HDR.comment = SIG.HDR.comment(1:256);
    end
    totalhdr = [totalhdr uint8([SIG.HDR.comment char(32*ones(1,256-length(SIG.HDR.comment)))])];

    if ~isnumeric(SIG.HDR.Fs),
        error('writeSIG: SIG.HDR.Fs is not numeric.');
    end
    totalhdr = [totalhdr typecast(SIG.HDR.Fs,'uint8')];

    if ~ischar(SIG.HDR.startTime),
        error('writeSIG: SIG.HDR.startTime is not a string.');
    end
    try
        if isempty(SIG.HDR.startTime)
            totalhdr = [totalhdr typecast(double(0),'uint8')];
        else
            totalhdr = [totalhdr typecast(datenum(SIG.HDR.startTime),'uint8')];
        end
    catch
        error('writeSIG: SIG.HDR.startTime can not be converted.');
    end

    if ~isnumeric(SIG.HDR.timecolpresent),
        error('writeSIG: SIG.HDR.timecolpresent is not numeric.');
    end
    totalhdr = [totalhdr uint8(SIG.HDR.timecolpresent)];

    if ~isnumeric(SIG.HDR.numOfCol),
        error('writeSIG: SIG.HDR.numOfCol is not numeric.');
    end
    totalhdr = [totalhdr uint8(SIG.HDR.numOfCol)];

    if ~ischar(SIG.HDR.label),
        error('writeSIG: SIG.HDR.label is not a string.');
    end
    if ~isnumeric(SIG.HDR.delay),
        error('writeSIG: SIG.HDR.delay is not numeric.');
    end
    if ~isnumeric(SIG.HDR.format),
        error('writeSIG: SIG.HDR.format is not numeric.');
    end
    if ~ischar(SIG.HDR.unit),
        error('writeSIG: SIG.HDR.unit is not a string.');
    end
    if ~isnumeric(SIG.HDR.gain),
        error('writeSIG: SIG.HDR.gain is not numeric.');
    end
    if ~isnumeric(SIG.HDR.offset),
        error('writeSIG: SIG.HDR.offset is not numeric.');
    end
    
    for i=1:int16(SIG.HDR.numOfCol),
        if ~ischar(SIG.HDR.label(i,:)),
            error('writeSIG: SIG.HDR.label is not a string.');
        end
        label = SIG.HDR.label(i,:);
        if length(label) > 32,
            display('writeSIG: SIG.HDR.label will be truncated.');
            label = label(1:32);
        end
        totalhdr = [totalhdr uint8([label char(32*ones(1,32-length(label)))])];

        totalhdr = [totalhdr typecast(SIG.HDR.delay(i),'uint8')];
        totalhdr = [totalhdr uint8(SIG.HDR.format(i))];

        if ~ischar(SIG.HDR.unit(i,:)),
            error('writeSIG: SIG.HDR.unit is not a string.');
        end
        unit = SIG.HDR.unit(i,:);
        if length(unit) > 32,
            display('writeSIG: SIG.HDR.unit will be truncated.');
            unit = unit(1:32);
        end
        totalhdr = [totalhdr uint8([unit char(32*ones(1,32-length(unit)))])];

        totalhdr = [totalhdr typecast(SIG.HDR.gain(i),'uint8')];
        totalhdr = [totalhdr typecast(SIG.HDR.offset(i),'uint8')];
    end

    % Overwrite size of header
    totalLengthOfHeader = uint16(length(totalhdr));
    totalhdr(2:3) = typecast(totalLengthOfHeader,'uint8');
    
%% Select formats
    for i = 1:int16(SIG.HDR.numOfCol),
        switch (SIG.HDR.format(i))
            case 0, FormatSamplesType{i} = '';
            case 1, FormatSamplesType{i} = 'uint8';
            case 2, FormatSamplesType{i} = 'int8';
            case 3, FormatSamplesType{i} = 'uint16';
            case 4, FormatSamplesType{i} = 'int16';
            case 5, FormatSamplesType{i} = 'uint32';
            case 6, FormatSamplesType{i} = 'int32';
            case 7, FormatSamplesType{i} = 'uint64';
            case 8, FormatSamplesType{i} = 'int64';
            case 9, FormatSamplesType{i} = 'float32';
            case 10,FormatSamplesType{i} = 'float64';
            case 11,FormatSamplesType{i} = 'uchar';
            case 12,FormatSamplesType{i} = 'schar';
            case 13,FormatSamplesType{i} = 'char16';
            otherwise
                error('loadSIG: Sample format is not supported');
        end
    end
    
    
%% Write to file    
    
    fpt = fopen(sigfile,'w');
    if fpt<0,
        error(['writeSIG: Can not open file: ' sigfile]);
    end
    
    % Write header to file
    fwrite(fpt,totalhdr,'uint8');
    
    % write data
    for i=1:int16(SIG.HDR.numOfCol),
        fwrite(fpt,SIG.DATA{i},FormatSamplesType{i});
    end
    
    % Close file
    fclose(fpt);
end