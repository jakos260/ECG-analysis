function writeANN(varargin)

% Version 1
% Patrick.scholten@medtronic.com
% March 3, 2008
% First version: write ANN to file

%% Check input parameters

    % check num of input arguments
    error(nargchk(1, 2, nargin));
    
    % check first argument (filename)
    if ~ischar(varargin{1}),
        error(['writeANN: Filename variable is not a string: ' num2str(varargin{1})]);
    end
    annfile = varargin{1};
    [pathstr, name, ext, versn] = fileparts(annfile);
    if isempty(ext)
        error(['writeANN: Filename variable is not a string: ' sigfile]);
    end

%% Check fieldnames of structure
    % Check version
    VALID_FIELD_NAMES_ANN_v1 = {'HDR','DATA'};
    VALID_FIELD_NAMES_ANN_HDR_v1 = {'version','sizeOfHeader',...
            'fileCreated','sigFilename','numOfAnnotations','comment'};    
    VALID_FIELD_NAMES_ANN_DATA_v1 = {'timeStart','timeStop',...
            'code','dataField','labelField'};
    VALID_FIELD_NAMES_ANN_DATA_sub_v1 = {'size','data'};
        
    % Check second argument
    ANN = varargin{2};
    
    fnames_ANN = fieldnames(ANN);
    if sum(ismember(fnames_ANN,VALID_FIELD_NAMES_ANN_v1))~=length(VALID_FIELD_NAMES_ANN_v1),
        error('writeANN: Structure ANN not valid, fieldnames are missing.');
    end
    fnames_ANN_HDR = fieldnames(ANN.HDR);
    if sum(ismember(fnames_ANN_HDR,VALID_FIELD_NAMES_ANN_HDR_v1))~=length(VALID_FIELD_NAMES_ANN_HDR_v1),
        error('writeANN: Structure ANN.HDR not valid, fieldnames are missing.');
    end
    fnames_ANN_DATA = fieldnames(ANN.DATA);
    if sum(ismember(fnames_ANN_DATA,VALID_FIELD_NAMES_ANN_DATA_v1))~=length(VALID_FIELD_NAMES_ANN_DATA_v1),
        error('writeANN: Structure ANN.DATA not valid, fieldnames are missing.');
    end

    % Check content of DATA
    if isempty(ANN.DATA),
        display('writeANN: ANN.DATA is empty.');
    end
    if length(ANN.DATA) ~= ANN.HDR.numOfAnnotations,
        error('writeANN: length of ANN.DATA is not equal to ANN.HDR.numOfAnnotations.');
    end
    
%% Create header of ANN file

    if ~isnumeric(ANN.HDR.version),
        error('writeANN: ANN.HDR.version is not numeric.');
    end
    totalhdr = uint8(ANN.HDR.version);
    totalhdr = [totalhdr typecast(uint16(0),'uint8')]; % placeholder size of header

    if ~ischar(ANN.HDR.fileCreated),
        error('writeANN: ANN.HDR.fileCreated is not a string.');
    end
    try
        if isempty(ANN.HDR.fileCreated)
            totalhdr = [totalhdr typecast(double(0),'uint8')];
        else
            totalhdr = [totalhdr typecast(datenum(ANN.HDR.fileCreated),'uint8')];
        end
    catch
        error('writeANN: ANN.HDR.fileCreated can not be converted.');
    end
    
    if ~isnumeric(ANN.HDR.numOfAnnotations),
        error('writeANN: ANN.HDR.numOfAnnotations is not numeric.');
    end
    totalhdr = [totalhdr typecast(uint16(ANN.HDR.numOfAnnotations),'uint8')];
    if ~ischar(ANN.HDR.sigFilename),
        error('writeANN: ANN.HDR.sigFilename is not a string.');
    end

    if length(ANN.HDR.sigFilename)>32,
        display('writeANN: ANN.HDR.sigFilename will be truncated.');
        totalhdr = [totalhdr uint8(ANN.HDR.sigFilename(length(ANN.HDR.sigFilename)-32:end))];
    end
    totalhdr = [totalhdr uint8([ANN.HDR.sigFilename char(32*ones(1,32-length(ANN.HDR.sigFilename)))])];
    
    if ~ischar(ANN.HDR.comment),
        error('writeANN: ANN.HDR.comment is not a string.');
    end
    totalhdr = [totalhdr uint8(ANN.HDR.comment)];
    
    % Overwrite size of header
    totalLengthOfHeader = uint16(length(totalhdr));
    totalhdr(2:3) = typecast(totalLengthOfHeader,'uint8');
 
%% Create data part of ANN file

    totaldata = uint8([]);
    for i=1:length(ANN.DATA)
        if ~isnumeric(ANN.DATA(i).timeStart),
            error(['writeANN: ANN.DATA.timeStart is not numeric, number: ' num2str(i)]);
        end
        totaldata = [totaldata typecast(uint64(ANN.DATA(i).timeStart),'uint8')];

        if ~isnumeric(ANN.DATA(i).timeStop),
            error('writeANN: ANN.DATA.timeStop is not numeric.');
        end
        totaldata = [totaldata typecast(uint64(ANN.DATA(i).timeStop),'uint8')];

        if ~isnumeric(ANN.DATA(i).code),
            error('writeANN: ANN.DATA.code is not numeric.');
        end
        totaldata = [totaldata typecast(uint16(ANN.DATA(i).code),'uint8')];

        if ~isnumeric(ANN.DATA(i).dataField.data),
            error('writeANN: ANN.DATA.dataField.data is not numeric.');
        end
        totaldata = [totaldata typecast(uint16(length(ANN.DATA(i).dataField.data)),'uint8')];
        totaldata = [totaldata ANN.DATA(i).dataField.data];

        if ~ischar(ANN.DATA(i).labelField.data),
            error('writeANN: ANN.DATA.labelField.data is not a string.');
        end
        totaldata = [totaldata typecast(uint16(length(ANN.DATA(i).labelField.data)),'uint8')];
        totaldata = [totaldata uint8(ANN.DATA(i).labelField.data)];
    end

%% Write to file    
    
    fpt = fopen(annfile,'w');
    if fpt<0,
        error(['writeANN: Can not open file: ' annfile]);
    end
    % Write header to file
    fwrite(fpt,totalhdr,'uint8');
    % write data
    fwrite(fpt,totaldata,'uint8');
    % Close file
    fclose(fpt);

end