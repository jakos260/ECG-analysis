function ANN = loadANN(varargin)

% Version 1
% Patrick.scholten@medtronic.com
% March 3, 2008
% First version: write ANN to file

%% Check input parameters

    % check num of input arguments
    error(nargchk(1, 3, nargin));
    
    % Check outputs
    if nargout > 1
    	error('loadANN: This fuction has one output.');
    end
    
    % check first argument
    if ~ischar(varargin{1}),
        error(['loadANN: Filename variable is not a string: ' num2str(varargin{1})]);
    elseif (exist(varargin{1},'file')~=2),
        error(['loadANN: File does not exist: ' varargin{1}]);
    end
    annfile = varargin{1};
    
%% Read Header

    fpt = fopen(annfile,'r');
    if fpt < 0,
        error('loadANN: Cannot open %s', annfile);
    end
    
    % read header version number
    ANN.HDR.version = uint8(fread(fpt,1,'uint8'));
    % Check version number
    if (ANN.HDR.version ~= 1),
        error(['loadANN: Header Verion not implemented: ' ANN.HDR.version]);
    end
    
    % read size of header
    ANN.HDR.size = typecast(uint8(fread(fpt,2,'uint8')),'uint16');
    % read header information
    hdr_arr = uint8(fread(fpt,double(ANN.HDR.size)-3,'uint8'));
    
    % read date    
    ANN.HDR.fileCreated = datestr(typecast(hdr_arr(1:8),'double'),'dd-mmm-yyyy HH:MM:SS.FFF');      % When is this file created (20 bytes)
    ANN.HDR.numOfAnnotations = typecast(hdr_arr(9:10),'uint16');
    ANN.HDR.sigFilename = char(hdr_arr(11:42))';
    ANN.HDR.comment = char(hdr_arr(43:end))';

%% Read Data

    for i=1:int32(ANN.HDR.numOfAnnotations)
        ANN.DATA(i).timeStart = typecast(uint8(fread(fpt,8,'uint8')),'uint64');
        ANN.DATA(i).timeStop = typecast(uint8(fread(fpt,8,'uint8')),'uint64');
        ANN.DATA(i).code = typecast(uint8(fread(fpt,2,'uint8')),'uint16');
        ANN.DATA(i).dataField.size = typecast(uint8(fread(fpt,2,'uint8')),'uint16');
        ANN.DATA(i).dataField.data = fread(fpt,double(ANN.DATA(i).dataField.size),'uint8')';
        ANN.DATA(i).labelField.size = typecast(uint8(fread(fpt,2,'uint8')),'uint16');
        ANN.DATA(i).labelField.data = fread(fpt,double(ANN.DATA(i).labelField.size),'*char')';
    end

% ANN.DATA(2).dataField.data      = [9 8 7 6 5 4 3 2 1]; % uint8
% ANN.DATA(2).labelField.data     = 'labelField2';

    fclose(fpt);
end