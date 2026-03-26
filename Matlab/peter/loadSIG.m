function SIG = loadSIG(varargin)	
% This function reads the SIG file and stores the information in one
% structure.
%
% Syntax
% SIG = loadSIG(filename,[n1 s n2],options)
% filename = .sig file (string)
% [n1 s n2] = rows to be read with interval 's'
% options: 'HDR'            : read header from sig file only
% options: 'DATA'           : read data from sig only
% options: 'BOTH' or empty  : read data from sig only
%
% Examples
% read the whole SIG file, header and data
% SIG = loadSIG('c:\demo.sig') or SIG = loadSIG('c:\demo.sig','BOTH');
% 
% read the whole SIG file, header only. The second inputvarianle will be
% neglected.
% SIG = loadSIG('c:\demo.sig','HDR') or SIG = loadSIG('c:\demo.sig',[],'HDR') ;
% 
% read the whole SIG file, data only
% SIG = loadSIG('c:\demo.sig','DATA');
% 
% read header and read only the first 1000 samples
% SIG = loadSIG('c:\demo.sig',1000) or ECG = loadECG('c:\demo.sig',1000,'BOTH');
%
% read header and read only the rows between 500 and 1000
% SIG = loadSIG('c:\data.ecg',[500 1000]) or ECG = loadECG('c:\demo.sig',[500 1000]) or
% SIG = loadSIG('c:\data.ecg',[500 1 1000]) or ECG = loadECG('c:\demo.sig',[500 1 1000]);
% 
% read header and read only the rows between 500 and 1000, skip 2 samples
% SIG = loadSIG('c:\data.ecg',[500 2 1000]) or ECG = loadECG('c:\demo.sig',[500 2 1000]);
%
% read header and read only the rows between 500 and 1000, skip 2 samples, start at sample 10
% SIG = loadSIG('c:\data.ecg',[500 2 1000 10]) or ECG = loadECG('c:\demo.sig',[500 2 1000 10]);
%
% read NO header and read only the first 1000 rows
% SIG = loadSIG('c:\data.ecg',1000,'DATA');

% Version 1: First version
% Patrick.scholten@medtronic.com
% February 27, 2007
%
% Version 2: Variables to double and 'datestr'
% February 29, 2007
%
% Version 3: skip read and offset read: [100 4 300] and [100 4 300 10]
% March 18, 2008

%% Check input parameters

    % check num of input arguments
    error(nargchk(1, 3, nargin));
    
    % Check outputs
    if nargout > 1
    	error('loadSIG: This fuction has one output.');
    end
    
    % check first argument
    if ~ischar(varargin{1}),
        error(['loadSIG: Filename variable is not a string: ' num2str(varargin{1})]);
    elseif (exist(varargin{1},'file')~=2),
        error(['loadSIG: File does not exist: ' varargin{1}]);
    end
    Input.sigfile = varargin{1};
        
    switch(nargin)
        case 1
            Input.start = uint32([]);
            Input.stop = uint32([]);
            Input.increment = 0;
            Input.offset = 0;
            Input.option = 2;
        case 2
            if isnumeric(varargin{2})
               if min(double(varargin{2}))<1,
                   if length(varargin{2})~=4
                       error('loadSIG: Minimum for second argument is 1.');
                   end
               end
               switch (length(varargin{2}))
                   case {1}
                       Input.start = uint32(0);
                       Input.stop = uint32(varargin{2}(1));
                       Input.increment = 0;
                       Input.offset = 0;
                   case {2}
                       if varargin{2}(2)<varargin{2}(1),
                           error('loadSIG: Second argument: second element is larger than first element.');
                       end
                       Input.start = uint32(double(varargin{2}(1)) - 1);
                       Input.stop = uint32(double(varargin{2}(2)));
                       Input.increment = 0;
                       Input.offset = 0;
                   case {3}
                       if varargin{2}(3)<varargin{2}(1),
                           error('loadSIG: Second argument: third element is larger than first element.');
                       end
                       if varargin{2}(3)<varargin{2}(2)
                           error('loadSIG: Second argument: increment element value is not correct.');
                       end
                       Input.start = uint32(varargin{2}(1) - 1);
                       Input.stop = uint32(varargin{2}(3));
                       Input.increment = uint32(varargin{2}(2));
                       Input.offset = 0;
                   case {4}
                       if varargin{2}(3)<varargin{2}(1),
                           error('loadSIG: Second argument: third element is larger than first element.');
                       end
                       if varargin{2}(3)<varargin{2}(2)
                           error('loadSIG: Second argument: increment element value is not correct.');
                       end
                       Input.start = uint32(varargin{2}(1) - 1);
                       Input.stop = uint32(varargin{2}(3));
                       Input.increment = uint32(varargin{2}(2));
                       Input.offset = uint32(varargin{2}(4));
                   otherwise
                       error('loadSIG: Length of array of second argument is not legal (1, 2 ,or 3).');
               end
               Input.option = 2;
            else
                if ~ismember(varargin{2},{'HDR','BOTH'})
                    error('loadSIG: Second argument is not a member of ''HDR'',''BOTH''.');
                end
                Input.option = find(strcmp(varargin{2},{'HDR','BOTH'}));
                Input.start = uint32([]);
                Input.stop = uint32([]);
                Input.offset = 0;
                Input.increment = 0;
            end
        case 3
            if ~isnumeric(varargin{2})
                error('loadSIG: Second argument is not a number (array).');
            end
            switch (length(varargin{2}))
                case {1}
                    Input.start = uint32(0);
                    Input.stop = uint32(varargin{2}(1));
                    Input.increment = 0;
                    Input.offset = 0;
                 case {2}
                     if varargin{2}(2)<varargin{2}(1),
                         error('loadSIG: Second argument: second element is larger than first element.');
                     end
                     Input.start = uint32(varargin{2}(1) - 1);
                     Input.stop = uint32(varargin{2}(2));
                     Input.increment = 0;
                     Input.offset = 0;
                 case {3}
                     if varargin{2}(3)<varargin{2}(1),
                         error('loadSIG: Second argument: third element is larger than first element.');
                     end
                     if varargin{2}(3)<varargin{2}(2)
                         error('loadSIG: Second argument: increment element value is not correct.');
                     end
                     Input.start = uint32(varargin{2}(1) - 1);
                     Input.stop = uint32(varargin{2}(3));
                     Input.increment = uint32(varargin{2}(2));
                     Input.offset = 0;
                   case {4}
                       if varargin{2}(3)<varargin{2}(1),
                           error('loadSIG: Second argument: third element is larger than first element.');
                       end
                       if varargin{2}(3)<varargin{2}(2)
                           error('loadSIG: Second argument: increment element value is not correct.');
                       end
                       if varargin{2}(3)<varargin{2}(4)
                           error('loadSIG: Second argument: offset element value is not correct.');
                       end
                       Input.start = uint32(varargin{2}(1) - 1);
                       Input.stop = uint32(varargin{2}(3));
                       Input.increment = uint32(varargin{2}(2));
                       Input.offset = uint32(varargin{2}(4));
                otherwise
                     error('loadSIG: Length of array of second argument is not legal (1, 2 ,or 3).');
            end
            if ~ismember(varargin{3},{'HDR','BOTH'})
                error('loadSIG: Third argument is not a member of ''HDR'',''BOTH''.');
            end
            Input.option = find(strcmp(varargin{3},{'HDR','BOTH'}));
    end
    
%% Read Header

    files = dir(Input.sigfile);
    sizeOfFile = files.bytes;
    
    fpt = fopen(Input.sigfile,'r');
    if fpt < 0,
        error('loadSIG: Cannot open %s', Input.sigfile);
    end
    
    % read header version number
    SIG.HDR.version = uint8(fread(fpt,1,'uint8'));
    
    % Check version number
    if (SIG.HDR.version ~= 1)
        error(['loadSIG: Header Verion not implemented: ' SIG.HDR.version]);
    end

    % read size of header
    SIG.HDR.size = typecast(uint8(fread(fpt,2,'uint8')),'uint16');

    % read header information
    hdr_arr = uint8(fread(fpt,double(SIG.HDR.size)-3,'uint8'));

    % Create header structure
    SIG.HDR.file_created        = datestr(typecast(hdr_arr(1:8),'double'),'dd-mmm-yyyy HH:MM:SS.FFF');      % When is this file created (20 bytes)
    SIG.HDR.origin              = char(hdr_arr(9:40))';                 % Origin of the data,e.g. study name (32 bytes)
    SIG.HDR.comment             = char(hdr_arr(41:296))';               % Comment field (256 bytes)
    SIG.HDR.Fs                  = typecast(hdr_arr(297:304),'double');  % Samplefrequency of the data (4 bytes)
    SIG.HDR.startTime           = datestr(typecast(hdr_arr(305:312),'double'),'dd-mmm-yyyy HH:MM:SS.FFF'); 
    SIG.HDR.timecolpresent      = hdr_arr(313);                         % Is first column time vector?
    SIG.HDR.numOfCol            = hdr_arr(314);                         % Number of data columns (1 bytes)
    
    if SIG.HDR.numOfCol==0
        SIG.DATA = {};
        return;
    end

    for i = 1:int16(SIG.HDR.numOfCol),
        SIG.HDR.label(i,:)  = char(hdr_arr(315+(i-1)*89:346+(i-1)*89))';            % Label (32 bytes)
        SIG.HDR.delay(i,:)  = typecast(hdr_arr(347+(i-1)*89:354+(i-1)*89),'double');% Delay
        SIG.HDR.format(i)   = hdr_arr(355+(i-1)*89);                                % Format of each data column (1 byte)
        SIG.HDR.unit(i,:)   = char(hdr_arr(356+(i-1)*89:387+(i-1)*89))';            % Unit of each data column (32 bytes)
        SIG.HDR.gain(i)     = typecast(hdr_arr(388+(i-1)*89:395+(i-1)*89),'double');% Gain of each data column (4 bytes)
        SIG.HDR.offset(i)   = typecast(hdr_arr(396+(i-1)*89:403+(i-1)*89),'double');% Offset of each data column (4 bytes)
    end

%% Select formats

    for i = 1:int16(SIG.HDR.numOfCol),
        switch (SIG.HDR.format(i))
            case 0, FormatSamplesType{i} = '';
                    ByteLengthSamples(i) = 0;
            case 1, FormatSamplesType{i} = 'uint8';
                    ByteLengthSamples(i) = 1;
            case 2, FormatSamplesType{i} = 'int8';
                    ByteLengthSamples(i) = 1;
            case 3, FormatSamplesType{i} = 'uint16'; % ushort
                    ByteLengthSamples(i) = 2;
            case 4, FormatSamplesType{i} = 'int16'; %short
                    ByteLengthSamples(i) = 2;
            case 5, FormatSamplesType{i} = 'uint32';
                    ByteLengthSamples(i) = 4;
            case 6, FormatSamplesType{i} = 'int32'; %int
                    ByteLengthSamples(i) = 4;
            case 7, FormatSamplesType{i} = 'uint64';
                    ByteLengthSamples(i) = 8;
            case 8, FormatSamplesType{i} = 'int64';
                    ByteLengthSamples(i) = 8;
            case 9, FormatSamplesType{i} = 'float32'; % single
                    ByteLengthSamples(i) = 4;
            case 10,FormatSamplesType{i} = 'float64'; % double
                    ByteLengthSamples(i) = 8;
            case 11,FormatSamplesType{i} = 'uchar'; % char
                    ByteLengthSamples(i) = 1;
            case 12,FormatSamplesType{i} = 'schar';
                    ByteLengthSamples(i) = 1;
            case 13,FormatSamplesType{i} = 'char16';
                    ByteLengthSamples(i) = 2;
            otherwise
                error('loadSIG: Sample format is not supported');
        end
    end
    
%% Calculated number of rows
    NumOfRows = (double(sizeOfFile) - double(SIG.HDR.size))/sum(double(ByteLengthSamples));
    
    % Set start and stop parameters
    if isempty(Input.start) && isempty(Input.stop)
        Input.start = double(0);
        Input.stop = double(NumOfRows);
    end
    
    if Input.stop > NumOfRows,
        error(['loadSIG: Second argument is not valid (too large): ' num2str(Input.stop)])
    end

%% Read Data

    if Input.option == 1,
        % read no data
        SIG.DATA = [];
    else % Input.option = 2
        block_offset = 0;
        for i=1:SIG.HDR.numOfCol,
            row_offset = double(Input.start)*double(ByteLengthSamples(i));
            fseek(fpt,double(Input.offset*ByteLengthSamples(i)) + double(row_offset) + double(block_offset) + double(SIG.HDR.size), 'bof');
            [SIG.DATA{i},count] = fread(fpt,[1,double(floor((double(Input.stop) - double(Input.start))/(Input.increment+1)))],['*' FormatSamplesType{i}],Input.increment*ByteLengthSamples(i));
            if count ~= double(floor((double(Input.stop) - double(Input.start))/(Input.increment+1))),
                error('loadSIG: Not all samples are read.');
            end
            block_offset = double(block_offset) + double(ByteLengthSamples(i)) * NumOfRows;
        end
    end
    fclose(fpt);
end