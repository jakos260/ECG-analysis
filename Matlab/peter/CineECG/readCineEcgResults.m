function varargout = readCineEcgResults(varargin)

analyseBeats = 0;

if length(varargin) < 1
    error('This routine needs at least one parameters');
else
    fn = varargin{1};
    pp=2;
    while pp<=nargin
        if ischar(varargin{pp})
            key=lower(varargin{pp});
            switch key
                case 'ecgdir'
                    ECGDIR = varargin{pp+1}; pp=pp+1;
                case 'domedian'
                    analyseBeats = 1; pp=pp+1;
                case 'doall'
                    analyseBeats = 2; pp=pp+1;
                otherwise
                    error('unknown parameter');
            end
        end
    end
end


[basedir,name,~] = fileparts(fn);
if ~exist('ECGDIR')
    ECGDIR = fullfile(basedir,'ECG_DATA');
end

XML=xml2struct(fn);
if isfield(XML,'ALVALE')
    ROOTDATA = XML.ALVALE{2}.ALVALEDATA;
    if isfield(XML.ALVALE{2}.ALVALEDATA.ECG_cases,'ECG_CASE')
        DATA = XML.ALVALE{2}.ALVALEDATA.ECG_cases.ECG_CASE.ECGfiles.ECGFiledata;
    elseif isfield(XML.ALVALE{2}.ALVALEDATA.ECG_cases,'CIPS_ECG_CASE')
        DATA = XML.ALVALE{2}.ALVALEDATA.ECG_cases.CIPS_ECG_CASE.ECGfiles.ECGFiledata;
    else
        error('not a proper case file')
    end
else
    ROOTDATA = XML.CINEECG{2}.CINEECGDATA;
    if isfield(XML.CINEECG{2}.CINEECGDATA.ECG_cases,'ECG_CASE')
        DATA = XML.CINEECG{2}.CINEECGDATA.ECG_cases.ECG_CASE.ECGfiles.ECGFiledata;
    else
        error('not a proper case file')
    end
end

if length(ROOTDATA.ECG_cases.ECG_CASE) >= 1    
    if all(size(DATA) == 1)
        D=DATA;
        clear DATA
        DATA{1} = D;
    end
    if analyseBeats==1 ||  analyseBeats==2
        DATACASES.MEDIANDATA = readMedianData(DATA,ECGDIR);
    end
    if analyseBeats==0 ||  analyseBeats==2
        DATACASES.BEATDATA = readCaseData(DATA,ECGDIR);
    end
    varargout{1} = DATACASES;
    if nargout == 2
        varargout{2} = readModelData(ROOTDATA);   
    end
end
end
%%
function val= removeText(str)
    if isfield( str,'Text')
        val = str2double(strsplit(str.Text,','));
    else
        val = str2double(struct2array(str));
    end
end
%%
function VER=read3DVertices(str)
    V=removeText(str);
    if mod(length(V),3)==0 
        VER = reshape(V,3,length(V)/3)';
    else
        VER=[];
    end
end
%%
function DATA = readMedianData(CASEDATA, ECGDIR)

iEcgv=0;
iEcga=0;
for i=1:length(CASEDATA)    
    disp(['reading median ECG file: ' CASEDATA{i}.ECGfilename.Text])   
    
    if exist(fullfile(ECGDIR,[CASEDATA{i}.ECGfilename.Text '.ecg.medianecg']),'file')
        ECG = loadmat(fullfile(ECGDIR,[CASEDATA{i}.ECGfilename.Text '.ecg.medianecg']));
    elseif exist(fullfile(ECGDIR,[CASEDATA{i}.ECGfilename.Text '.medianecg']),'file')
        ECG = loadmat(fullfile(ECGDIR,[CASEDATA{i}.ECGfilename.Text '.medianecg']));
    % elseif isfield(CASEDATA{i},'mediandata')
    %     M=mkzip;
    %     M.unzip(base64decode(CASEDATA{i}.mediandata.Text))
    %     ECG = (base64decode(CASEDATA{i}.mediandata.Text));
    end
    zerotims = [];
    if isfield(CASEDATA{i},'medianbaselinePoints')
        zerotims = removeText( CASEDATA{i}.medianbaselinePoints);  
    end
   
    
    if isfield(CASEDATA{i},'medianvresults')    
        if isfield(CASEDATA{i}.medianvresults,'MEANTSIDIAGNOSTIC')
            BEATS = CASEDATA{i}.medianvresults.MEANTSIDIAGNOSTIC;
        else
            BEATS = CASEDATA{i}.medianvresults.CIPS_MEANTSIRESULT;
        end        
        if isscalar(BEATS)
            if isfield(BEATS,'beatCineECG') && length(BEATS.beatCineECG.Text) > 1
                iEcgv=iEcgv+1;
                DATAV{iEcgv}.ECGfilename = CASEDATA{i}.ECGfilename.Text;   
                if exist('ECG')
                    BASEL=spline(zerotims,ECG(:,zerotims),1:size(ECG,2));
                    DATAV{iEcgv}.ECG = ECG - BASEL;
                    DATAV{iEcgv}.beats{1} = readResult(BEATS,DATAV{iEcgv}.ECG );
                else
                    DATAV{iEcgv}.beats{1} = readResult(BEATS,[] );
                end
                
            end
        else
            for j= 1:length(BEATS)
                if isfield(BEATS{j},'medianecg') && ~isempty(BEATS{j}.medianbaselinePoints.Text)
                    iEcgv = iEcgv +1;
                    if exist(fullfile(ECGDIR,[CASEDATA{i}.ECGfilename.Text '_' num2str(j-1) '.vmedianecg']),'file')
                        DATAV{iEcgv}.ECG = loadmat(fullfile(ECGDIR,[CASEDATA{i}.ECGfilename.Text '_' num2str(j-1) '.vmedianecg']));
                    % else
                    %     ECG = base64decode(BEATS{j}.medianecg.Text);
                    %     DATAV{iEcgv}.ECG = readECGDATA(BEATS{j}.medianecg.Text);
                    end
                    % if strcmp(CASEDATA{i}.ECGfilename.Text, '19356_2021.json')
                    %     stop=1;
                    % end
                    DATAV{iEcgv}.ECGfilename = CASEDATA{i}.ECGfilename.Text; 
                    zerotims = removeText(BEATS{j}.medianbaselinePoints);
                    zerotims(zerotims <0)=[];
                    BASEL=spline(zerotims,DATAV{iEcgv}.ECG(:,zerotims),1:size(DATAV{iEcgv}.ECG,2));
                    DATAV{iEcgv}.ECG = DATAV{iEcgv}.ECG - BASEL;                    
                    DATAV{iEcgv}.beats = readResult(BEATS{j},DATAV{iEcgv}.ECG ); 
                end
            end
        end
    end
    if isfield(CASEDATA{i},'medianaresults')
        if isfield(CASEDATA{i}.medianaresults,'MEANTSIDIAGNOSTIC')
            BEATS = CASEDATA{i}.medianaresults.MEANTSIDIAGNOSTIC;
        else
            BEATS = CASEDATA{i}.medianaresults.CIPS_MEANTSIRESULT;
        end

        if isscalar(BEATS)
            if isfield(BEATS,'beatCineECG') && length(BEATS.beatCineECG.Text) > 1
                iEcga=iEcga+1;   
                if exist('ECG')
                    if length(zerotims) >=2
                        BASEL=spline(zerotims,ECG(:,zerotims),1:size(ECG,2));
                    else
                        BASEL=zeros(size(ECG));
                    end
                    DATAA{iEcga}.ECG = ECG - BASEL;       
                    DATAA{iEcga}.beats{1} = readResult(BEATS,DATAA{iEcga}.ECG );
                    if DATAA{iEcga}.beats{1}.onsetP < zerotims(1) - 40
                        pzerotims = sort([DATAA{iEcga}.beats{1}.onsetP zerotims]);
                        BASEL=spline(pzerotims,ECG(:,pzerotims),1:size(ECG,2));
                        DATAA{iEcga}.ECG = ECG - BASEL;                
                    end
                else
                    DATAA{iEcga}.beats{1} = readResult(BEATS,[] );
                end
                DATAA{iEcga}.ECGfilename = CASEDATA{i}.ECGfilename.Text;    
                
                
            end
        else
            for j= 1:length(BEATS)
                if isfield(BEATS{j},'medianecg') && ~isempty(BEATS{j}.medianbaselinePoints.Text)
                    iEcga = iEcga + 1;
                    if exist(fullfile(ECGDIR,[CASEDATA{i}.ECGfilename.Text '_' num2str(j-1) '.amedianecg']),'file')
                        DATAA{iEcga}.ECG = loadmat(fullfile(ECGDIR,[CASEDATA{i}.ECGfilename.Text '_' num2str(j-1) '.amedianecg']));
                    % else
                    %     ECG = base64decode(BEATS{j}.medianecg.Text);
                    %     DATAA{iEcga}.ECG = readECGDATA(BEATS{j}.medianecg.Text);
                    end
                    DATAA{iEcga}.ECGfilename = CASEDATA{i}.ECGfilename.Text; 
                    zerotims = removeText(BEATS{j}.medianbaselinePoints);
                    BASEL=spline(zerotims,DATAA{iEcga}.ECG(:,zerotims),1:size(DATAA{iEcga}.ECG,2));
                    DATAA{iEcga}.ECG = DATAA{iEcga}.ECG - BASEL;                    
                    DATAA{iEcga}.beats = readResult(BEATS{j},DATAA{iEcga}.ECG ); 
                end
            end
        end
    end
%     if ~strcmp(DATAA{iEcga}.ECGfilename,DATAV{iEcgv}.ECGfilename)
%          disp(['missing beat in median ECG file: ' CASEDATA{i}.ECGfilename.Text])
%     end
end

if ~exist('DATAA') || ( ~exist('DATAA') && isempty(DATAA))
    DATA.VENTRICULAR = DATAV;
else
    DATA.ATRIAL = DATAA;
    DATA.VENTRICULAR = DATAV;
end
end
%%
function DATA = readCaseData(CASEDATA, ECGDIR)

iEcgv=0;
iEcga=0;

for i=1:length(CASEDATA)    
%     disp(['reading median ECG file: ' CASEDATA{i}.ECGfilename.Text])   
    
    if exist(fullfile(ECGDIR,[CASEDATA{i}.ECGfilename.Text '.ecg']),'file')
        ECG = loadmat(fullfile(ECGDIR,[CASEDATA{i}.ECGfilename.Text '.ecg']));
    elseif exist(fullfile(ECGDIR,CASEDATA{i}.ECGfilename.Text),'file')
        ECG = loadmat(fullfile(ECGDIR, CASEDATA{i}.ECGfilename.Text ));
    end
    zerotims = removeText( CASEDATA{i}.baselinePoints);
    if isfield(CASEDATA{i},'selectedVentricularBeats')
    
        if isfield(CASEDATA{i}.selectedVentricularBeats,'MEANTSIDIAGNOSTIC')
            BEATS = CASEDATA{i}.selectedVentricularBeats.MEANTSIDIAGNOSTIC;
        else
            BEATS = CASEDATA{i}.selectedVentricularBeats.CIPS_MEANTSIRESULT;
        end

        if length(BEATS) ==1
            if isfield(BEATS,'beatCineECG') && length(BEATS.beatCineECG.Text) > 1
                iEcgv=iEcgv+1;
                DATAV{iEcgv}.ECGfilename = CASEDATA{i}.ECGfilename.Text;    
                BASEL=spline(zerotims,ECG(:,zerotims),1:size(ECG,2));
                DATAV{iEcgv}.ECG = ECG - BASEL;
                DATAV{iEcgv}.beats{1} = readResult(BEATS,DATAV{iEcgv}.ECG );
            end
        else
            k= 1;
            iEcgv=iEcgv+1;
            DATAV{iEcgv}.ECGfilename = CASEDATA{i}.ECGfilename.Text;    
            BASEL=spline(zerotims,ECG(:,zerotims),1:size(ECG,2));
            DATAV{iEcgv}.ECG = ECG - BASEL;
            for j= 1:length(BEATS)
                if isfield(BEATS{j},'beatCineECG') && length(BEATS{j}.beatCineECG.Text) > 1
                    DATAV{iEcgv}.beats{k} = readResult(BEATS{j},DATAV{iEcgv}.ECG );
                    k=k+1;
                end
            end
        end
    end
    if isfield(CASEDATA{i},'selectedAtrialBeats')
        if isfield(CASEDATA{i}.selectedAtrialBeats,'MEANTSIDIAGNOSTIC')
            BEATS = CASEDATA{i}.selectedAtrialBeats.MEANTSIDIAGNOSTIC;
        else
            BEATS = CASEDATA{i}.selectedAtrialBeats.CIPS_MEANTSIRESULT;
        end

        if length(BEATS) ==1
            if isfield(BEATS,'beatCineECG') && length(BEATS.beatCineECG.Text) > 1
                iEcga=iEcga+1;    
                BASEL=spline(zerotims,ECG(:,zerotims),1:size(ECG,2));
                DATAA{iEcga}.ECG = ECG - BASEL;       
                DATAA{iEcga}.ECGfilename = CASEDATA{i}.ECGfilename.Text;    
                DATAA{iEcga}.beats{1} = readResult(BEATS,DATAA{iEcga}.ECG );
                if DATAA{iEcga}.beats{1}.onsetP < zerotims(1) - 40
                    pzerotims = [DATAA{iEcga}.beats{1}.onsetP zerotims];
                    BASEL=spline(pzerotims,ECG(:,pzerotims),1:size(ECG,2));
                    DATAA{iEcga}.ECG = ECG - BASEL;                
                end
            end
        else
            k= 1;
            iEcga=iEcga+1;    
            DATAA{iEcga}.ECGfilename = CASEDATA{i}.ECGfilename.Text;    
            BASEL=spline(zerotims,ECG(:,zerotims),1:size(ECG,2));
            DATAA{iEcga}.ECG = ECG - BASEL;
            pzerotims = zerotims;
            for j= 1:length(BEATS)
                if isfield(BEATS{j},'beatCineECG') && length(BEATS{j}.beatCineECG.Text) > 1 && isfield(DATAA{iEcga},'beats')
                    pzerotims = sort([DATAA{iEcga}.beats{1}.onsetP pzerotims]);
                    DATAA{iEcga}.beats{k} = readResult(BEATS{j},DATAA{iEcga}.ECG );
                    k=k+1;
                end
            end
            pzerotims(diff(pzerotims) < 10) = [];
            BASEL=spline(zerotims,ECG(:,pzerotims),1:size(ECG,2));
            DATAA{iEcga}.ECG = ECG - BASEL;
        end
    end
end

if ~exist('DATAA') || ( ~exist('DATAA') && isempty(DATAA))
    DATA = DATAV;
else
    DATA.ATRIAL = DATAA;
    DATA.VENTRICULAR = DATAV;
end
end

%%
function SPECS = readResult(BEAT,ECG)

SPECS=[];
if isfield(BEAT, 'Pwave_onset' ), SPECS.onsetP           = removeText(BEAT.Pwave_onset); end
if isfield(BEAT, 'endPwave' ),    SPECS.endP             = removeText(BEAT.endPwave);    end
if isfield(BEAT, 'onsetQRS' ),    SPECS.onsetqrs         = removeText(BEAT.onsetQRS);    end
if isfield(BEAT, 'Twave_end' ),   SPECS.endtwave         = removeText(BEAT.Twave_end);   end
if isfield(BEAT, 'EndQRS' ),      SPECS.time_Jpoint      = removeText(BEAT.EndQRS);      end
if isfield(BEAT, 'PeakTwave' ),   SPECS.time_apexT       = removeText(BEAT.PeakTwave);   end
if isfield(SPECS,'onsetqrs') && SPECS.onsetqrs>0
    SPECS.qrsduration           = SPECS.time_Jpoint - SPECS.onsetqrs;
    SPECS.qrstduration          = SPECS.endtwave - SPECS.onsetqrs;
    if ~isempty(ECG)
        SPECS.ECGbeat               = ECG(:,SPECS.onsetqrs:SPECS.endtwave);
        SPECS.ECGqrs                = ECG(:,SPECS.onsetqrs:SPECS.time_Jpoint);
        SPECS.ECGtwave              = ECG(:,SPECS.time_Jpoint:SPECS.endtwave);
    end
else
    SPECS.PwaveDuration         = SPECS.endP - SPECS.onsetP;
end

if isfield(BEAT, 'beatCineECG'), SPECS.beatCineECG         = read3DVertices(BEAT.beatCineECG); end
if isfield(BEAT, 'beatVcg'),SPECS.beatVcg                  = read3DVertices(BEAT.beatVcg);end
if isfield(BEAT, 'beatCineECGHeart'),SPECS.beatCineECGHeart= read3DVertices(BEAT.beatCineECGHeart);end
if isfield(BEAT, 'beatVcgHeart'),SPECS.beatVcgHeart        = read3DVertices(BEAT.beatVcgHeart);end
if isfield(BEAT, 'meanQRSaxis'),SPECS.meanQRSaxis          = read3DVertices(BEAT.meanQRSaxis);end
if isfield(BEAT, 'meanQRSaxisPosition'),SPECS.meanQRSaxisPosition = read3DVertices(BEAT.meanQRSaxisPosition);end
if isfield(BEAT, 'truthorigin'),SPECS.truthorigin          = read3DVertices(BEAT.truthorigin);end
% if isfield(BEAT, 'velocityprofile'), SPECS.velocityprofile = removeText(BEAT.velocityprofile);end

if isfield(BEAT,'area') 
    a= removeText(BEAT.area);
    if ~isnan(a)
        SPECS.normalized_area       = removeText(BEAT.normalized_area);
        SPECS.area                  = removeText(BEAT.area);
    end
end

end
%%
function MODELDATA = readModelData(ROOTDATA)

% Model extraction from iECG file.
if isfield(ROOTDATA, 'volumemodel_group')
    if isfield(ROOTDATA.volumemodel_group, 'ventricles')
        
        % 1. Base 64 decode.
        ventricles = base64decode(ROOTDATA.volumemodel_group.ventricles.Text);
        
        
        % # vertices, 32 bit bin array.
        binarr = double([de2bi(ventricles(5),8, 'left-msb') de2bi(ventricles(6),8, 'left-msb') ...
            de2bi(ventricles(7),8, 'left-msb') de2bi(ventricles(8),8, 'left-msb') ]);
        lenghtVER = bi2de(binarr(10:end), 'left-msb');
        
        % Retrieve vertices, 64 bit (double precision). Saved: X1,
        % Y1, Z1, X2, Y2, Z2, .... XlengthVER, YlengthVER,
        % ZlengthVER.
        possiex = 8;
        possiey = 16;
        possiez = 24;
        VER = zeros(lenghtVER,3);
        for ii = 1:lenghtVER
            % X value
            binarr = double([de2bi(ventricles(possiex+1),8, 'left-msb') de2bi(ventricles(possiex+2),8, 'left-msb') ...
                de2bi(ventricles(possiex+3),8, 'left-msb') de2bi(ventricles(possiex+4),8, 'left-msb') ...
                de2bi(ventricles(possiex+5),8, 'left-msb') de2bi(ventricles(possiex+6),8, 'left-msb') ...
                de2bi(ventricles(possiex+7),8, 'left-msb') de2bi(ventricles(possiex+8),8, 'left-msb')]);
            [~,~,~, VER(ii,1)]=ieee754(binarr);
            possiex = possiex + 24;
            
            % Y value
            binarr = double([de2bi(ventricles(possiey+1),8, 'left-msb') de2bi(ventricles(possiey+2),8, 'left-msb') ...
                de2bi(ventricles(possiey+3),8, 'left-msb') de2bi(ventricles(possiey+4),8, 'left-msb') ...
                de2bi(ventricles(possiey+5),8, 'left-msb') de2bi(ventricles(possiey+6),8, 'left-msb') ...
                de2bi(ventricles(possiey+7),8, 'left-msb') de2bi(ventricles(possiey+8),8, 'left-msb')]);
            [~,~,~, VER(ii,2)]=ieee754(binarr);
            possiey = possiey + 24;
            
            % Z value
            binarr = double([de2bi(ventricles(possiez+1),8, 'left-msb') de2bi(ventricles(possiez+2),8, 'left-msb') ...
                de2bi(ventricles(possiez+3),8, 'left-msb') de2bi(ventricles(possiez+4),8, 'left-msb') ...
                de2bi(ventricles(possiez+5),8, 'left-msb') de2bi(ventricles(possiez+6),8, 'left-msb') ...
                de2bi(ventricles(possiez+7),8, 'left-msb') de2bi(ventricles(possiez+8),8, 'left-msb')]);
            [~,~,~, VER(ii,3)]=ieee754(binarr);
            possiez = possiez + 24;
            
        end
        
        startITRI = possiex;
        
        % # triangles, 32 bit array
        binarr = double([de2bi(ventricles(startITRI+1),8, 'left-msb') de2bi(ventricles(startITRI+2),8, 'left-msb') ...
            de2bi(ventricles(startITRI+3),8, 'left-msb') de2bi(ventricles(startITRI+4),8, 'left-msb') ]);
        lenghtITRI = bi2de(binarr(10:end), 'left-msb');
        
        ITRIloc = possiex + 4;
        ITRI = zeros(lenghtITRI,3);
        for ii = 1:lenghtITRI
            % ITRI1
            binarr = double([de2bi(ventricles(ITRIloc+1),8, 'left-msb') de2bi(ventricles(ITRIloc+2),8, 'left-msb') ...
                de2bi(ventricles(ITRIloc+3),8, 'left-msb') de2bi(ventricles(ITRIloc+4),8, 'left-msb')]);
            [~,~,ITRI(ii,1),~]=ieee754(binarr);
            ITRIloc = ITRIloc + 4;
            
            % ITRI2
            binarr = double([de2bi(ventricles(ITRIloc+1),8, 'left-msb') de2bi(ventricles(ITRIloc+2),8, 'left-msb') ...
                de2bi(ventricles(ITRIloc+3),8, 'left-msb') de2bi(ventricles(ITRIloc+4),8, 'left-msb')]);
            [~,~,ITRI(ii,2),~]=ieee754(binarr);
            ITRIloc = ITRIloc + 4;
            
            % ITRI3
            binarr = double([de2bi(ventricles(ITRIloc+1),8, 'left-msb') de2bi(ventricles(ITRIloc+2),8, 'left-msb') ...
                de2bi(ventricles(ITRIloc+3),8, 'left-msb') de2bi(ventricles(ITRIloc+4),8, 'left-msb') ]);
            [~,~,ITRI(ii,3),~]=ieee754(binarr);
            ITRIloc = ITRIloc + 4;
            
        end
        ITRI = ITRI+1;
        
        % Modeloutput ventricles
        MODELDATA.VER = VER;
        MODELDATA.ITRI = ITRI;
        clearvars VER ITRI ventricles
    end
    if isfield(ROOTDATA.volumemodel_group, 'ventriclestypes')
        MODELDATA.ventTYP = strsplit(ROOTDATA.volumemodel_group.ventriclestypes.Text, ',');
        MODELDATA.ventTYP = str2double(MODELDATA.ventTYP);
    end
    if isfield(ROOTDATA.volumemodel_group, 'leftblood')
        
        % 1. Base 64 decode.
        leftblood = base64decode(ROOTDATA.volumemodel_group.leftblood.Text);
        
        % # vertices, 32 bit bin array.
        binarr = double([de2bi(leftblood(5),8, 'left-msb') de2bi(leftblood(6),8, 'left-msb') ...
            de2bi(leftblood(7),8, 'left-msb') de2bi(leftblood(8),8, 'left-msb') ]);
        lenghtVER = bi2de(binarr(10:end), 'left-msb');
        
        % Retrieve vertices, 64 bit (double precision). Saved: X1,
        % Y1, Z1, X2, Y2, Z2, .... XlengthVER, YlengthVER,
        % ZlengthVER.
        possiex = 8;
        possiey = 16;
        possiez = 24;
        VER = zeros(lenghtVER,3);
        for ii = 1:lenghtVER
            % X value
            binarr = double([de2bi(leftblood(possiex+1),8, 'left-msb') de2bi(leftblood(possiex+2),8, 'left-msb') ...
                de2bi(leftblood(possiex+3),8, 'left-msb') de2bi(leftblood(possiex+4),8, 'left-msb') ...
                de2bi(leftblood(possiex+5),8, 'left-msb') de2bi(leftblood(possiex+6),8, 'left-msb') ...
                de2bi(leftblood(possiex+7),8, 'left-msb') de2bi(leftblood(possiex+8),8, 'left-msb')]);
            [~,~,~, VER(ii,1)]=ieee754(binarr);
            possiex = possiex + 24;
            
            % Y value
            binarr = double([de2bi(leftblood(possiey+1),8, 'left-msb') de2bi(leftblood(possiey+2),8, 'left-msb') ...
                de2bi(leftblood(possiey+3),8, 'left-msb') de2bi(leftblood(possiey+4),8, 'left-msb') ...
                de2bi(leftblood(possiey+5),8, 'left-msb') de2bi(leftblood(possiey+6),8, 'left-msb') ...
                de2bi(leftblood(possiey+7),8, 'left-msb') de2bi(leftblood(possiey+8),8, 'left-msb')]);
            [~,~,~, VER(ii,2)]=ieee754(binarr);
            possiey = possiey + 24;
            
            % Z value
            binarr = double([de2bi(leftblood(possiez+1),8, 'left-msb') de2bi(leftblood(possiez+2),8, 'left-msb') ...
                de2bi(leftblood(possiez+3),8, 'left-msb') de2bi(leftblood(possiez+4),8, 'left-msb') ...
                de2bi(leftblood(possiez+5),8, 'left-msb') de2bi(leftblood(possiez+6),8, 'left-msb') ...
                de2bi(leftblood(possiez+7),8, 'left-msb') de2bi(leftblood(possiez+8),8, 'left-msb')]);
            [~,~,~, VER(ii,3)]=ieee754(binarr);
            possiez = possiez + 24;
            
        end
        
        startITRI = possiex;
        
        % # triangles, 32 bit array
        binarr = double([de2bi(leftblood(startITRI+1),8, 'left-msb') de2bi(leftblood(startITRI+2),8, 'left-msb') ...
            de2bi(leftblood(startITRI+3),8, 'left-msb') de2bi(leftblood(startITRI+4),8, 'left-msb') ]);
        lenghtITRI = bi2de(binarr(10:end), 'left-msb');
        
        ITRIloc = possiex + 4;
        ITRI = zeros(lenghtITRI,3);
        for ii = 1:lenghtITRI
            % ITRI1
            binarr = double([de2bi(leftblood(ITRIloc+1),8, 'left-msb') de2bi(leftblood(ITRIloc+2),8, 'left-msb') ...
                de2bi(leftblood(ITRIloc+3),8, 'left-msb') de2bi(leftblood(ITRIloc+4),8, 'left-msb')]);
            [~,~,ITRI(ii,1),~]=ieee754(binarr);
            ITRIloc = ITRIloc + 4;
            
            % ITRI2
            binarr = double([de2bi(leftblood(ITRIloc+1),8, 'left-msb') de2bi(leftblood(ITRIloc+2),8, 'left-msb') ...
                de2bi(leftblood(ITRIloc+3),8, 'left-msb') de2bi(leftblood(ITRIloc+4),8, 'left-msb')]);
            [~,~,ITRI(ii,2),~]=ieee754(binarr);
            ITRIloc = ITRIloc + 4;
            
            % ITRI3
            binarr = double([de2bi(leftblood(ITRIloc+1),8, 'left-msb') de2bi(leftblood(ITRIloc+2),8, 'left-msb') ...
                de2bi(leftblood(ITRIloc+3),8, 'left-msb') de2bi(leftblood(ITRIloc+4),8, 'left-msb') ]);
            [~,~,ITRI(ii,3),~]=ieee754(binarr);
            ITRIloc = ITRIloc + 4;
            
        end
        ITRI = ITRI+1;
        
        % Modeloutput ventricles
        MODELDATA.LcavVER = VER;
        MODELDATA.LcavITRI = ITRI;
        clearvars VER ITRI leftblood
        
    end
    if isfield(ROOTDATA.volumemodel_group, 'leftbloodtypes')
        MODELDATA.LTYP = strsplit(ROOTDATA.volumemodel_group.leftbloodtypes.Text, ',');
        MODELDATA.LTYP = str2double(MODELDATA.LTYP);
    end
    if isfield(ROOTDATA.volumemodel_group, 'rightblood')
        
        % 1. Base 64 decode.
        rightblood = base64decode(ROOTDATA.volumemodel_group.rightblood.Text);
        
        % # vertices, 32 bit bin array.
        binarr = double([de2bi(rightblood(5),8, 'left-msb') de2bi(rightblood(6),8, 'left-msb') ...
            de2bi(rightblood(7),8, 'left-msb') de2bi(rightblood(8),8, 'left-msb') ]);
        lenghtVER = bi2de(binarr(10:end), 'left-msb');
        
        % Retrieve vertices, 64 bit (double precision). Saved: X1,
        % Y1, Z1, X2, Y2, Z2, .... XlengthVER, YlengthVER,
        % ZlengthVER.
        possiex = 8;
        possiey = 16;
        possiez = 24;
        VER = zeros(lenghtVER,3);
        for ii = 1:lenghtVER
            % X value
            binarr = double([de2bi(rightblood(possiex+1),8, 'left-msb') de2bi(rightblood(possiex+2),8, 'left-msb') ...
                de2bi(rightblood(possiex+3),8, 'left-msb') de2bi(rightblood(possiex+4),8, 'left-msb') ...
                de2bi(rightblood(possiex+5),8, 'left-msb') de2bi(rightblood(possiex+6),8, 'left-msb') ...
                de2bi(rightblood(possiex+7),8, 'left-msb') de2bi(rightblood(possiex+8),8, 'left-msb')]);
            [~,~,~, VER(ii,1)]=ieee754(binarr);
            possiex = possiex + 24;
            
            % Y value
            binarr = double([de2bi(rightblood(possiey+1),8, 'left-msb') de2bi(rightblood(possiey+2),8, 'left-msb') ...
                de2bi(rightblood(possiey+3),8, 'left-msb') de2bi(rightblood(possiey+4),8, 'left-msb') ...
                de2bi(rightblood(possiey+5),8, 'left-msb') de2bi(rightblood(possiey+6),8, 'left-msb') ...
                de2bi(rightblood(possiey+7),8, 'left-msb') de2bi(rightblood(possiey+8),8, 'left-msb')]);
            [~,~,~, VER(ii,2)]=ieee754(binarr);
            possiey = possiey + 24;
            
            % Z value
            binarr = double([de2bi(rightblood(possiez+1),8, 'left-msb') de2bi(rightblood(possiez+2),8, 'left-msb') ...
                de2bi(rightblood(possiez+3),8, 'left-msb') de2bi(rightblood(possiez+4),8, 'left-msb') ...
                de2bi(rightblood(possiez+5),8, 'left-msb') de2bi(rightblood(possiez+6),8, 'left-msb') ...
                de2bi(rightblood(possiez+7),8, 'left-msb') de2bi(rightblood(possiez+8),8, 'left-msb')]);
            [~,~,~, VER(ii,3)]=ieee754(binarr);
            possiez = possiez + 24;
            
        end
        
        startITRI = possiex;
        
        % # triangles, 32 bit array
        binarr = double([de2bi(rightblood(startITRI+1),8, 'left-msb') de2bi(rightblood(startITRI+2),8, 'left-msb') ...
            de2bi(rightblood(startITRI+3),8, 'left-msb') de2bi(rightblood(startITRI+4),8, 'left-msb') ]);
        lenghtITRI = bi2de(binarr(10:end), 'left-msb');
        
        ITRIloc = possiex + 4;
        ITRI = zeros(lenghtITRI,3);
        for ii = 1:lenghtITRI
            % ITRI1
            binarr = double([de2bi(rightblood(ITRIloc+1),8, 'left-msb') de2bi(rightblood(ITRIloc+2),8, 'left-msb') ...
                de2bi(rightblood(ITRIloc+3),8, 'left-msb') de2bi(rightblood(ITRIloc+4),8, 'left-msb')]);
            [~,~,ITRI(ii,1),~]=ieee754(binarr);
            ITRIloc = ITRIloc + 4;
            
            % ITRI2
            binarr = double([de2bi(rightblood(ITRIloc+1),8, 'left-msb') de2bi(rightblood(ITRIloc+2),8, 'left-msb') ...
                de2bi(rightblood(ITRIloc+3),8, 'left-msb') de2bi(rightblood(ITRIloc+4),8, 'left-msb')]);
            [~,~,ITRI(ii,2),~]=ieee754(binarr);
            ITRIloc = ITRIloc + 4;
            
            % ITRI3
            binarr = double([de2bi(rightblood(ITRIloc+1),8, 'left-msb') de2bi(rightblood(ITRIloc+2),8, 'left-msb') ...
                de2bi(rightblood(ITRIloc+3),8, 'left-msb') de2bi(rightblood(ITRIloc+4),8, 'left-msb') ]);
            [~,~,ITRI(ii,3),~]=ieee754(binarr);
            ITRIloc = ITRIloc + 4;
            
        end
        ITRI = ITRI+1;
        
        % Modeloutput ventricles
        MODELDATA.RcavVER = VER;
        MODELDATA.RcavITRI = ITRI;
        clearvars VER ITRI rightblood
    end
    if isfield(ROOTDATA.volumemodel_group, 'rightbloodtypes')
        MODELDATA.RTYP = strsplit(ROOTDATA.volumemodel_group.rightbloodtypes.Text, ',');
        MODELDATA.RTYP = str2double(MODELDATA.RTYP);
    end
    if isfield(ROOTDATA.volumemodel_group, 'thorax')
        % 1. Base 64 decode.
        thorax = base64decode(ROOTDATA.volumemodel_group.thorax.Text);
        
        % # vertices, 32 bit bin array.
        binarr = double([de2bi(thorax(5),8, 'left-msb') de2bi(thorax(6),8, 'left-msb') ...
            de2bi(thorax(7),8, 'left-msb') de2bi(thorax(8),8, 'left-msb') ]);
        lenghtVER = bi2de(binarr(10:end), 'left-msb');
        
        % Retrieve vertices, 64 bit (double precision). Saved: X1,
        % Y1, Z1, X2, Y2, Z2, .... XlengthVER, YlengthVER,
        % ZlengthVER.
        possiex = 8;
        possiey = 16;
        possiez = 24;
        VER = zeros(lenghtVER,3);
        for ii = 1:lenghtVER
            % X value
            binarr = double([de2bi(thorax(possiex+1),8, 'left-msb') de2bi(thorax(possiex+2),8, 'left-msb') ...
                de2bi(thorax(possiex+3),8, 'left-msb') de2bi(thorax(possiex+4),8, 'left-msb') ...
                de2bi(thorax(possiex+5),8, 'left-msb') de2bi(thorax(possiex+6),8, 'left-msb') ...
                de2bi(thorax(possiex+7),8, 'left-msb') de2bi(thorax(possiex+8),8, 'left-msb')]);
            [~,~,~, VER(ii,1)]=ieee754(binarr);
            possiex = possiex + 24;
            
            % Y value
            binarr = double([de2bi(thorax(possiey+1),8, 'left-msb') de2bi(thorax(possiey+2),8, 'left-msb') ...
                de2bi(thorax(possiey+3),8, 'left-msb') de2bi(thorax(possiey+4),8, 'left-msb') ...
                de2bi(thorax(possiey+5),8, 'left-msb') de2bi(thorax(possiey+6),8, 'left-msb') ...
                de2bi(thorax(possiey+7),8, 'left-msb') de2bi(thorax(possiey+8),8, 'left-msb')]);
            [~,~,~, VER(ii,2)]=ieee754(binarr);
            possiey = possiey + 24;
            
            % Z value
            binarr = double([de2bi(thorax(possiez+1),8, 'left-msb') de2bi(thorax(possiez+2),8, 'left-msb') ...
                de2bi(thorax(possiez+3),8, 'left-msb') de2bi(thorax(possiez+4),8, 'left-msb') ...
                de2bi(thorax(possiez+5),8, 'left-msb') de2bi(thorax(possiez+6),8, 'left-msb') ...
                de2bi(thorax(possiez+7),8, 'left-msb') de2bi(thorax(possiez+8),8, 'left-msb')]);
            [~,~,~, VER(ii,3)]=ieee754(binarr);
            possiez = possiez + 24;
            
        end
        
        startITRI = possiex;
        
        % # triangles, 32 bit array
        binarr = double([de2bi(thorax(startITRI+1),8, 'left-msb') de2bi(thorax(startITRI+2),8, 'left-msb') ...
            de2bi(thorax(startITRI+3),8, 'left-msb') de2bi(thorax(startITRI+4),8, 'left-msb') ]);
        lenghtITRI = bi2de(binarr(10:end), 'left-msb');
        
        ITRIloc = possiex + 4;
        ITRI = zeros(lenghtITRI,3);
        for ii = 1:lenghtITRI
            % ITRI1
            binarr = double([de2bi(thorax(ITRIloc+1),8, 'left-msb') de2bi(thorax(ITRIloc+2),8, 'left-msb') ...
                de2bi(thorax(ITRIloc+3),8, 'left-msb') de2bi(thorax(ITRIloc+4),8, 'left-msb')]);
            [~,~,ITRI(ii,1),~]=ieee754(binarr);
            ITRIloc = ITRIloc + 4;
            
            % ITRI2
            binarr = double([de2bi(thorax(ITRIloc+1),8, 'left-msb') de2bi(thorax(ITRIloc+2),8, 'left-msb') ...
                de2bi(thorax(ITRIloc+3),8, 'left-msb') de2bi(thorax(ITRIloc+4),8, 'left-msb')]);
            [~,~,ITRI(ii,2),~]=ieee754(binarr);
            ITRIloc = ITRIloc + 4;
            
            % ITRI3
            binarr = double([de2bi(thorax(ITRIloc+1),8, 'left-msb') de2bi(thorax(ITRIloc+2),8, 'left-msb') ...
                de2bi(thorax(ITRIloc+3),8, 'left-msb') de2bi(thorax(ITRIloc+4),8, 'left-msb') ]);
            [~,~,ITRI(ii,3),~]=ieee754(binarr);
            ITRIloc = ITRIloc + 4;
            
        end
        ITRI = ITRI+1;
        
        % Modeloutput ventricles
        MODELDATA.TVER = VER;
        MODELDATA.TITRI = ITRI;
        clearvars VER ITRI thorax
    end
    if isfield(ROOTDATA, 'leadsystemmodel_group')
        if isfield(ROOTDATA.leadsystemmodel_group, 'leadsystems')
            if isfield(ROOTDATA.leadsystemmodel_group.leadsystems, 'geomleadsystemgroup')
                if isfield(ROOTDATA.leadsystemmodel_group.leadsystems.geomleadsystemgroup, 'electrodepositions')
                    MODELDATA.elecPOS = strsplit(ROOTDATA.leadsystemmodel_group.leadsystems.geomleadsystemgroup.electrodepositions.Text, ',');
                    MODELDATA.elecPOS = str2double(MODELDATA.elecPOS);
                    MODELDATA.elecPOS = reshape(MODELDATA.elecPOS, [3,9])';
                end
            end
        end
    end
    
end
end
%%
function ECG = readECGDATA(DATA)


tempFile = 'tmp.t1';
fid = fopen(tempFile, 'w');decompressedData = zlibdecode(base64decode(DATA));
fwrite(fid, decompressedData,'int8');
fclose(fid);
ECG = loadmat(tempFile );
delete(tempFile);

end
%%
function y = base64decode(x, outfname, alg)
%BASE64DECODE Perform base64 decoding on a string.
%
% INPUT:
%   x    - block of data to be decoded.  Can be a string or a numeric  
%          vector containing integers in the range 0-255. Any character
%          not part of the 65-character base64 subset set is silently
%          ignored.  Characters occuring after a '=' padding character are 
%          never decoded. If the length of the string to decode (after 
%          ignoring non-base64 chars) is not a multiple of 4, then a 
%          warning is generated.
%
%   outfname - if provided the binary date from decoded string will be
%          saved into a file. Since Base64 coding is often used to embbed
%          binary data in xml files, this option can be used to extract and
%          save them.
%
%   alg  - Algorithm to use: can take values 'java' or 'matlab'. Optional
%          variable defaulting to 'java' which is a little faster. If 
%          'java' is chosen than core of the code is performed by a call to
%          a java library. Optionally all operations can be performed using
%          matleb code. 
%
% OUTPUT:
%   y    - array of binary data returned as uint8 
%
%   This function is used to decode strings from the Base64 encoding specified
%   in RFC 2045 - MIME (Multipurpose Internet Mail Extensions).  The Base64
%   encoding is designed to represent arbitrary sequences of octets in a form
%   that need not be humanly readable.  A 65-character subset ([A-Za-z0-9+/=])
%   of US-ASCII is used, enabling 6 bits to be represented per printable
%   character.
%
%   See also BASE64ENCODE.
%
%   Written by Jarek Tuszynski, SAIC, jaroslaw.w.tuszynski_at_saic.com
%
%   Matlab version based on 2004 code by Peter J. Acklam
%   E-mail:      pjacklam@online.no
%   URL:         http://home.online.no/~pjacklam
%   http://home.online.no/~pjacklam/matlab/software/util/datautil/base64encode.m

if nargin<3, alg='java';  end
if nargin<2, outfname=''; end

%% if x happen to be a filename than read the file
if (numel(x)<256)
  if (exist(x, 'file')==2)
    fid = fopen(x,'rb');
    x = fread(fid, 'uint8');   
    fclose(fid);
  end
end
x = uint8(x(:)); % unify format

%% Perform conversion
switch (alg)
  case 'java' 
    base64 = org.apache.commons.codec.binary.Base64;
    y = base64.decode(x);
    y = mod(int16(y),256); % convert from int8 to uint8
  case 'matlab'
    %%  Perform the mapping
    %   A-Z  ->  0  - 25
    %   a-z  ->  26 - 51
    %   0-9  ->  52 - 61
    %   + -  ->  62       '-' is URL_SAFE alternative
    %   / _  ->  63       '_' is URL_SAFE alternative
    map = uint8(zeros(1,256)+65);
    map(uint8(['A':'Z', 'a':'z', '0':'9', '+/=']))= 0:64;
    map(uint8('-_'))= 62:63;  % URL_SAFE alternatives
    x = map(x);  % mapping
    
    x(x>64)=[]; % remove non-base64 chars
    if rem(numel(x), 4)
      warning('Length of base64 data not a multiple of 4; padding input.');
    end
    x(x==64)=[]; % remove padding characters
    
    %% add padding and reshape
    nebytes = length(x);         % number of encoded bytes
    nchunks = ceil(nebytes/4);   % number of chunks/groups
    if rem(nebytes, 4)>0
      x(end+1 : 4*nchunks) = 0;  % add padding
    end
    x = reshape(uint8(x), 4, nchunks);
    y = repmat(uint8(0), 3, nchunks);            % for the decoded data
    
    %% Rearrange every 4 bytes into 3 bytes
    %    00aaaaaa 00bbbbbb 00cccccc 00dddddd
    % to form
    %    aaaaaabb bbbbcccc ccdddddd
    y(1,:) = bitshift(x(1,:), 2);                 % 6 highest bits of y(1,:)
    y(1,:) = bitor(y(1,:), bitshift(x(2,:), -4)); % 2 lowest  bits of y(1,:)
    y(2,:) = bitshift(x(2,:), 4);                 % 4 highest bits of y(2,:)
    y(2,:) = bitor(y(2,:), bitshift(x(3,:), -2)); % 4 lowest  bits of y(2,:)
    y(3,:) = bitshift(x(3,:), 6);                 % 2 highest bits of y(3,:)
    y(3,:) = bitor(y(3,:), x(4,:));               % 6 lowest  bits of y(3,:)
    
    %% remove extra padding
    switch rem(nebytes, 4)
      case 2
        y = y(1:end-2);
      case 3
        y = y(1:end-1);
    end
end

%% reshape to a row vector and make it a character array
y = uint8(reshape(y, 1, numel(y)));

%% save to file if needed
if ~isempty(outfname)
  fid = fopen(outfname,'wb');
  fwrite(fid, y, 'uint8');  
  fclose(fid);
end

end
%%
function output = zlibdecode(input)
%ZLIBDECODE Decompress input bytes using ZLIB.
%
%    output = zlibdecode(input)
%
% The function takes a compressed byte array INPUT and returns inflated
% bytes OUTPUT. The INPUT is a result of GZIPENCODE function. The OUTPUT
% is always an 1-by-N uint8 array. JAVA must be enabled to use the function.
%
% See also zlibencode typecast

error(nargchk(1, 1, nargin));
error(javachk('jvm'));
if ischar(input)
  warning('zlibdecode:inputTypeMismatch', ...
          'Input is char, but treated as uint8.');
  input = uint8(input);
end
if ~isa(input, 'int8') && ~isa(input, 'uint8')
    error('Input must be either int8 or uint8.');
end

buffer = java.io.ByteArrayOutputStream();
zlib = java.util.zip.InflaterOutputStream(buffer);
zlib.write(input, 0, numel(input));
zlib.close();
output = typecast(buffer.toByteArray(), 'uint8')';

end
