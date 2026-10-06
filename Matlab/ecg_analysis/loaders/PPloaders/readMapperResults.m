function DATACASES = readMapperResults( varargin )

analyseBeats = 0;

if length(varargin) < 1
    error('This routine needs at least two parameters');
else
    fn = varargin{1};
    pp=2;
    while pp<=nargin
        if ischar(varargin{pp})
            key=lower(varargin{pp});
            switch key
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
if isempty(basedir)
    basedir = pwd;
end

XML=xml2struct(fn);

if ~isempty(XML.POLLYMAPPER{2}.POLLYDATA.ECG_cases.ECG_CASE) 
    POLLYDATA = XML.POLLYMAPPER{2}.POLLYDATA;
    POLLYCASE = POLLYDATA.ECG_cases.ECG_CASE;
    DATACASES.patientid=name;
    if isfield(POLLYDATA,'modelname')
        DATACASES.modelName = POLLYDATA.modelname.Text;
    end
    if isfield(POLLYDATA,'leadsystemmodel_group')
        if size(POLLYDATA.leadsystemmodel_group.leadsystems,2) ==1
            leads = removeText(POLLYDATA.leadsystemmodel_group.leadsystems.geomleadsystemgroup.electrodepositions);
            DATACASES.elecs = reshape(leads,3,length(leads)/3)';
        else
            useLead = removeText(POLLYDATA.leadsystemmodel_group.currentLeadSystem)+1;            
            leads = removeText(POLLYDATA.leadsystemmodel_group.leadsystems{useLead}.geomleadsystemgroup.electrodepositions);
            DATACASES.elecs = reshape(leads,3,length(leads)/3)';          
        end
    end
    
    if analyseBeats == 0    
        DATACASES.DATA = readCaseData(basedir,POLLYCASE);
    elseif analyseBeats == 1
        DATACASES.DATA = readMedianData(basedir,POLLYCASE);
    end
    
    DATACASES.name = POLLYCASE.ECGcasename.Text;
    
else
    for k=1:length(XML.POLLYMAPPER{2}.POLLYDATA.ECG_cases.ECG_CASE)
        POLLYCASE = XML.POLLYMAPPER{2}.POLLYDATA.ECG_cases.ECG_CASE{k};
        DATA = readCaseData(fullfile(basedir,XML.POLLY{2}.POLLYDATA.patientId.Text),POLLYCASE);
        if ~isempty(DATA)
            DATACASES.DATA = DATA;
            DATACASES.patientid=name;
            DATACASES.name = POLLYCASE.ECGcasename.Text;            
        end
    end
end


%%
function VER = read3DVertices(str)
V=removeText(str);
if length(V) < 3 
    VER=[];
else
    VER = reshape(V,3,length(V)/3)';
end



%%
function val= removeText(str)
if isfield( str,'Text')
    val = str2double(strsplit(str.Text,','));
else
    val = str2double(struct2array(str));
end

%% selected beats
function DATA = readCaseData(basedir,POLLYCASE)

DATA = struct(); 
if isfield(POLLYCASE,'ECGfiles')
    for iEcg=1:length(POLLYCASE.ECGfiles.ECGFiledata)
        if length(POLLYCASE.ECGfiles.ECGFiledata) < 2
            ecg=POLLYCASE.ECGfiles.ECGFiledata;
        else
            ecg=POLLYCASE.ECGfiles.ECGFiledata{iEcg};
        end
        
        % Zamiana nazwy pliku na bezpieczną nazwę zmiennej w strukturze MATLABa
        rawFileName = ecg.ECGfilename.Text;
        [~, fnameBase, ~] = fileparts(rawFileName);
        fieldName = regexprep(rawFileName, '\W', '_');
        if ~isletter(fieldName(1))
            fieldName = ['sig_', fieldName];
        end
        
        if length(rawFileName) >= 4 && (strcmp(rawFileName(end-3:end),'.ecg') || strcmp(rawFileName(end-3:end),'.bsm'))
            if exist( fullfile(basedir,'ECG_DATA',rawFileName),'file')
                ECG = loadmat(fullfile(basedir,'ECG_DATA',rawFileName));
            else
                ECG=[];
            end
        else
            if exist(fullfile(basedir,'ECG_DATA',[rawFileName '.ecg']), 'file')
                ECG = loadmat(fullfile(basedir,'ECG_DATA',[rawFileName '.ecg']));
            else
                ECG=[];
            end
        end
        
        % Przypisanie bezpośrednio do fieldName (obok beats)
        DATA.(fieldName).ECG = ECG;
        DATA.(fieldName).filename = rawFileName;
        
        if isfield(ecg,'autointerpret')
            DATA.(fieldName).autointerpret = ecg.autointerpret.Text;
        end
        
        if isfield(ecg,'selectedVentricularBeats')
            beats = ecg.selectedVentricularBeats.CyncRESULT;
            ecg = rmfield(ecg,[{'useECGSignal'} {'selectedVentricularBeats'} {'Attributes'}]);
            
            i=0;
            for k=1:length(beats)
                if length(beats) < 2
                    selbeat = beats;
                else
                    selbeat = beats{k};
                end
                
                if isfield(selbeat,'initdep') || isfield(selbeat,'finaldep')
                    i=i+1;
                    SPECS = struct();
                    
                    if isfield(selbeat,'Pwave_onset')
                        SPECS.onsetP = removeText(selbeat.Pwave_onset);
                    end
                    if isfield(selbeat,'onsetQRS')
                        SPECS.onsetqrs = removeText(selbeat.onsetQRS);
                    else
                        SPECS.onsetqrs = 1;
                    end
                    if isfield(selbeat,'Twave_end')
                        SPECS.endtwave = removeText(selbeat.Twave_end);
                    else
                        SPECS.endtwave = size(ECG, 2);
                    end
                    if isfield(selbeat,'EndQRS')
                        SPECS.time_Jpoint = removeText(selbeat.EndQRS);
                    else
                        SPECS.time_Jpoint = SPECS.onsetqrs;
                    end
                    if isfield(selbeat,'PeakTwave')
                        SPECS.time_apexT = removeText(selbeat.PeakTwave);
                    end
                    
                    SPECS.qrsduration   = SPECS.time_Jpoint - SPECS.onsetqrs;
                    SPECS.qrstduration  = SPECS.endtwave - SPECS.onsetqrs;
                    
                    DATA.(fieldName).beats{i}.SPECS = SPECS;
                    
                    if isfield(selbeat,'initdep')
                        DATA.(fieldName).beats{i}.initdep = removeText(selbeat.initdep)';
                    end
                    if isfield(selbeat,'initdepvelocity')
                        DATA.(fieldName).beats{i}.initdepvelocity = removeText(selbeat.initdepvelocity);
                    end
                    
                    if isfield(selbeat,'finaldep')
                        DATA.(fieldName).beats{i}.finaldep = removeText(selbeat.finaldep)';
                    end
                    if isfield(selbeat,'initialcorrelation')
                        DATA.(fieldName).beats{i}.initialcorrelation = removeText(selbeat.initialcorrelation)';
                    end
                    if isfield(selbeat,'finalcorrelation')
                        DATA.(fieldName).beats{i}.finalcorrelation = removeText(selbeat.finalcorrelation)';
                    end
                    
                    if isfield(selbeat,'finalrep')
                        DATA.(fieldName).beats{i}.finalrep = removeText(selbeat.finalrep)';
                        DATA.(fieldName).beats{i}.rep      = removeText(selbeat.finalrep)';
                    end
                    
                    if ( isfield(selbeat,'massActivated') )
                        DATA.(fieldName).beats{i}.massActivated = removeText(selbeat.massActivated)';
                    end
                    if ( isfield(selbeat,'massLVseptum') )
                        DATA.(fieldName).beats{i}.massLVseptum = removeText(selbeat.massLVseptum)';
                    end
                    if ( isfield(selbeat,'massLVanterior') )
                        DATA.(fieldName).beats{i}.massLVanterior = removeText(selbeat.massLVanterior)';
                    end
                    if ( isfield(selbeat,'massLVposterior') )
                        DATA.(fieldName).beats{i}.massLVposterior = removeText(selbeat.massLVposterior)';
                    end
                    if ( isfield(selbeat,'massRVseptum') )
                        DATA.(fieldName).beats{i}.massRVseptum = removeText(selbeat.massRVseptum)';
                    end
                    if ( isfield(selbeat,'massRVfreewall') )
                        DATA.(fieldName).beats{i}.massRVfreewall = removeText(selbeat.massRVfreewall)';
                    end
                    if ( isfield(selbeat,'massLVActivated') )
                        DATA.(fieldName).beats{i}.massLVActivated = removeText(selbeat.massLVActivated)';
                    end
                    if ( isfield(selbeat,'massRVActivated') )
                        DATA.(fieldName).beats{i}.massRVActivated = removeText(selbeat.massRVActivated)';
                    end
                    if ( isfield(selbeat,'massSEActivated') )
                        DATA.(fieldName).beats{i}.massSEActivated = removeText(selbeat.massSEActivated)';
                    end
                    
                    if isfield(selbeat,'truthorigin')
                        DATA.(fieldName).beats{i}.truthorigin = removeText(selbeat.truthorigin)';
                    elseif isfield(selbeat,'DEPREPtruthLocation') 
                        DATA.(fieldName).beats{i}.truthorigin = removeText(selbeat.DEPREPtruthLocation)';
                    end

                    if ( isfield(selbeat,'DEPREPbeatMeanTSI') )
                        DATA.(fieldName).beats{i}.beatMeanTSI = read3DVertices(selbeat.DEPREPbeatMeanTSI);
                    elseif ( isfield(selbeat,'beatMeanTSI') )
                        DATA.(fieldName).beats{i}.beatMeanTSI = read3DVertices(selbeat.beatMeanTSI);
                    end
                    if ( isfield(selbeat,'DEPREPbeatMeanTSIHeart') )
                        DATA.(fieldName).beats{i}.beatMeanTSIHeart = read3DVertices(selbeat.DEPREPbeatMeanTSIHeart);
                    elseif ( isfield(selbeat,'beatMeanTSIHeart') )
                        DATA.(fieldName).beats{i}.beatMeanTSIHeart = read3DVertices(selbeat.beatMeanTSIHeart);
                    end
                    if ( isfield(selbeat,'DEPREPbeatVcgHeart') )
                        DATA.(fieldName).beats{i}.beatVcgHeart = read3DVertices(selbeat.DEPREPbeatVcgHeart);
                    elseif ( isfield(selbeat,'beatVcgHeart') )
                        DATA.(fieldName).beats{i}.beatVcgHeart = read3DVertices(selbeat.beatVcgHeart);
                    end
                    if ( isfield(selbeat,'DEPREPmeanQRSaxisPosition') )
                        DATA.(fieldName).beats{i}.meanQRSaxis = read3DVertices(selbeat.DEPREPmeanQRSaxisPosition);
                    elseif ( isfield(selbeat,'meanQRSaxis') )
                        DATA.(fieldName).beats{i}.meanQRSaxis = read3DVertices(selbeat.meanQRSaxis);
                    end
                    if ( isfield(selbeat,'DEPREPmeanQRSaxisPosition') )
                        DATA.(fieldName).beats{i}.meanQRSaxisPosition = read3DVertices(selbeat.DEPREPmeanQRSaxisPosition);
                    elseif ( isfield(selbeat,'meanQRSaxisPosition') )
                        DATA.(fieldName).beats{i}.meanQRSaxisPosition = read3DVertices(selbeat.meanQRSaxisPosition);
                    end
                    
                    if ~isempty(ECG)
                        sIdx = max(1, SPECS.onsetqrs);
                        eIdx = min(SPECS.endtwave, size(ECG,2));
                        if isempty(eIdx) || eIdx < sIdx
                            eIdx = size(ECG,2);
                        end
                        DATA.(fieldName).beats{i}.ECG = ECG(:, sIdx:eIdx);
                    end
                end
            end
        end
        
        if isfield(ecg,'selectedAtrialBeats')
            beats = ecg.selectedAtrialBeats.CyncRESULT;
            i=0;
            for k=1:length(beats)
                if isfield(beats{k},'finaldep') && isfield(beats{k},'finalrep')
                    i=i+1;
                    selbeat = beats{k};
                    
                    SPECS = struct();
                    if isfield(selbeat, 'Pwave_onset')
                        SPECS.onsetP = removeText(selbeat.Pwave_onset);
                    else
                        SPECS.onsetP = 1;
                    end
                    if isfield(selbeat, 'onsetQRS')
                        SPECS.onsetqrs = removeText(selbeat.onsetQRS);
                    else
                        SPECS.onsetqrs = 1;
                    end
                    if isfield(selbeat, 'Twave_end')
                        SPECS.endtwave = removeText(selbeat.Twave_end);
                    else
                        SPECS.endtwave = size(ECG, 2);
                    end
                    if isfield(selbeat, 'EndQRS')
                        SPECS.time_Jpoint = removeText(selbeat.EndQRS);
                    else
                        SPECS.time_Jpoint = SPECS.onsetqrs;
                    end
                    if isfield(selbeat, 'PeakTwave')
                        SPECS.time_apexT = removeText(selbeat.PeakTwave);
                    end
                    
                    SPECS.time_apexU    = -1;
                    SPECS.depSlope      = 2;
                    
                    if isfield(selbeat, 'initialrepslope')
                        SPECS.initialSlope  = removeText(selbeat.initialrepslope);
                    end
                    if isfield(selbeat, 'platslope')
                        SPECS.plateauslope  = removeText(selbeat.platslope);
                    else
                        SPECS.plateauslope = 0.014;
                    end
                    if SPECS.plateauslope == 0
                        SPECS.plateauslope = 0.014;
                    end
                    
                    if isfield(selbeat, 'repslope')
                        SPECS.repslope  = removeText(selbeat.repslope);
                    else
                        SPECS.repslope = 0.045;
                    end
                    if SPECS.repslope == 0
                        SPECS.repslope = 0.045;
                    end
                    
                    SPECS.repCorrection = 0;
                    SPECS.useCumsum     = 0;
                    SPECS.qrsduration   = SPECS.time_Jpoint - SPECS.onsetqrs;
                    SPECS.qrstduration  = SPECS.endtwave - SPECS.onsetqrs;
                    
                    DATA.(fieldName).beats{i}.finaldep = removeText(beats{k}.finaldep)';
                    DATA.(fieldName).beats{i}.dep      = removeText(beats{k}.finaldep)'; 
                    DATA.(fieldName).beats{i}.finalrep = removeText(beats{k}.finalrep)';
                    DATA.(fieldName).beats{i}.rep      = removeText(beats{k}.finalrep)';
                    DATA.(fieldName).beats{i}.SPECS    = SPECS;
                    
                    sIdx = max(1, SPECS.onsetP);
                    eIdx = min(SPECS.endtwave, size(ECG,2));
                    if isempty(eIdx) || eIdx < sIdx
                        eIdx = size(ECG,2);
                    end
                    DATA.(fieldName).beats{i}.ECG = ECG(:, sIdx:eIdx);
                end
            end
        end
    end
else
    disp('no analysed data found')
end

%% selected beats (median)
function DATA = readMedianData(basedir,POLLYCASE)

DATA = struct(); 
if isfield(POLLYCASE,'ECGfiles')
    for iE=1:length(POLLYCASE.ECGfiles.ECGFiledata)
        if  length(POLLYCASE.ECGfiles.ECGFiledata) < 2
            ecg=POLLYCASE.ECGfiles.ECGFiledata;
        else
            ecg=POLLYCASE.ECGfiles.ECGFiledata{iE};
        end
        
        rawFileName = ecg.ECGfilename.Text;
        fieldName = regexprep(rawFileName, '\W', '_');
        if ~isletter(fieldName(1))
            fieldName = ['sig_', fieldName];
        end
        
        % Logika wczytywania sygnałów wg uaktualnionych wytycznych (tryb median)
        if length(rawFileName) >= 4 && strcmp(rawFileName(end-3:end), '.bsm')
            targetFileName = [rawFileName '.medianecg'];
        elseif length(rawFileName) >= 4 && strcmp(rawFileName(end-3:end), '.ecg')
            targetFileName = rawFileName;
        else
            targetFileName = [rawFileName '.ecg'];
        end
        
        if exist(fullfile(basedir, 'ECG_DATA', targetFileName), 'file')
            ECG = loadmat(fullfile(basedir, 'ECG_DATA', targetFileName));
        else
            ECG = [];
        end

        if isfield(ecg,'medianvresults')
            beats = ecg.medianvresults.CyncRESULT;
            ecg = rmfield(ecg,[{'useECGSignal'} {'medianvresults'} {'Attributes'}]);

            % Przypisanie bezpośrednio do fieldName (obok beats)
            DATA.(fieldName).ECG = ECG;
            DATA.(fieldName).filename = targetFileName;

            if ( isfield(ecg,'forpeter') )
                DATA.(fieldName).forpeter = removeText(ecg.forpeter);
            end

            i=0;
            for k=1:length(beats)
                if length(beats)<2
                    selbeat = beats;
                else
                    selbeat = beats{k};
                end
                
                if isfield(selbeat,'initdep') || isfield(selbeat,'finaldep')
                    i=i+1;
                    SPECS = struct();
                    
                    if isfield(selbeat,'Pwave_onset')
                        SPECS.onsetP = removeText(selbeat.Pwave_onset);
                    end
                    if isfield(selbeat,'onsetQRS')
                        SPECS.onsetqrs = removeText(selbeat.onsetQRS);
                    else
                        SPECS.onsetqrs = 1;
                    end
                    if isfield(selbeat,'Twave_end')
                        SPECS.endtwave = removeText(selbeat.Twave_end);
                    else
                        SPECS.endtwave = size(ECG, 2);
                    end
                    if isfield(selbeat,'EndQRS')
                        SPECS.time_Jpoint = removeText(selbeat.EndQRS);
                    else
                        SPECS.time_Jpoint = SPECS.onsetqrs;
                    end
                    if isfield(selbeat,'PeakTwave')
                        SPECS.time_apexT = removeText(selbeat.PeakTwave);
                    end
                    
                    SPECS.qrsduration   = SPECS.time_Jpoint - SPECS.onsetqrs;
                    SPECS.qrstduration  = SPECS.endtwave - SPECS.onsetqrs;

                    DATA.(fieldName).beats{i}.SPECS = SPECS;
                    
                    if isfield(selbeat,'initdep')
                        DATA.(fieldName).beats{i}.initdep = removeText(selbeat.initdep)';
                    end
                    if isfield(selbeat,'initdepvelocity')
                        DATA.(fieldName).beats{i}.initdepvelocity = removeText(selbeat.initdepvelocity);
                    end

                    if isfield(selbeat,'finaldep')
                        DATA.(fieldName).beats{i}.finaldep = removeText(selbeat.finaldep)';                        
                    end
                    if isfield(selbeat,'initialcorrelation')
                        DATA.(fieldName).beats{i}.initialcorrelation = removeText(selbeat.initialcorrelation)';
                    end
                    if isfield(selbeat,'finalcorrelation')
                        DATA.(fieldName).beats{i}.finalcorrelation = removeText(selbeat.finalcorrelation)';
                    end
                    if isfield(selbeat,'finalrep')
                        DATA.(fieldName).beats{i}.finalrep = removeText(selbeat.finalrep)';
                        DATA.(fieldName).beats{i}.rep      = removeText(selbeat.finalrep)';
                    end
                    
                    if ( isfield(selbeat,'massActivated') )
                        DATA.(fieldName).beats{i}.massActivated = removeText(selbeat.massActivated)';
                    end
                    if ( isfield(selbeat,'massLVseptum') )
                        DATA.(fieldName).beats{i}.massLVseptum = removeText(selbeat.massLVseptum)';
                    end
                    if ( isfield(selbeat,'massLVanterior') )
                        DATA.(fieldName).beats{i}.massLVanterior = removeText(selbeat.massLVanterior)';
                    end
                    if ( isfield(selbeat,'massLVposterior') )
                        DATA.(fieldName).beats{i}.massLVposterior = removeText(selbeat.massLVposterior)';
                    end
                    if ( isfield(selbeat,'massRVseptum') )
                        DATA.(fieldName).beats{i}.massRVseptum = removeText(selbeat.massRVseptum)';
                    end
                    if ( isfield(selbeat,'massRVfreewall') )
                        DATA.(fieldName).beats{i}.massRVfreewall = removeText(selbeat.massRVfreewall)';
                    end
                    if ( isfield(selbeat,'massLVActivated') )
                        DATA.(fieldName).beats{i}.massLVActivated = removeText(selbeat.massLVActivated)';
                    end
                    if ( isfield(selbeat,'massRVActivated') )
                        DATA.(fieldName).beats{i}.massRVActivated = removeText(selbeat.massRVActivated)';
                    end
                    if ( isfield(selbeat,'massSEActivated') )
                        DATA.(fieldName).beats{i}.massSEActivated = removeText(selbeat.massSEActivated)';
                    end
                    if isfield(selbeat,'truthorigin')
                        DATA.(fieldName).beats{i}.truthorigin = removeText(selbeat.truthorigin)';
                    elseif isfield(selbeat,'DEPREPtruthLocation')    
                        DATA.(fieldName).beats{i}.truthorigin = removeText(selbeat.DEPREPtruthLocation)';
                    end
                    
                    if ( isfield(selbeat,'beatMeanTSI') )
                        DATA.(fieldName).beats{i}.beatMeanTSI = read3DVertices(selbeat.beatMeanTSI);
                    end
                    if ( isfield(selbeat,'beatMeanTSIHeart') )
                        DATA.(fieldName).beats{i}.beatMeanTSIHeart = read3DVertices(selbeat.beatMeanTSIHeart);
                    end
                    if ( isfield(selbeat,'DEPREPbeatMeanTSI') )
                        DATA.(fieldName).beats{i}.beatMeanTSI = read3DVertices(selbeat.DEPREPbeatMeanTSI);
                    end
                    if ( isfield(selbeat,'DEPREPbeatMeanTSIHeart') )
                        DATA.(fieldName).beats{i}.beatMeanTSIHeart = read3DVertices(selbeat.DEPREPbeatMeanTSIHeart);
                    end
                    if ( isfield(selbeat,'beatVcgHeart') )
                        DATA.(fieldName).beats{i}.beatVcgHeart = read3DVertices(selbeat.beatVcgHeart);
                    end
                    if ( isfield(selbeat,'meanQRSaxis') )
                        DATA.(fieldName).beats{i}.meanQRSaxis = read3DVertices(selbeat.meanQRSaxis);
                    end
                    if ( isfield(selbeat,'meanQRSaxisPosition') )
                        DATA.(fieldName).beats{i}.meanQRSaxisPosition = read3DVertices(selbeat.meanQRSaxisPosition);
                    end
                    if ( isfield(selbeat,'DEPREPmeanQRSaxis') )
                        DATA.(fieldName).beats{i}.meanQRSaxis = read3DVertices(selbeat.DEPREPmeanQRSaxis);
                    end
                    if ( isfield(selbeat,'DEPREPmeanQRSaxisPosition') )
                        DATA.(fieldName).beats{i}.meanQRSaxisPosition = read3DVertices(selbeat.DEPREPmeanQRSaxisPosition);
                    end
                    if ( isfield(selbeat,'initialAngle') )
                        DATA.(fieldName).beats{i}.initialAngle = removeText(selbeat.initialAngle);
                    end
                    
                    if size(ECG,1)>1
                        sIdx = max(1, SPECS.onsetqrs);
                        eIdx = min(SPECS.endtwave, size(ECG,2));
                        if isempty(eIdx) || eIdx < sIdx
                            eIdx = size(ECG,2);
                        end
                        DATA.(fieldName).beats{i}.ECG = ECG(:, sIdx:eIdx);
                    end
                end
            end
        end
    end
else
    disp('no analysed data found')
end