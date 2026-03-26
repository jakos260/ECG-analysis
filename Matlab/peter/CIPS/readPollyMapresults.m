function DATACASES = readPollyMapresults( varargin )

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

DATA=[];
if isfield(POLLYCASE,'ECGfiles')
    DATA = cell(length(POLLYCASE.ECGfiles.ECGFiledata));
    for iEcg=1:length(POLLYCASE.ECGfiles.ECGFiledata)
        if length(POLLYCASE.ECGfiles.ECGFiledata) < 2
            ecg=POLLYCASE.ECGfiles.ECGFiledata;
        else
            ecg=POLLYCASE.ECGfiles.ECGFiledata{iEcg};
        end
        if strcmp(ecg.ECGfilename.Text(end-3:end),'.ecg')
            if exist( fullfile(basedir,'ECG_DATA',[ecg.ECGfilename.Text]),'file')
                ECG = loadmat(fullfile(basedir,'ECG_DATA',[ecg.ECGfilename.Text]));
            else
                ECG=[];
            end
        else
            ECG = loadmat(fullfile(basedir,'ECG_DATA',[ecg.ECGfilename.Text '.ecg']));
        end
        
        DATA{iEcg}.ecg.ECG = ECG;
        DATA{iEcg}.ecg.filename = ecg.ECGfilename.Text;
        
        if isfield(ecg,'autointerpret')
            DATA{iEcg}.ecg.autointerpret = ecg.autointerpret.Text;
%             DATA{iEcg}.ecg.autointerpret = base64decode(ecg.autointerpret.Text);
        end
        
        
        if isfield(ecg,'selectedVentricularBeats')
            beats = ecg.selectedVentricularBeats.CyncRESULT;
            ecg =rmfield(ecg,[{'useECGSignal'} {'selectedVentricularBeats'} {'Attributes'}]);
            
            
            i=0;
            for k=1:length(beats)
                if length(beats) < 2
                    selbeat = beats;
                else
                    selbeat = beats{k};
                end
                if isfield(selbeat,'initdep') %&& isfield(beats{k},'finalrep')
                    i=i+1;
                    if isfield(selbeat,'onsetP')
                        SPECS.onsetP        = removeText(selbeat.Pwave_onset);
                    end
                    SPECS.onsetqrs      = removeText(selbeat.onsetQRS);
                    SPECS.endtwave      = removeText(selbeat.Twave_end);
                    SPECS.time_Jpoint   = removeText(selbeat.EndQRS);
                    SPECS.time_apexT    = removeText(selbeat.PeakTwave);
                    SPECS.qrsduration   = SPECS.time_Jpoint - SPECS.onsetqrs;
                    SPECS.qrstduration  = SPECS.endtwave - SPECS.onsetqrs;
                    
                    DATA{iEcg}.beats{i}.SPECS           = SPECS;
                    DATA{iEcg}.beats{i}.initdep         = removeText(selbeat.initdep)';
                    DATA{iEcg}.beats{i}.initdepvelocity = removeText(selbeat.initdepvelocity);
                    
                    if isfield(selbeat,'finaldep')
                        DATA{iEcg}.beats{i}.finaldep    = removeText(selbeat.finaldep)';
                    end
                    if isfield(selbeat,'initialcorrelation')
                        DATA{iEcg}.beats{i}.initialcorrelation    = removeText(selbeat.initialcorrelation)';
                    end
                    if isfield(selbeat,'finalcorrelation')
                        DATA{iEcg}.beats{i}.finalcorrelation    = removeText(selbeat.finalcorrelation)';
                    end
                    
                    if isfield(selbeat,'finalrep')
                        DATA{iEcg}.beats{i}.rep         = removeText(selbeat.finalrep)';
                    end
                    if ( isfield(selbeat,'massActivated') )
                        DATA{iEcg}.beats{i}.massActivated = removeText(selbeat.massActivated)';
                    end
                    
                    if ( isfield(selbeat,'massLVseptum') )
                        DATA{iEcg}.beats{i}.massLVseptum = removeText(selbeat.massLVseptum)';
                    end
                    if ( isfield(selbeat,'massLVanterior') )
                        DATA{iEcg}.beats{i}.massLVanterior = removeText(selbeat.massLVanterior)';
                    end
                    if ( isfield(selbeat,'massLVposterior') )
                        DATA{iEcg}.beats{i}.massLVposterior = removeText(selbeat.massLVposterior)';
                    end
                    if ( isfield(selbeat,'massRVseptum') )
                        DATA{iEcg}.beats{i}.massRVseptum = removeText(selbeat.massRVseptum)';
                    end
                    if ( isfield(selbeat,'massActivated') )
                        DATA{iEcg}.beats{i}.massRVfreewall = removeText(selbeat.massRVfreewall)';
                    end
                   
                    if ( isfield(selbeat,'massLVActivated') )
                        DATA{iEcg}.beats{i}.massLVActivated = removeText(selbeat.massLVActivated)';
                    end
                    if ( isfield(selbeat,'massRVActivated') )
                        DATA{iEcg}.beats{i}.massRVActivated = removeText(selbeat.massRVActivated)';
                    end
                    if ( isfield(selbeat,'massSEActivated') )
                        DATA{iEcg}.beats{i}.massSEActivated = removeText(selbeat.massSEActivated)';
                    end
                    if isfield(selbeat,'truthorigin')
                        DATA{iEcg}.beats{i}.truthorigin = removeText(selbeat.truthorigin)';
                    elseif isfield(selbeat,'DEPREPtruthLocation') 
                        DATA{iEcg}.beats{i}.truthorigin = removeText(selbeat.DEPREPtruthLocation)';
                    end

                    if ( isfield(selbeat,'DEPREPbeatMeanTSI') )
                        DATA{iEcg}.beats{i}.beatMeanTSI         = read3DVertices(selbeat.DEPREPbeatMeanTSI);
                    elseif ( isfield(selbeat,'beatMeanTSI') )
                        DATA{iEcg}.beats{i}.beatMeanTSI         = read3DVertices(selbeat.beatMeanTSI);
                    end
                    if ( isfield(selbeat,'DEPREPbeatMeanTSIHeart') )
                        DATA{iEcg}.beats{i}.beatMeanTSIHeart         = read3DVertices(selbeat.DEPREPbeatMeanTSIHeart);
                    elseif ( isfield(selbeat,'beatMeanTSIHeart') )
                        DATA{iEcg}.beats{i}.beatMeanTSIHeart         = read3DVertices(selbeat.beatMeanTSIHeart);
                    end
                    if ( isfield(selbeat,'DEPREPbeatVcgHeart') )
                        DATA{iEcg}.beats{i}.beatVcgHeart         = read3DVertices(selbeat.DEPREPbeatVcgHeart);
                    elseif ( isfield(selbeat,'beatVcgHeart') )
                        DATA{iEcg}.beats{i}.beatVcgHeart        = read3DVertices(selbeat.beatVcgHeart);
                    end
                    if ( isfield(selbeat,'DEPREPmeanQRSaxisPosition') )
                        DATA{iEcg}.beats{i}.meanQRSaxis         = read3DVertices(selbeat.DEPREPmeanQRSaxisPosition);
                    elseif ( isfield(selbeat,'meanQRSaxis') )
                        DATA{iEcg}.beats{i}.meanQRSaxis         = read3DVertices(selbeat.meanQRSaxis);
                    end
                    if ( isfield(selbeat,'DEPREPmeanQRSaxisPosition') )
                        DATA{iEcg}.beats{i}.meanQRSaxisPosition = read3DVertices(selbeat.DEPREPmeanQRSaxisPosition);
                    elseif ( isfield(selbeat,'meanQRSaxisPosition') )
                        DATA{iEcg}.beats{i}.meanQRSaxisPosition = read3DVertices(selbeat.meanQRSaxisPosition);
                    end
                    if ~isempty(ECG)
                        DATA{iEcg}.beats{i}.ECG                 = ECG(:,max(1,SPECS.onsetqrs):min(SPECS.endtwave,size(ECG,2)));
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
                    
                    SPECS.onsetP        = removeText(selbeat.Pwave_onset);
                    SPECS.onsetqrs      = removeText(selbeat.onsetQRS);
                    SPECS.endtwave      = removeText(selbeat.Twave_end);
                    SPECS.time_Jpoint   = removeText(selbeat.EndQRS);
                    SPECS.time_apexT    = removeText(selbeat.PeakTwave);
                    SPECS.time_apexU    = -1;
                    SPECS.depSlope      = 2;
                    SPECS.initialSlope  = removeText(selbeat.initialrepslope);
                    
                    SPECS.plateauslope  = removeText(selbeat.platslope);
                    if SPECS.plateauslope == 0
                        SPECS.plateauslope = 0.014;
                    end
                    SPECS.repslope      = removeText(selbeat.repslope);
                    if SPECS.repslope == 0
                        SPECS.repslope = 0.045;
                    end
                    SPECS.repCorrection = 0;
                    SPECS.useCumsum     = 0;
                    SPECS.qrsduration   = SPECS.time_Jpoint - SPECS.onsetqrs;
                    SPECS.qrstduration  = SPECS.endtwave - SPECS.onsetqrs;
                    
                    DATA{iEcg}.beats{i}.dep = removeText(beats{k}.finaldep)';
                    DATA{iEcg}.beats{i}.rep = removeText(beats{k}.finalrep)';
                    DATA{iEcg}.beats{i}.SPECS        = SPECS;
                    DATA{iEcg}.beats{i}.ECG         = ECG(:,max(1,SPECS.onsetP):min(SPECS.endtwave,size(ECG,2)));
                end
            end
        end
    end
else
    disp('no analysed data found')
end

%% selected beats
function DATA = readMedianData(basedir,POLLYCASE)

DATA=[];
iEcg=0;
if isfield(POLLYCASE,'ECGfiles')
    for iE=1:length(POLLYCASE.ECGfiles.ECGFiledata)
        if  length(POLLYCASE.ECGfiles.ECGFiledata) < 2
            ecg=POLLYCASE.ECGfiles.ECGFiledata;
        else
            ecg=POLLYCASE.ECGfiles.ECGFiledata{iE};
        end
        if exist( fullfile(basedir,'ECG_DATA',[ecg.ECGfilename.Text '.medianecg']),'file')
            if exist( fullfile(basedir,'ECG_DATA',[ecg.ECGfilename.Text '.medianecg']),'file')
                ECG = loadmat(fullfile(basedir,'ECG_DATA',[ecg.ECGfilename.Text '.medianecg']));
            else
                ECG=[];
            end
        else
            ECG = loadmat(fullfile(basedir,'ECG_DATA',[ecg.ECGfilename.Text '.ecg.medianecg']));
        end


        if isfield(ecg,'medianvresults')
            beats = ecg.medianvresults.CyncRESULT;
            ecg =rmfield(ecg,[{'useECGSignal'} {'medianvresults'} {'Attributes'}]);
            iEcg = iEcg + 1;

            DATA{iEcg}.ecg.ECG = ECG;
            DATA{iEcg}.ecg.filename = [ecg.ECGfilename.Text '.medianecg'];

            if ( isfield(ecg,'forpeter') )
                DATA{iEcg}.ecg.forpeter        = removeText(ecg.forpeter);
            end

            i=0;
            for k=1:length(beats)
                if length(beats)<2
                    selbeat = beats;
                else
                    selbeat = beats{k};
                end
                if isfield(selbeat,'initdep') %&& isfield(beats{k},'finalrep')
                    i=i+1;
                    if isfield(selbeat,'onsetP')
                        SPECS.onsetP        = removeText(selbeat.Pwave_onset);
                    end
                    SPECS.onsetqrs      = removeText(selbeat.onsetQRS);
                    SPECS.endtwave      = removeText(selbeat.Twave_end);
                    SPECS.time_Jpoint   = removeText(selbeat.EndQRS);
                    SPECS.time_apexT    = removeText(selbeat.PeakTwave);
                    SPECS.qrsduration   = SPECS.time_Jpoint - SPECS.onsetqrs;
                    SPECS.qrstduration  = SPECS.endtwave - SPECS.onsetqrs;

                    DATA{iEcg}.beats{i}.SPECS           = SPECS;
                    DATA{iEcg}.beats{i}.initdep         = removeText(selbeat.initdep)';
                    DATA{iEcg}.beats{i}.initdepvelocity = removeText(selbeat.initdepvelocity);

                    if isfield(selbeat,'finaldep')
                        DATA{iEcg}.beats{i}.finaldep    = removeText(selbeat.finaldep)';                        
                    end
                    if isfield(selbeat,'initialcorrelation')
                        DATA{iEcg}.beats{i}.initialcorrelation    = removeText(selbeat.initialcorrelation)';
                    end
                    if isfield(selbeat,'finalcorrelation')
                        DATA{iEcg}.beats{i}.finalcorrelation    = removeText(selbeat.finalcorrelation)';
                    end
                    if isfield(selbeat,'finalrep')
                        DATA{iEcg}.beats{i}.rep         = removeText(selbeat.finalrep)';
                    end
                    if ( isfield(selbeat,'massActivated') )
                        DATA{iEcg}.beats{i}.massActivated = removeText(selbeat.massActivated)';
                    end
                    
                    if ( isfield(selbeat,'massLVseptum') )
                        DATA{iEcg}.beats{i}.massLVseptum = removeText(selbeat.massLVseptum)';
                    end
                    if ( isfield(selbeat,'massLVanterior') )
                        DATA{iEcg}.beats{i}.massLVanterior = removeText(selbeat.massLVanterior)';
                    end
                    if ( isfield(selbeat,'massLVposterior') )
                        DATA{iEcg}.beats{i}.massLVposterior = removeText(selbeat.massLVposterior)';
                    end
                    if ( isfield(selbeat,'massRVseptum') )
                        DATA{iEcg}.beats{i}.massRVseptum = removeText(selbeat.massRVseptum)';
                    end
                    if ( isfield(selbeat,'massActivated') )
                        DATA{iEcg}.beats{i}.massRVfreewall = removeText(selbeat.massRVfreewall)';
                    end

                    
                    if ( isfield(selbeat,'massLVActivated') )
                        DATA{iEcg}.beats{i}.massLVActivated = removeText(selbeat.massLVActivated)';
                    end
                    if ( isfield(selbeat,'massRVActivated') )
                        DATA{iEcg}.beats{i}.massRVActivated = removeText(selbeat.massRVActivated)';
                    end
                    if ( isfield(selbeat,'massSEActivated') )
                        DATA{iEcg}.beats{i}.massSEActivated = removeText(selbeat.massSEActivated)';
                    end
                    if isfield(selbeat,'truthorigin')
                        DATA{iEcg}.beats{i}.truthorigin = removeText(selbeat.truthorigin)';
                    elseif isfield(selbeat,'DEPREPtruthLocation')    
                        DATA{iEcg}.beats{i}.truthorigin = removeText(selbeat.DEPREPtruthLocation)';
                    end
                    if ( isfield(selbeat,'beatMeanTSI') )
                        DATA{iEcg}.beats{i}.beatMeanTSI         = read3DVertices(selbeat.beatMeanTSI);
                    end
                    if ( isfield(selbeat,'beatMeanTSIHeart') )
                        DATA{iEcg}.beats{i}.beatMeanTSIHeart    = read3DVertices(selbeat.beatMeanTSIHeart);
                    end
                    if ( isfield(selbeat,'DEPREPbeatMeanTSI') )
                        DATA{iEcg}.beats{i}.beatMeanTSI         = read3DVertices(selbeat.DEPREPbeatMeanTSI);
                    end
                    if ( isfield(selbeat,'DEPREPbeatMeanTSIHeart') )
                        DATA{iEcg}.beats{i}.beatMeanTSIHeart    = read3DVertices(selbeat.DEPREPbeatMeanTSIHeart);
                    end

                    
                    if ( isfield(selbeat,'beatMeanTSIHeart') )
                        DATA{iEcg}.beats{i}.beatVcgHeart        = read3DVertices(selbeat.beatVcgHeart);
                    end
                    if ( isfield(selbeat,'meanQRSaxis') )
                        DATA{iEcg}.beats{i}.meanQRSaxis         = read3DVertices(selbeat.meanQRSaxis);
                    end
                    if ( isfield(selbeat,'meanQRSaxisPosition') )
                        DATA{iEcg}.beats{i}.meanQRSaxisPosition = read3DVertices(selbeat.meanQRSaxisPosition);
                    end
                    if ( isfield(selbeat,'DEPREPmeanQRSaxis') )
                        DATA{iEcg}.beats{i}.meanQRSaxis         = read3DVertices(selbeat.DEPREPmeanQRSaxis);
                    end
                    if ( isfield(selbeat,'DEPREPmeanQRSaxisPosition') )
                        DATA{iEcg}.beats{i}.meanQRSaxisPosition = read3DVertices(selbeat.DEPREPmeanQRSaxisPosition);
                    end
                    
                    if ( isfield(selbeat,'initialAngle') )
                        DATA{iEcg}.beats{i}.initialAngle        = removeText(selbeat.initialAngle);
                    end
                    if size(ECG,1)>1
                        DATA{iEcg}.beats{i}.ECG             = ECG(:,max(1,SPECS.onsetqrs):min(SPECS.endtwave,size(ECG,2)));
                    end
                end
            end
        end
    end
else
    disp('no analysed data found')
end