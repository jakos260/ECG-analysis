function DATACASES = readAlvaleResults(fn,ECGDIR)

[basedir,name,~] = fileparts(fn);
if ~exist('ECGDIR')
    ECGDIR = fullfile(basedir,'ECG_DATA');
end

XML=xml2struct(fn);
if isfield(XML.ALVALE{2}.ALVALEDATA.ECG_cases,'ECG_CASE')
    DATA = XML.ALVALE{2}.ALVALEDATA.ECG_cases.ECG_CASE.ECGfiles.ECGFiledata;
elseif isfield(XML.ALVALE{2}.ALVALEDATA.ECG_cases,'CIPS_ECG_CASE')
    DATA = XML.ALVALE{2}.ALVALEDATA.ECG_cases.CIPS_ECG_CASE.ECGfiles.ECGFiledata;
else
    error('not a proper case file')
end
    
if 1%length(XML.ALVALE{2}.ALVALEDATA.ECG_cases.ECG_CASE) == 1
        
    if all(size(DATA) == 1)
        D=DATA;
        clear DATA
        DATA{1} = D;
    end
    
    for i=1:length(DATA)
        if exist(fullfile(ECGDIR,DATA{i}.ECGfilename.Text),'file') 
            ECG = loadmat(fullfile(ECGDIR, DATA{i}.ECGfilename.Text ));
        elseif exist(fullfile(ECGDIR,[DATA{i}.ECGfilename.Text '.ecg']),'file')
            ECG = loadmat(fullfile(ECGDIR,[DATA{i}.ECGfilename.Text '.ecg']));
        elseif exist(fullfile(ECGDIR,[DATA{i}.ECGfilename.Text '.peacg']),'file')
            ECG = loadmat(fullfile(ECGDIR,[DATA{i}.ECGfilename.Text '.peacg']));
        end
        DATACASES{i}.ECG = ECG;
        
        DATACASES{i}.ECGfilename = DATA{i}.ECGfilename.Text;
%         disp([DATA{i}.ECGfilename.Text ' ' num2str([i size(DATACASES{i}.ECG)])])
        if isfield(DATA{i},'selectedVentricularBeats')
            if isfield(DATA{i}.selectedVentricularBeats,'MEANTSIDIAGNOSTIC')
                BEATS = DATA{i}.selectedVentricularBeats.MEANTSIDIAGNOSTIC;
            else
                BEATS = DATA{i}.selectedVentricularBeats.CIPS_MEANTSIRESULT;
            end
                
            if length(BEATS) ==1
                if isfield(BEATS,'meanTSI')
                    DATACASES{i}.BEAT{1} = readResult(BEATS,DATACASES{i}.ECG );
                end
            else
                k= 1;
                for j= 1:length(BEATS)
                    if isfield(BEATS{j},'meanTSI')
                        DATACASES{i}.BEAT{k} = readResult(BEATS{j},DATACASES{i}.ECG );
                        k=k+1;
                    end
                end
            end
        end
    end
    
%     ALVALECASE = XML.ALVALE{2}.ALVALEDATA.ECG_cases.CIPS_ECG_CASE;
%     MODELDATA = XML.ALVALE{2}.ALVALEDATA.volumemodel_group;
%     LEADDATA= XML.ALVALE{2}.ALVALEDATA.ECG_cases.leadsystemmodel_group;
%     METADATA = XML.ALVALE{2}.ALVALEDATA;
%     DATACASES.patientid=name;
end

DATA = DATA{1};
clear XML


%%
function val= removeText(str)
val = str2num(struct2array(str));

function MAT = readMatrix(A)
MAT = reshape(A(3:end),A(2), A(1))';

function VER=read3DVertices(str)
V=removeText(str);
VER = reshape(V,3,length(V)/3)';


function DATA = readLeadsData(LEADDATA)

DATA.leadpos = read3DVertices(LEADDATA.leadsystems.geomleadsystemgroup.electrodepositions);

%%
function SPECS = readResult(BEAT,ECG)

SPECS=[];
if isfield(BEAT, 'Pwave_onset' ), SPECS.onsetP             = removeText(BEAT.Pwave_onset); end
if isfield(BEAT, 'onsetQRS' ), SPECS.onsetqrs              = removeText(BEAT.onsetQRS);end
if isfield(BEAT, 'Twave_end' ), SPECS.endtwave             = removeText(BEAT.Twave_end);end
if isfield(BEAT, 'EndQRS' ), SPECS.time_Jpoint             = removeText(BEAT.EndQRS);end
if isfield(BEAT, 'PeakTwave' ), SPECS.time_apexT           = removeText(BEAT.PeakTwave); end
if isfield(BEAT, 'onsetqrs' ), SPECS.qrsduration           = SPECS.time_Jpoint - SPECS.onsetqrs;end
if isfield(BEAT, 'onsetqrs' ), SPECS.qrstduration          = SPECS.endtwave - SPECS.onsetqrs;    end
if isfield(BEAT, 'onsetqrs' ), SPECS.ECGbeat               = ECG(:,SPECS.onsetqrs:SPECS.endtwave);end
if isfield(BEAT, 'onsetqrs' ), SPECS.ECGqrs                = ECG(:,SPECS.onsetqrs:SPECS.time_Jpoint);end
if isfield(BEAT, 'onsetqrs' ), SPECS.ECGtwave              = ECG(:,SPECS.time_Jpoint:SPECS.endtwave);end
SPECS.meanTSI                                              = read3DVertices(BEAT.meanTSI);
SPECS.meanTSIHeart                                         = read3DVertices(BEAT.meanTSIHeart);
SPECS.beatMeanTSI                                          = read3DVertices(BEAT.beatMeanTSI);
if isfield(BEAT, 'beatMeanTSIHeart'),SPECS.beatMeanTSIHeart= read3DVertices(BEAT.beatMeanTSIHeart);end
SPECS.vcg                                                  = read3DVertices(BEAT.vcg);
SPECS.vcgHeart                                             = read3DVertices(BEAT.vcgHeart);
if isfield(BEAT, 'beatVcgHeart'),SPECS.beatVcgHeart        = read3DVertices(BEAT.beatVcgHeart);end
SPECS.meanQRSaxis                                          = read3DVertices(BEAT.meanQRSaxis);
SPECS.meanQRSaxisPosition                                  = read3DVertices(BEAT.meanQRSaxisPosition);
if isfield(BEAT, 'truthorigin'),SPECS.truthorigin          = read3DVertices(BEAT.truthorigin);end
if isfield(BEAT, 'velocityprofile'), SPECS.velocityprofile = removeText(BEAT.velocityprofile);end

SPECS.focusLVratio          = removeText(BEAT.focusLVratio);
SPECS.initialAngle          = removeText(BEAT.initialAngle);
SPECS.normalized_area       = removeText(BEAT.normalized_area);
SPECS.area                  = removeText(BEAT.area);
SPECS.initialAngle          = removeText(BEAT.initialAngle);