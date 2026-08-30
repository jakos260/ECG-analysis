function DATA = readGeomPeacsModelDataset(dirname,subject,varargin)

DATA.subject = subject;

modeldir = fullfile(dirname, subject, 'model');
modeldir = [modeldir filesep];

signaldir = fullfile(dirname, subject, 'signals');
signaldir = [signaldir filesep];

leaddir = fullfile(dirname, subject, 'leads');
leaddir = [leaddir filesep];

if exist([modeldir 'atria.adj2d'], 'file')
    [DATA.GEOM.atria.VER, DATA.GEOM.atria.ITRI] = loadtri([modeldir  'atria.tri']);
    [DATA.ATRIA.geom.VER, DATA.ATRIA.geom.ITRI] = loadtri([modeldir  'atria.tri']);
    DATA.ATRIA.ADJsurf = loadmat([modeldir  'atria.adj2d']);
    DATA.ATRIA.ADJ3D = loadmat([modeldir  'atria.adj3d']);
    DATA.ATRIA.DISTsurf = loadmat([modeldir  'atria.dst2d']);
    DATA.ATRIA.DIST3D = loadmat([modeldir  'ATRIA.dst3d']);
    DATA.ATRIA.ADJANIS = loadmat([modeldir  'ATRIA.adjanis']);
    DATA.ATRIA.DISTANIS = loadmat([modeldir  'ATRIA.dstanis']);
    DATA.GEOM.atria.geomtyp = loadmat([modeldir  'atria.typ']);
end

[DATA.VENTR.geom.VER, DATA.VENTR.geom.ITRI] = loadtri([modeldir  'ventricles.tri']);

DATA.VENTR.ADJsurf          = loadmat([modeldir  'ventricles.adj2d']);
DATA.VENTR.ADJ3D            = loadmat([modeldir  'ventricles.adj3d']);
DATA.VENTR.ADJanis          = loadmat([modeldir  'ventricles.adjanis']);
DATA.VENTR.ADJneigh         = loadmat([modeldir  'ventricles.adjneigh']);
DATA.VENTR.DISTsurf         = loadmat([modeldir  'ventricles.dst2d']);
DATA.VENTR.DIST3D           = loadmat([modeldir  'ventricles.dst3d']);
DATA.VENTR.DISTanis         = loadmat([modeldir  'ventricles.dstanis']);
DATA.VENTR.HEARTDIST        = loadmat([modeldir  'ventricles.heartdist']);
DATA.GEOM.ventr.geomtyp     = loadmat([modeldir  'ventricles.typ']);
DATA.GEOM.ventr.scar        = loadmat([modeldir  'ventricles.scar']);
DATA.GEOM.ventr.walls       = loadmat([modeldir  'ventricles.walls']);
DATA.GEOM.ventr.segments    = loadmat([modeldir  'ventricles.segments']);

[DATA.GEOM.ventr.VER, DATA.GEOM.ventr.ITRI]    = loadtri([modeldir  'ventricles.tri']);
[DATA.GEOM.lcav.VER, DATA.GEOM.lcav.ITRI]      = loadtri([modeldir  'lcav.tri']);
[DATA.GEOM.rcav.VER, DATA.GEOM.rcav.ITRI]      = loadtri([modeldir  'rcav.tri']);

if exist([modeldir  'llung.tri'], 'file')
    [DATA.GEOM.llung.VER,DATA.GEOM.llung.ITRI] = loadtri([modeldir  'llung.tri']);
end
if exist( [modeldir  'rlung.tri'], 'file')
    [DATA.GEOM.rlung.VER,DATA.GEOM.rlung.ITRI] = loadtri([modeldir  'rlung.tri']);
end

[DATA.GEOM.thorax.VER,DATA.GEOM.thorax.ITRI]   = loadtri([modeldir  'thorax.tri']);

if exist([modeldir  'liver.tri'], 'file')
    [DATA.GEOM.liver.VER,DATA.GEOM.liver.ITRI] = loadtri([modeldir  'liver.tri']);
end
if exist([modeldir  'ribcage.tri'], 'file')
    [DATA.GEOM.ribcage.VER,DATA.GEOM.ribcage.ITRI] = loadtri([modeldir  'ribcage.tri']);
end
if exist([modeldir  'fatpad_1.tri'], 'file')
    [DATA.GEOM.fatpad_1.VER,DATA.GEOM.fatpad_1.ITRI] = loadtri([modeldir  'fatpad_1.tri']);
end
if exist([modeldir  'fatpad_2.tri'], 'file')
    [DATA.GEOM.fatpad_2.VER,DATA.GEOM.fatpad_2.ITRI] = loadtri([modeldir  'fatpad_2.tri']);
end
if exist([modeldir  'coil.tri'], 'file')
    [DATA.GEOM.coil.VER,DATA.GEOM.coil.ITRI]   = loadtri([modeldir  'coil.tri']);
end

if nargin ==3 && ischar(varargin{end})
    lead_path = varargin{end};
else
    lead_path = leaddir;
end

if isempty(dir(fullfile(lead_path, '*.lead'))) && isempty(dir(fullfile(lead_path, '*.leads'))) && isempty(dir(fullfile(lead_path, '*lead12*'))) && exist(fullfile(lead_path, subject), 'dir')
    lead_path = fullfile(lead_path, subject);
end

if nargin == 3 && ~ischar(varargin{end})
    wct = varargin{end};
    wct_arg = true;
else
    wct_arg = false;
end
        
if ~wct_arg
    dd = [dir(fullfile(lead_path, '*.lead')); ...
          dir(fullfile(lead_path, '*.leads')); ...
          dir(fullfile(lead_path, '*lead12*'))];
    
    if ~isempty(dd)
        idx_12lead = 1; % default if first file if there is no standard12lead
        for i = 1:length(dd)
            if contains(dd(i).name, 'lead12', 'IgnoreCase', true) || ...
               contains(dd(i).name, 'standard12lead', 'IgnoreCase', true) || ...
               contains(dd(i).name, '12lead', 'IgnoreCase', true)
                idx_12lead = i;
                break;
            end
        end
        
        wct_file = fullfile(dd(idx_12lead).folder, dd(idx_12lead).name);
        
        ver = loadmat(wct_file);
        
        if size(ver, 2) >= 4
            ver = ver(:, 2:4);
        end
        
        if isfield(DATA, 'GEOM') && isfield(DATA.GEOM, 'thorax')
            meanTh = mean(DATA.GEOM.thorax.VER);
            wct = [];
            for i = 1:3
                TRIS = linetris(DATA.GEOM.thorax.VER, DATA.GEOM.thorax.ITRI, ver(i,:), [meanTh(1) meanTh(2) ver(i,3)]);
                if isempty(TRIS)
                    TRIS = linetris(DATA.GEOM.thorax.VER, DATA.GEOM.thorax.ITRI, ver(i,:), meanTh);
                    TRIS(TRIS(:,2)<0, :) = [];
                end
                TRIS(abs(TRIS(:,end)) > min(abs(TRIS(:,end))), :) = [];
                wct = [wct DATA.GEOM.thorax.ITRI(TRIS(1),1)];
            end
        end
    else
        fprintf('[WARNING] No leads files found for WCT.\n');
    end
end

if exist('wct', 'var') && ~isempty(wct)
    if exist([modeldir 'all.aedl'], 'file')
        AA = 40 * loadmat([modeldir 'all.aedl']);
    elseif exist([modeldir 'thorax.aedl'], 'file') 
        AA = 40 * loadmat([modeldir 'thorax.aedl']);
    else
        AA=[];
    end
    
    if exist([modeldir 'all.vedl'], 'file')
        AV = 40*loadmat([modeldir 'all.vedl']);
    elseif exist([modeldir 'thorax.vedl'], 'file')
        AV = 40*loadmat([modeldir 'thorax.vedl']);
    else
        AV=[];
    end
    
    if ~isempty(AV)
        n = length(DATA.GEOM.thorax.VER);
        DATA.VENTR.THORAX = AV(1:n,:);
        AV = doWCT(AV,calcAwct(DATA.VENTR.THORAX,wct));

        DATA.VENTR.THORAX = AV(1:n,:);
        if size(AV,1)>n
            DATA.VENTR.LCAV= (AV(n+1: n + length(DATA.GEOM.lcav.VER),:));
            n = n + length(DATA.GEOM.lcav.VER);
        end
        if size(AV,1)>n
            DATA.VENTR.RCAV= (AV(n+1: n + length(DATA.GEOM.rcav.VER),:));
            n = n + length(DATA.GEOM.rcav.VER);
        end

        if size(AV,1)>n && isfield(DATA.GEOM,'rlung')
            DATA.VENTR.RLUNG= (AV(n+1: n + length(DATA.GEOM.rlung.VER),:));
            n = n + length(DATA.GEOM.rlung.VER);
        end
        if size(AV,1)>n && isfield(DATA.GEOM,'llung')
            DATA.VENTR.LLUNG= (AV(n+1: n + length(DATA.GEOM.llung.VER),:));
            n = n + length(DATA.GEOM.llung.VER);
        end

        if isfield(DATA.GEOM,'liver') 
            DATA.VENTR.LIVER= (AV(n+1: n + length(DATA.GEOM.liver.VER),:));
            n = n + length(DATA.GEOM.liver.VER);
        end

        if isfield(DATA.GEOM,'ribcage') 
            DATA.VENTR.RIBCAGE= (AV(n+1: n + length(DATA.GEOM.ribcage.VER),:));
            n = n + length(DATA.GEOM.ribcage.VER);
        end
        if isfield(DATA.GEOM,'fatpad_1') 

            DATA.VENTR.FATPAD_1= (AV(n+1: n + length(DATA.GEOM.fatpad_1.VER),:));
            n = n + length(DATA.GEOM.fatpad_1.VER);
        end

        if isfield(DATA.GEOM,'fatpad_2') 

            DATA.VENTR.FATPAD_1= (AV(n+1: n + length(DATA.GEOM.fatpad_2.VER),:));
            n = n + length(DATA.GEOM.fatpad_2.VER);
        end
        if isfield(DATA.GEOM,'coil') 
            DATA.VENTR.COIL= (AV(n+1: n + length(DATA.GEOM.coil.VER),:));
            n = n + length(DATA.GEOM.coil.VER);
        end
        
        if size(AV,1)>n
            DATA.VENTR.VENTRICLES = (AV(end-length(DATA.GEOM.ventr.VER)+1:end,:));
        end
    end
    
    if size(AA,1) > n
        n = length(DATA.GEOM.thorax.VER);
        DATA.ATRIA.THORAX = AA(1:n,:);
        AA = doWCT(AA,calcAwct(DATA.ATRIA.THORAX,wct));
        DATA.ATRIA.THORAX = AA(1:n,:);
        DATA.ATRIA.LCAV= (AA(n+1: n + length(DATA.GEOM.lcav.VER),:));
        n = n + length(DATA.GEOM.lcav.VER);

        DATA.ATRIA.RCAV= (AA(n+1: n + length(DATA.GEOM.rcav.VER),:));
        n = n + length(DATA.GEOM.rcav.VER);
        if isfield(DATA.GEOM,'rlung')
            DATA.ATRIA.RLUNG= (AA(n+1: n + length(DATA.GEOM.rlung.VER),:));
            n = n + length(DATA.GEOM.rlung.VER);
        end
        if isfield(DATA.GEOM,'llung')                    
            DATA.ATRIA.LLUNG= (AA(n+1: n + length(DATA.GEOM.llung.VER),:));
            n = n + length(DATA.GEOM.llung.VER);
        end
        if isfield(DATA.GEOM,'liver')
            DATA.VENTR.LIVER= (AV(n+1: n + length(DATA.GEOM.liver.VER),:));
            n = n + length(DATA.GEOM.liver.VER);
        end

        if isfield(DATA.GEOM,'ribcage')
            DATA.ATRIA.RIBCAGE= (AA(n+1: n + length(DATA.GEOM.ribcage.VER),:));
            n = n + length(DATA.GEOM.ribcage.VER);
        end
        if isfield(DATA.GEOM,'fatpad_1')
            DATA.ATRIA.FATPAD_1= (AA(n+1: n + length(DATA.GEOM.fatpad_1.VER),:));
            n = n + length(DATA.GEOM.fatpad_1.VER);
        end
        if isfield(DATA.GEOM,'fatpad_2')
            DATA.ATRIA.FATPAD_2= (AA(n+1: n + length(DATA.GEOM.fatpad_2.VER),:));
            n = n + length(DATA.GEOM.fatpad_2.VER);
        end
        if isfield(DATA.GEOM,'coil') 
            DATA.ATRIA.COIL= (AV(n+1: n + length(DATA.GEOM.coil.VER),:));
            n = n + length(DATA.GEOM.coil.VER);
        end
        DATA.ATRIA.ATRIA = (AA(end-length(DATA.GEOM.atria.VER)+1: end,:));
    end
end

dd_leads = [dir(fullfile(lead_path, '*.lead')); dir(fullfile(lead_path, '*.leads'))];
for i = 1:length(dd_leads)
    filepath = fullfile(dd_leads(i).folder, dd_leads(i).name);
    
    ver = loadmat(filepath);
    if isempty(ver)
        continue;
    end

    if size(ver, 2) >= 4
        ver = ver(:, 2:4);
    end

    name = clean_fieldname(dd_leads(i).name, subject);

    [AMA_A, AMA_V] = getAMA(DATA, ver);
    
    DATA.VENTR.LEADPOS.(name) = ver;
    DATA.VENTR.AMA.(name) = AMA_V;

    if ~isempty(AMA_A) && isfield(DATA, 'ATRIA')
        DATA.ATRIA.LEADPOS.(name) = ver;
        DATA.ATRIA.AMA.(name) = AMA_A;
    end        
end

% load signal info file into a TEMPORARY variable
jsonFileName = 'info.json';
raw_info = struct();
if exist([signaldir, jsonFileName], "file")
    jsonStr = fileread([signaldir, jsonFileName]);
    raw_info = jsondecode(jsonStr);    
end

% Initialize the final INFO structure as empty
DATA.VENTR.SIGNALS.INFO = struct();

% Load signals and categorize by extensions
dd_ecgs = [dir(fullfile(signaldir, '*.bsm*')); dir(fullfile(signaldir, '*.ecg*'))];

% Inicjalizacja liczników i "pamięci" mapowania nazw
s_count = 1;
bsl_count = 1;
name_map = struct();

for i = 1:length(dd_ecgs)
    filepath = fullfile(dd_ecgs(i).folder, dd_ecgs(i).name);
    filename = dd_ecgs(i).name;
    
    if contains(filepath, '.imap') || contains(filepath, '.iecg')
        continue;
    end
    
    sig = loadmat(filepath);
    if isempty(sig)
        continue;
    end

    % Wyciągamy oryginalną nazwę bazową (bez rozszerzeń)
    base_name = filename;
    base_name = strrep(base_name, '.bsm.medianecg', '');
    base_name = strrep(base_name, '.bsm', '');
    base_name = strrep(base_name, '.ecg', '');
    
    % Używamy funkcji clean_fieldname, aby uzyskać znormalizowany klucz (stara nazwa)
    map_key = clean_fieldname(filename, subject);
    
    % Sprawdzamy, czy ten sygnał dostał już nową nazwę w tej iteracji pacjenta
    if isfield(name_map, map_key)
        new_name = name_map.(map_key);
    else
        % Generujemy nową nazwę i zapisujemy ją w mapie
        if startsWith(base_name, 'B', 'IgnoreCase', true) || ...
           contains(base_name, 'Baseline', 'IgnoreCase', true) || ...
           contains(base_name, 'BSL', 'IgnoreCase', true)
            new_name = sprintf('sBSL%d', bsl_count);
            bsl_count = bsl_count + 1;
        else
            new_name = sprintf('s%d', s_count);
            s_count = s_count + 1;
        end
        
        name_map.(map_key) = new_name;
        
        % --- Zapis do struktury INFO (wykonuje się tylko RAZ dla każdego przypadku) ---
        DATA.VENTR.SIGNALS.INFO.(new_name) = struct();
        DATA.VENTR.SIGNALS.INFO.(new_name).name = base_name;
        
        % Szukamy pola 'beats' w tymczasowej strukturze z JSONa
        % UWAGA: jsondecode czasami zachowuje oryginalne nazwy literowe bez dodawania 'x',
        % więc bezpieczniej jest sprawdzić dwie wersje klucza (z 'x' i bez 'x').
        matched_key = '';
        if isfield(raw_info, map_key)
            matched_key = map_key;
        elseif isfield(raw_info, map_key(2:end))
            matched_key = map_key(2:end);
        end
        
        if ~isempty(matched_key) && isfield(raw_info.(matched_key), 'beats')
            DATA.VENTR.SIGNALS.INFO.(new_name).beats = raw_info.(matched_key).beats;
        end
    end
    
    % Route data to specific substructures based on extension
    if contains(filename, '.bsm.medianecg', 'IgnoreCase', true)
        DATA.VENTR.SIGNALS.BSMmedECG.(new_name) = removeBaseline(sig);
    elseif contains(filename, '.bsm', 'IgnoreCase', true)
        DATA.VENTR.SIGNALS.BSM.(new_name) = removeBaseline(sig);
    elseif contains(filename, '.ecg', 'IgnoreCase', true)
        DATA.VENTR.SIGNALS.ECG.(new_name) = removeBaseline(sig);
    end
end


% CLASSIFICATION LOGIC
if isfield(DATA, 'VENTR') && isfield(DATA.GEOM, 'ventr') && isfield(DATA.GEOM.ventr, 'VER')
    DATA.GEOM.ventr.endoVER=zeros(1,length(DATA.GEOM.ventr.VER));
    for i=1:length(DATA.GEOM.ventr.VER)
        if any(DATA.GEOM.lcav.VER(:,1)==DATA.GEOM.ventr.VER(i,1) & ...
                DATA.GEOM.lcav.VER(:,2)==DATA.GEOM.ventr.VER(i,2) & ...
                DATA.GEOM.lcav.VER(:,3)==DATA.GEOM.ventr.VER(i,3))
            DATA.GEOM.ventr.endoVER(i) = 2;
        end
    end

    for i=1:length(DATA.GEOM.ventr.VER)
        if any(DATA.GEOM.rcav.VER(:,1)==DATA.GEOM.ventr.VER(i,1) & ...
                DATA.GEOM.rcav.VER(:,2)==DATA.GEOM.ventr.VER(i,2) & ...
                DATA.GEOM.rcav.VER(:,3)==DATA.GEOM.ventr.VER(i,3))
            DATA.GEOM.ventr.endoVER(i) = 1;
        end
    end

    DATA.GEOM.ventr.type = DATA.GEOM.ventr.endoVER;

    if isfield(DATA.VENTR, 'ADJ3D')
        ADJ = DATA.VENTR.ADJ3D;
        if size(ADJ,1) == length(DATA.GEOM.ventr.VER)
            ADJ(DATA.VENTR.ADJsurf>0)=0;
            for i=1:length(DATA.GEOM.ventr.VER)
                lvD = min(ADJ(i,DATA.GEOM.ventr.endoVER==2 & ADJ(i,:)>0 & ADJ(i,:)< 35));
                rvD = min(ADJ(i,DATA.GEOM.ventr.endoVER==1 & ADJ(i,:)>0 & ADJ(i,:)< 35));
                evD = min(ADJ(i,DATA.GEOM.ventr.endoVER==0 & ADJ(i,:)>0 & ADJ(i,:)< 35));
                if isempty(evD) && DATA.GEOM.ventr.endoVER(i) ~= 0
                    if DATA.GEOM.ventr.endoVER(i) ==1
                        DATA.GEOM.ventr.type(i) = 3;% rvseptum endo
                    else
                        DATA.GEOM.ventr.type(i) = 4;% lvseptum endo
                    end
                elseif isempty(lvD) && DATA.GEOM.ventr.endoVER(i) ~= 2
                    if DATA.GEOM.ventr.endoVER(i) == 0
                        DATA.GEOM.ventr.type(i) = 5;% rvfree epi
                    else
                        DATA.GEOM.ventr.type(i) = 6;% rvfree endo
                    end
                elseif isempty(rvD) && DATA.GEOM.ventr.endoVER(i) ~= 1
                    if DATA.GEOM.ventr.endoVER(i) == 0
                        DATA.GEOM.ventr.type(i) = 7;% lvfree epi
                    else
                        DATA.GEOM.ventr.type(i) = 8;% lvfree endo
                    end
                else
                    if lvD < rvD * 3 % this is lv
                        if DATA.GEOM.ventr.endoVER(i) == 0 %epi
                            DATA.GEOM.ventr.type(i) = 7;% lvfree epi
                        elseif DATA.GEOM.ventr.endoVER(i) == 1 %rv endo
                            DATA.GEOM.ventr.type(i) = 3;% rvseptum endo
                        else
                            DATA.GEOM.ventr.type(i) = 4;% lvseptum endo
                        end
                    else
                        if DATA.GEOM.ventr.endoVER(i) == 0 %epi
                            DATA.GEOM.ventr.type(i) = 5;% rvfree epi
                        elseif DATA.GEOM.ventr.endoVER(i) == 1 %rv endo
                            if isempty(lvD)
                                DATA.GEOM.ventr.type(i) = 6;% rvfree endo
                            else
                                DATA.GEOM.ventr.type(i) = 3;% rvfree endo
                            end
                        else
                            if isempty(rvD)
                                DATA.GEOM.ventr.type(i) = 8;% lvfree endo
                            else
                                DATA.GEOM.ventr.type(i) = 4;% lvseptum
                            end
                        end
                    end
                end
            end
            
        end
    end
end

if isfield(DATA.GEOM,'atria') && isfield(DATA.GEOM.atria, 'VER')
    DATA.GEOM.atria.endoVER=zeros(1,length(DATA.GEOM.atria.VER));
    for i=1:length(DATA.GEOM.atria.VER)
        if any(DATA.GEOM.lcav.VER(:,1)==DATA.GEOM.atria.VER(i,1) & ...
                DATA.GEOM.lcav.VER(:,2)==DATA.GEOM.atria.VER(i,2) & ...
                DATA.GEOM.lcav.VER(:,3)==DATA.GEOM.atria.VER(i,3))
            DATA.GEOM.atria.endoVER(i) = 2;
        end
    end
    for i=1:length(DATA.GEOM.atria.VER)
        if any(DATA.GEOM.rcav.VER(:,1)==DATA.GEOM.atria.VER(i,1) & ...
                DATA.GEOM.rcav.VER(:,2)==DATA.GEOM.atria.VER(i,2) & ...
                DATA.GEOM.rcav.VER(:,3)==DATA.GEOM.atria.VER(i,3))
            DATA.GEOM.atria.endoVER(i) = 1;
        end
    end
end

end

%%=======================================================================
%  HELPERS
%%=======================================================================

function s = removeBaseline(signal)
    s = baselinecor(signal);
end

function name = clean_fieldname(filename, subject)
    name = filename;
    name = strrep(name, '.leads', '');
    name = strrep(name, '.lead', '');
    name = strrep(name, '.medianecg', '');
    name = strrep(name, '.bsm', '');
    name = strrep(name, '.ecg', '');
    name = strrep(name, ' ', '');
    name = strrep(name, subject, '');
    name = strrep(name, '-', '_');
    name = strrep(name, '.', '_');
    
    if isempty(name)
        name = 'Signal';
    end
    
    % Prepends 'x' to keep the field name legal in MATLAB
    name = ['x' name];
end

function A = doWCT(Ain,Awct)
if ~isempty(Ain)
    A = Ain - ones(size(Ain,1),1) * Awct;
else
    A= Ain;
end
end

function Awct = calcAwct(AthorsoIn,wct)
Awct = mean(AthorsoIn(wct,:));
end

%%=======================================================================
function [M, extraresult]=loadmat(name)
f=fopen(name);
if (f==-1)
    M=[];
    extraresult='';
    return;
end
[N,nr]=fscanf(f,'%d',2);
if (nr~=2)
    fclose(f);
    f=fopen(name);
    [magic ,nr]=fread(f,8,'char');
    if (char(magic')==';;mbfmat')
        fread(f,1,'char');
        hs=fread(f,1,'long');
        fread(f,1,'char');
        fread(f,1,'char');
        fread(f,1,'char');
        N=fread(f,2,'long');
        M=fread(f,[N(2),N(1)],'double');
    else
        fclose(f);
        f=fopen(name);
        N=fread(f,2,'long');
        M=fread(f,[N(2),N(1)],'float');
    end
else
    M=fscanf(f,'%f',[N(2) N(1)]);
end
[extra,nextra]=fread(f,1000,'char');
fclose(f);
extra = char(extra);
if ~all(isspace(extra))
else
    extra=[];
end
M=M';
extraresult=extra;
end

%%========================================================================
function [pnt, dhk] = loadtri(fn)
fid = fopen(fn, 'rt');
if fid~=-1
    Npnt = fscanf(fid, '%d', 1);
    pnt  = fscanf(fid, '%f', [4, Npnt]);
    pnt  = pnt(2:4,:)';
    if (~(feof(fid)))
        [Ndhk, count] = fscanf(fid, '%d', 1);
        if (count ~= 0)
            dhk = fscanf(fid, '%d', [4, Ndhk]);
            dhk = dhk(2:4,:)';
        end
    else
        dhk = [];
    end
    fclose(fid);
else
    error(['unable to open file: ' fn]);
end
end

%%========================================================================
function [AMA_A,AMA_V] = getAMA(DATA,ver)
    if ~isfield(DATA, 'GEOM') || ~isfield(DATA.GEOM, 'thorax') || ~isfield(DATA.GEOM.thorax, 'VER')
        AMA_A = []; AMA_V = []; return;
    end
    
    meanTh = mean(DATA.GEOM.thorax.VER);
    if isfield(DATA,'ATRIA') && isfield(DATA.ATRIA,'THORAX')
        AMA_A = zeros(size(ver,1),size(DATA.ATRIA.THORAX,2));
    else
        AMA_A = [];
    end
    
    if isfield(DATA, 'VENTR') && isfield(DATA.VENTR, 'THORAX')
        AMA_V = zeros(size(ver,1),size(DATA.VENTR.THORAX,2));
        for i=1:size(ver,1)
            TRIS= linetris(DATA.GEOM.thorax.VER,DATA.GEOM.thorax.ITRI,ver(i,:),[meanTh(1) meanTh(2) ver(i,3)]);
            TRIS(abs(TRIS(:,end))> min(abs(TRIS(:,end))),:)=[];
            itri = DATA.GEOM.thorax.ITRI(TRIS(1),:);
            if isfield(DATA,'ATRIA') && isfield(DATA.ATRIA,'THORAX')
                AMA_A(i,:) = DATA.ATRIA.THORAX(itri(1),:) + ...
                    (DATA.ATRIA.THORAX(itri(2),:) - DATA.ATRIA.THORAX(itri(1),:)) * TRIS(1,3) + ...
                    (DATA.ATRIA.THORAX(itri(3),:) - DATA.ATRIA.THORAX(itri(1),:)) * TRIS(1,4);
            end
            AMA_V(i,:) = DATA.VENTR.THORAX(itri(1),:) + ...
                (DATA.VENTR.THORAX(itri(2),:) - DATA.VENTR.THORAX(itri(1),:)) * TRIS(1,3) + ...
                (DATA.VENTR.THORAX(itri(3),:) - DATA.VENTR.THORAX(itri(1),:)) * TRIS(1,4);
        end
    else
        AMA_V = [];
    end
end