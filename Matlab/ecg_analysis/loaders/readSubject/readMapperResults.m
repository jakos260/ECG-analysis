function DATA = readMapperResults(DATA, subject_folder)
    % Zdefiniowanie ścieżki do folderu mapper
    mapperdir = fullfile(subject_folder, 'mapper');
    
    % Inicjalizacja pustego pola MAPPER
    if ~isfield(DATA, 'MAPPER')
        DATA.MAPPER = struct();
    end

    % Szukamy plików *.imap
    imap_files = dir(fullfile(mapperdir, '*.imap'));
    
    if isempty(imap_files)
        warning('Nie znaleziono plików *.imap w folderze %s', mapperdir);
        return;
    end

    for i = 1:length(imap_files)
        fn = fullfile(imap_files(i).folder, imap_files(i).name);
        
        % Wyciągnięcie nazwy pliku imap (np. '12ECG' lub 'BSM')
        [~, mapper_name, ~] = fileparts(imap_files(i).name);
        
        % Zabezpieczenie nazwy pola (np. '12ECG' -> 'x12ECG')
        if ~isempty(mapper_name) && ~isletter(mapper_name(1))
            mapper_cat = ['x' mapper_name];
        else
            mapper_cat = mapper_name;
        end
        
        if isempty(mapper_cat)
            continue;
        end

        try
            imapStr = fileread(fn);
        catch
            continue;
        end

        % Wyodrębnienie poszczególnych sekcji ECGFiledata
        fileBlocks = regexp(imapStr, '(?s)<ECGFiledata.*?</ECGFiledata>', 'match');
        
        for b = 1:length(fileBlocks)
            block = fileBlocks{b};
            
            % Pobranie oryginalnej nazwy pliku EKG z bloku MAPPERA
            fname_token = regexp(block, '<ECGfilename>\s*([^<]+)\s*</ECGfilename>', 'tokens', 'once');
            if isempty(fname_token)
                continue;
            end
            rawFileName = strtrim(fname_token{1});

            % Rozbicie nazwy i usunięcie głównego rozszerzenia (.ecg / .bsm)
            type_idx = regexpi(rawFileName, '\.(ecg|bsm)');
            if isempty(type_idx)
                continue;
            end

            orig_base = rawFileName(1:type_idx(1)-1);
            rest = rawFileName(type_idx(1):end);

            if startsWith(rest, '.ecg', 'IgnoreCase', true)
                if length(rest) > 4, rest = rest(5:end); else rest = ''; end
            elseif startsWith(rest, '.bsm', 'IgnoreCase', true)
                if length(rest) > 4, rest = rest(5:end); else rest = ''; end
            end

            beat_str = '';
            ext_str = 'ALL';

            if startsWith(rest, '_')
                dot_idx = strfind(rest, '.');
                if isempty(dot_idx)
                    beat_str = rest(2:end);
                else
                    beat_str = rest(2:dot_idx(1)-1);
                    ext_str = upper(rest(dot_idx(1)+1:end));
                end
            elseif startsWith(rest, '.')
                ext_str = upper(rest(2:end));
            end

            % Zmapowanie nazwy używając słownika SIGNALS.INFO
            out_name = '';
            if isfield(DATA, 'SIGNALS') && isfield(DATA.SIGNALS, 'INFO')
                out_names = fieldnames(DATA.SIGNALS.INFO);
                for k = 1:length(out_names)
                    if isfield(DATA.SIGNALS.INFO.(out_names{k}), 'original_name') && ...
                       strcmpi(DATA.SIGNALS.INFO.(out_names{k}).original_name, orig_base)
                        out_name = out_names{k};
                        break;
                    end
                end
            end
            
            if isempty(out_name)
                safe_base = orig_base;
                if ~isletter(safe_base(1)), safe_base = ['s' safe_base]; end
                if isfield(DATA.SIGNALS.INFO, safe_base)
                    out_name = safe_base;
                else
                    continue; % Pomiń, jeśli sygnał nie istnieje
                end
            end

            % --- Wydobycie właściwych wyników (Beats) ---
            beatBlocks = regexp(block, '(?s)<(CyncRESULT|MEANTSIDIAGNOSTIC).*?</\1>', 'match');
            valid_beat_idx = 1; 
            
            for d = 1:length(beatBlocks)
                dBlock = beatBlocks{d};
                
                % Weryfikacja 1: Pomijanie pustych bloków z flagami = -1
                peak_qrs_val = parseTag(dBlock, 'PeakQRS', NaN);
                onset_qrs_val = parseTag(dBlock, 'onsetQRS', NaN);
                
                if (peak_qrs_val == -1) || (onset_qrs_val == -1)
                    continue; 
                end
                
                MAPRES = struct();
                MAPRES.initdep   = parseTagArray(dBlock, 'initdep');
                MAPRES.finaldep  = parseTagArray(dBlock, 'finaldep');
                MAPRES.initrep   = parseTagArray(dBlock, 'initialcorrelation'); 
                MAPRES.finalrep  = parseTagArray(dBlock, 'finalrep');
                MAPRES.mudep     = parseTagArray(dBlock, 'initdepvelocity'); 
                MAPRES.murep     = parseTagArray(dBlock, 'repslope'); 
                MAPRES.initfoci  = parseTagArray(dBlock, 'truthorigin'); 
                
                if isempty(MAPRES.initfoci)
                     MAPRES.initfoci = parseTagArray(dBlock, 'DEPREPtruthLocation');
                end

                % Weryfikacja 2: Odrzucenie bloku, jeśli wszystkie pola są puste
                if isempty(MAPRES.initdep) && isempty(MAPRES.finaldep) && ...
                   isempty(MAPRES.initrep) && isempty(MAPRES.finalrep) && ...
                   isempty(MAPRES.mudep)   && isempty(MAPRES.murep)    && ...
                   isempty(MAPRES.initfoci)
                    continue; % Pomijamy zapis, idziemy do kolejnego bloku
                end

                % Zapisanie metryk do odpowiedniej kategorii mappera
                if ~isempty(beat_str)
                    beat_num = str2double(beat_str);
                    if ~isnan(beat_num)
                        DATA.MAPPER.(mapper_cat).(out_name).BEATS.(ext_str)(beat_num).RES = MAPRES;
                    end
                else
                    DATA.MAPPER.(mapper_cat).(out_name).BEATS.(ext_str)(valid_beat_idx).RES = MAPRES;
                    valid_beat_idx = valid_beat_idx + 1;
                end
            end
        end
    end
end

% --- Funkcje Pomocnicze ---
function val = parseTag(xmlBlock, tagName, defaultVal)
    pattern = ['<' tagName '>\s*([^<]+)\s*</' tagName '>'];
    tok = regexp(xmlBlock, pattern, 'tokens', 'once');
    if ~isempty(tok)
        parts = strsplit(tok{1}, ',');
        num = str2double(strtrim(parts{1}));
        if ~isnan(num)
            val = num;
            return;
        end
    end
    val = defaultVal;
end

function valArray = parseTagArray(xmlBlock, tagName)
    valArray = [];
    pattern = ['<' tagName '>\s*([^<]+)\s*</' tagName '>'];
    tok = regexp(xmlBlock, pattern, 'tokens', 'once');
    if ~isempty(tok)
        parts = strsplit(tok{1}, ',');
        nums = str2double(strtrim(parts));
        if ~all(isnan(nums))
            valArray = nums';
        end
    end
end