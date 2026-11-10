function SIGNALS = readIecgMorphology(SIGNALS, signaldir)
    % Szukamy plików .iecg w folderze sygnałów i jego głównych podfolderach
    iecg_files = [dir(fullfile(signaldir, '*.iecg')); ...
                  dir(fullfile(signaldir, 'ECG_DATA', '*.iecg')); ...
                  dir(fullfile(signaldir, 'ECG_DATA', '12ECG', '*.iecg')); ...
                  dir(fullfile(signaldir, 'ECG_DATA', 'BSM', '*.iecg'))];

    % Odfiltrowanie ewentualnych folderów i duplikatów
    if ~isempty(iecg_files)
        iecg_files = iecg_files(~[iecg_files.isdir]);
        [~, unique_idx] = unique(cellfun(@(folder, name) fullfile(folder, name), {iecg_files.folder}, {iecg_files.name}, 'UniformOutput', false));
        iecg_files = iecg_files(unique_idx);
    end

    for i = 1:length(iecg_files)
        fn = fullfile(iecg_files(i).folder, iecg_files(i).name);
        
        % Wczytanie pliku jako surowy tekst (omija problem ekstremalnie wolnego xml2struct)
        try
            iecgStr = fileread(fn);
        catch
            continue;
        end

        % Wyodrębnienie poszczególnych sekcji ECGFiledata za pomocą szybkiego regexa
        fileBlocks = regexp(iecgStr, '(?s)<ECGFiledata.*?</ECGFiledata>', 'match');
        
        for b = 1:length(fileBlocks)
            block = fileBlocks{b};
            
            % Pobranie oryginalnej nazwy pliku
            fname_token = regexp(block, '<ECGfilename>\s*([^<]+)\s*</ECGfilename>', 'tokens', 'once');
            if isempty(fname_token)
                continue;
            end
            rawFileName = strtrim(fname_token{1});

            % Rozbicie nazwy (analogicznie do logiki ładowania sygnałów)
            type_idx = regexpi(rawFileName, '\.(ecg|bsm)');
            if isempty(type_idx)
                continue;
            end

            orig_base = rawFileName(1:type_idx(1)-1);
            rest = rawFileName(type_idx(1):end);

            % Bezpieczne ucinanie głównego rozszerzenia
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

            % Poszukiwanie nowej (zmapowanej) nazwy sygnału w słowniku INFO
            out_name = '';
            if isfield(SIGNALS, 'INFO')
                out_names = fieldnames(SIGNALS.INFO);
                for k = 1:length(out_names)
                    if isfield(SIGNALS.INFO.(out_names{k}), 'original_name') && ...
                       strcmpi(SIGNALS.INFO.(out_names{k}).original_name, orig_base)
                        out_name = out_names{k};
                        break;
                    end
                end
            end

            % Fallback: jeżeli plik nie był na liście mapowania, a sygnał wczytano pod oryginalną nazwą
            if isempty(out_name)
                safe_base = orig_base;
                if ~isletter(safe_base(1)), safe_base = ['s' safe_base]; end
                if isfield(SIGNALS.INFO, safe_base)
                    out_name = safe_base;
                else
                    continue; % Sygnał nie istnieje w strukturze, pomijamy metryki
                end
            end

            % --- Wydobycie wartości morfologicznych ---
            SPECS = struct();
            
            % Wydobycie wszystkich bloków MEANTSIDIAGNOSTIC wewnątrz danego ECGFiledata
            diagnosticBlocks = regexp(block, '(?s)<MEANTSIDIAGNOSTIC.*?</MEANTSIDIAGNOSTIC>', 'match');
            validBlock = '';
            
            % Szukamy bloku, w którym PeakQRS NIE WYNOSI -1
            for d = 1:length(diagnosticBlocks)
                dBlock = diagnosticBlocks{d};
                peak_qrs_val = parseTag(dBlock, 'PeakQRS', NaN);
                
                if ~isnan(peak_qrs_val) && peak_qrs_val ~= -1
                    validBlock = dBlock;
                    break; % Znaleziono poprawny blok
                end
            end
            
            % Jeśli żaden poprawny blok nie został znaleziony (lub w ogóle ich nie ma), używamy całego blocku jako fallback
            if isempty(validBlock)
                validBlock = block;
            end
            
            % Funkcja przeszukująca ZWERYFIKOWANY blok tekstu
            getVal = @(tag, defaultVal) parseTag(validBlock, tag, defaultVal);
            
            SPECS.onsetP      = getVal('Pwave_onset', NaN);
            SPECS.onsetqrs    = getVal('onsetQRS', 1);
            SPECS.endtwave    = getVal('Twave_end', NaN);
            SPECS.time_Jpoint = getVal('EndQRS', SPECS.onsetqrs);
            SPECS.time_apexT  = getVal('PeakTwave', NaN);
            SPECS.PeakQRS     = getVal('PeakQRS', NaN);

            % Obliczenia czasów trwania
            SPECS.qrsduration  = SPECS.time_Jpoint - SPECS.onsetqrs;
            SPECS.qrstduration = SPECS.endtwave - SPECS.onsetqrs;

            % --- Zapisanie metryk wprost do SIGNALS.INFO ---
            if ~isempty(beat_str)
                beat_num = str2double(beat_str);
                if ~isnan(beat_num)
                    SIGNALS.INFO.(out_name).MORPHOLOGY_BEATS.(ext_str)(beat_num).SPECS = SPECS;
                end
            else
                SIGNALS.INFO.(out_name).MORPHOLOGY.(ext_str).SPECS = SPECS;
            end
        end
    end
end

% --- Zagnieżdżona funkcja pomocnicza ---
function val = parseTag(xmlBlock, tagName, defaultVal)
    pattern = ['<' tagName '>\s*([^<]+)\s*</' tagName '>'];
    tok = regexp(xmlBlock, pattern, 'tokens', 'once');
    if ~isempty(tok)
        % Bierzemy pierwszą liczbę w przypadku ciągów oddzielonych przecinkami
        parts = strsplit(tok{1}, ',');
        num = str2double(strtrim(parts{1}));
        if ~isnan(num)
            val = num;
            return;
        end
    end
    val = defaultVal;
end