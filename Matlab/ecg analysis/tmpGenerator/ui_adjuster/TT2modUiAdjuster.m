classdef TT2modUiAdjuster < handle
    % Zaktualizowana wersja UI z czystymi suwakami (tylko min, 1, max)
    
    properties
        fig
        tmpGeneratorHandle, qValuesHandle
        L12_names = ["I", "II", "III", "aVR", "aVL", "aVF", "V1", "V2", "V3", "V4", "V5", "V6"];
        computeTimer
        sliderConfigs
        sliders = {}
        
        % Uchwyty UI
        ax_epi, ax_endo, ax_ecg, ax_corr
        cb_epi, cb_endo, dd_lead, tgl_twave
        
        % Obiekty graficzne
        line_epi_ref, line_epi_mod
        line_endo_ref, line_endo_mod
        line_ecg_ref, line_ecg_sim
        shade_ecg
        im_corr
        corr_texts
    end

    methods
        function obj = TT2modUiAdjuster(tmpGeneratorHandle, patient, offset, qValuesHandle, sliderConfigs)
            obj.tmpGeneratorHandle = tmpGeneratorHandle;
            obj.qValuesHandle = qValuesHandle;
            obj.sliderConfigs = sliderConfigs;
            
            delete(findall(groot,'Type','uifigure', 'Name', 'TT2mod Tuner V2'))
            
            obj.fig = uifigure('Name', 'TT2mod Tuner V2', 'Position', [100 100 700 850]);
            
            obj.computeTimer = timer('ExecutionMode', 'singleShot', 'StartDelay', 0.05, ...
                'TimerFcn', @(~,~) obj.runComputation());

            obj.loadEcgSimData(patient, offset);
            obj.createUI();
            obj.runComputation(); 
        end
    end

    methods (Access = private)
        function loadEcgSimData(obj, patient, offset)
            obj.fig.UserData.A12 = loadmat(append(DATA_PATH, 'ECGsim_data/', patient, '/model/ventricles2standard_12.mat'));
            STD_ref = loadmat(append(DATA_PATH, 'ECGsim_data/', patient, '/ecgs/standard_12.refECG'));
            obj.fig.UserData.STD_ref = STD_ref(:, offset:end);
            obj.fig.UserData.dep = loadmat(append(DATA_PATH, 'ECGsim_data/', patient, '/ventricular_beats/beat1/user.dep'));
            obj.fig.UserData.rep = loadmat(append(DATA_PATH, 'ECGsim_data/', patient, '/ventricular_beats/beat1/user.rep'));
            obj.fig.UserData.L = size(obj.fig.UserData.STD_ref, 2);
        end

        function createUI(obj)
            % 1. Panel kontrolny (Góra)
            uilabel(obj.fig, 'Text', 'Lead:', 'Position', [20 810 40 22]);
            obj.dd_lead = uidropdown(obj.fig, 'Items', obj.L12_names, 'Value', 'V2', ...
                'Position', [60 810 70 22], 'ValueChangedFcn', @(src,~) obj.triggerUpdate());
            
            % Przycisk fali T
            obj.tgl_twave = uibutton(obj.fig, 'state', 'Text', 'T-wave Corr', ...
                'Position', [150 810 100 22], 'Value', false, ...
                'ValueChangedFcn', @(~,~) obj.triggerUpdate());
            
            % 2. Wykres EKG
            obj.ax_ecg = uiaxes(obj.fig, 'Position', [50 600 600 200]);
            title(obj.ax_ecg, 'Selected ECG Lead'); 
            hold(obj.ax_ecg, 'on');
            
            % Cieniowanie dla obszaru 0-100 (początkowo ukryte, w tle)
            obj.shade_ecg = patch(obj.ax_ecg, [0 100 100 0], [-5 -5 5 5], 'k', ...
                'FaceAlpha', 0.15, 'EdgeColor', 'none', 'Visible', 'off');
            
            obj.line_ecg_ref = plot(obj.ax_ecg, 0, 0, 'k--', 'LineWidth', 1);
            obj.line_ecg_sim = plot(obj.ax_ecg, 0, 0, 'r', 'LineWidth', 1.5);

            % 3. Mapa korelacji (Czysta)
            obj.ax_corr = uiaxes(obj.fig, 'Position', [50 440 600 140]);
            obj.im_corr = imagesc(obj.ax_corr, zeros(2,6));
            cmap = [linspace(1,0,256)', linspace(0,1,256)', zeros(256,1)];
            colormap(obj.ax_corr, cmap);
            clim(obj.ax_corr, [0 1]);
            obj.ax_corr.XTick = []; obj.ax_corr.YTick = [];
            obj.ax_corr.XColor = 'none'; obj.ax_corr.YColor = 'none';
            
            obj.corr_texts = gobjects(2,6);
            for r = 1:2
                for c = 1:6
                    obj.corr_texts(r,c) = text(obj.ax_corr, c, r, '', ...
                        'HorizontalAlignment', 'center', 'FontWeight', 'bold');
                end
            end

            % 4. Checkboxy AP
            obj.cb_epi = uicheckbox(obj.fig, 'Text', 'Mod.', 'Value', 1, ...
                'Position', [160 400 60 22], 'ValueChangedFcn', @(~,~) obj.triggerUpdate());
            obj.cb_endo = uicheckbox(obj.fig, 'Text', 'Mod.', 'Value', 1, ...
                'Position', [480 400 60 22], 'ValueChangedFcn', @(~,~) obj.triggerUpdate());

            % 5. Subploty AP
            obj.ax_epi = uiaxes(obj.fig, 'Position', [50 200 280 180]);
            title(obj.ax_epi, 'AP Epicardium'); 
            obj.ax_endo = uiaxes(obj.fig, 'Position', [370 200 280 180]);
            title(obj.ax_endo, 'AP Endocardium'); 
            hold([obj.ax_epi, obj.ax_endo], 'on');
            
            obj.line_epi_ref = plot(obj.ax_epi, 0, 0, 'Color', [0.6 0.6 0.6], 'LineWidth', 2);
            obj.line_epi_mod = plot(obj.ax_epi, 0, 0, 'r', 'LineWidth', 1.5);
            obj.line_endo_ref = plot(obj.ax_endo, 0, 0, 'Color', [0.6 0.6 0.6], 'LineWidth', 2);
            obj.line_endo_mod = plot(obj.ax_endo, 0, 0, 'b', 'LineWidth', 1.5);

            % 6. Slidery
            nS = length(obj.sliderConfigs);
            for i = 1:nS
                cfg = obj.sliderConfigs(i);
                yPos = 140 - (i-1)*45;
                uilabel(obj.fig, 'Text', cfg.Name, 'Position', [20 yPos 100 22]);
                
                % Wyznaczenie skrajnych ticków oraz wartości 1 (unikalne dla bezpieczeństwa)
                m_ticks = unique([cfg.Limits(1), 1, cfg.Limits(2)]);
                
                sld = uislider(obj.fig, 'Limits', cfg.Limits, 'Value', cfg.Default, ...
                    'MajorTicks', m_ticks, 'MinorTicks', [], ...  % Ograniczenie widocznych podziałek
                    'Position', [120 yPos+10 450 3], 'ValueChangedFcn', @(~,~) obj.triggerUpdate());
                
                lblVal = uilabel(obj.fig, 'Text', num2str(cfg.Default, '%.2f'), 'Position', [580 yPos 50 22]);
                obj.sliders{i} = struct('handle', sld, 'label', lblVal);
            end
        end

        function triggerUpdate(obj)
            if strcmp(obj.computeTimer.Running, 'on'), stop(obj.computeTimer); end
            start(obj.computeTimer);
        end

        function runComputation(obj)
            ud = obj.fig.UserData;
            p = cellfun(@(s) s.handle.Value, obj.sliders);
            for i=1:length(obj.sliders), obj.sliders{i}.label.Text = num2str(p(i), '%.2f'); end

            norm_fn = @(v) (v - min(v)) / (max(v) - min(v));
            
            % AP Templates
            [t_ref_epi, v_ref_epi] = wrapper_TenTusscher2mod(0.1, 500, 1, [1,1,1]);
            [t_ref_endo, v_ref_endo] = wrapper_TenTusscher2mod(0.1, 500, 3, [1,1,1]);
            set(obj.line_epi_ref, 'XData', t_ref_epi, 'YData', norm_fn(v_ref_epi));
            set(obj.line_endo_ref, 'XData', t_ref_endo, 'YData', norm_fn(v_ref_endo));

            p_epi = [1,1,1]; if obj.cb_epi.Value, p_epi = p; end
            p_endo = [1,1,1]; if obj.cb_endo.Value, p_endo = p; end

            if obj.cb_epi.Value
                [t_m_epi, v_m_epi] = wrapper_TenTusscher2mod(0.1, 500, 1, p_epi);
                set(obj.line_epi_mod, 'Visible', 'on', 'XData', t_m_epi, 'YData', norm_fn(v_m_epi));
            else, set(obj.line_epi_mod, 'Visible', 'off'); end

            if obj.cb_endo.Value
                [t_m_endo, v_m_endo] = wrapper_TenTusscher2mod(0.1, 500, 3, p_endo);
                set(obj.line_endo_mod, 'Visible', 'on', 'XData', t_m_endo, 'YData', norm_fn(v_m_endo));
            else, set(obj.line_endo_mod, 'Visible', 'off'); end

            % Generator EKG
            [tmp_all, ecg_sim_all] = obj.tmpGeneratorHandle(p_epi, p_endo, ud.dep, ud.rep, ud.L, ud.A12);

            % Update wybranego odprowadzenia
            idx = obj.dd_lead.ValueIndex;
            set(obj.line_ecg_ref, 'XData', 1:ud.L, 'YData', ud.STD_ref(idx,:));
            set(obj.line_ecg_sim, 'XData', 1:ud.L, 'YData', ecg_sim_all(idx,:));

            % Dopasowanie limitów cieniowania
            yl = ylim(obj.ax_ecg);
            set(obj.shade_ecg, 'YData', [yl(1) yl(1) yl(2) yl(2)]);

            % LOGIKA KORELACJI I WIDOCZNOŚĆ CIENIOWANIA
            if obj.tgl_twave.Value
                range = 100:ud.L;
                title(obj.ax_corr, '12-Lead Correlation (T-wave)');
                set(obj.shade_ecg, 'Visible', 'on');
            else
                range = 1:ud.L;
                title(obj.ax_corr, '12-Lead Correlation');
                set(obj.shade_ecg, 'Visible', 'off');
            end

            corrs = zeros(1, 12);
            for j = 1:12
                r = corrcoef(ud.STD_ref(j, range), ecg_sim_all(j, range));
                corrs(j) = r(1,2);
            end
            
            corrs_mat = [corrs(1:6); corrs(7:12)];
            obj.im_corr.CData = corrs_mat;
            lead_names = ["I", "II", "III", "aVR", "aVL", "aVF"; "V1", "V2", "V3", "V4", "V5", "V6"];
            for r = 1:2
                for c = 1:6
                    val = corrs_mat(r,c);
                    tColor = 'w'; if val > 0.5, tColor = 'k'; end
                    obj.corr_texts(r,c).String = sprintf('%s\n%.2f', lead_names(r,c), val);
                    obj.corr_texts(r,c).Color = tColor;
                end
            end
            
            if ~isempty(obj.qValuesHandle), obj.qValuesHandle(obj.getRepolarizationTimes(tmp_all)); end
            drawnow limitrate
        end

        function rep = getRepolarizationTimes(obj, S)
            rep = zeros(size(S,1), 1);
            for i = 1:size(S,1)
                idx = find(S(i,:) < 0.5, 1, 'first');
                if isempty(idx), rep(i) = 0; else, rep(i) = idx; end
            end
        end
    end
end