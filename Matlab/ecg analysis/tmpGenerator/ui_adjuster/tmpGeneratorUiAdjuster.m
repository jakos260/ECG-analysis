classdef tmpGeneratorUiAdjuster < handle
    % service class of UI for adjusting parameters
    % for gets() -> the TMP generator

    properties
        fig
        tmpGeneratorHandle, qValuesHandle
        L12_matrix, L12_idx, L12_names;
        computeTimer
        colorMapGenerator
    end

    methods (Access = public)
        function obj = tmpGeneratorUiAdjuster(tmpGeneratorHandle, patient, offset, qValuesHandle)
            obj.colorMapGenerator = ColorMapGenerator();
            
            % close all
            close all force
            delete(findall(groot,'Type','uifigure'))
            drawnow
            
            obj.qValuesHandle = qValuesHandle;
            obj.fig = uifigure('Position',[100 100 600 1000]);

            obj.computeTimer = timer( ...
                'ExecutionMode','singleShot', ...
                'StartDelay',0.015, ...   % debounce delay (150 ms)
                'TimerFcn',@(~,~) obj.runComputation());
            obj.computeTimer.UserData = obj;


            obj.tmpGeneratorHandle = tmpGeneratorHandle;
            obj = obj.setL12Matrix();
            obj.loadEcgSimData(patient, offset);
            obj.createUI();
        end
    end

    methods (Access = private)

        % DATA LOADERS
        function loadEcgSimData(obj, patient, offset)
            obj.fig.UserData.dst        = loadmat(append(DATA_PATH, 'ECGsim_data/', patient, '/model/ventricle.dst3d'));
            obj.fig.UserData.V          = loadmat(append(DATA_PATH, 'ECGsim_data/', patient, '/model/ventricles2Thorax.mat'));
            obj.fig.UserData.A64        = loadmat(append(DATA_PATH, 'ECGsim_data/', patient, '/model/ventricles2BSM_(nijmegen_64).mat'));
            obj.fig.UserData.A12        = loadmat(append(DATA_PATH, 'ECGsim_data/', patient, '/model/ventricles2standard_12.mat'));
            BSM_ref                     = loadmat(append(DATA_PATH, 'ECGsim_data/', patient, '/ecgs/BSM_(nijmegen_64).refECG'));
            obj.fig.UserData.BSM_ref    = BSM_ref(:, offset:end);
            obj.fig.UserData.BSM_elc    = loadmat(append(DATA_PATH, 'ECGsim_data/', patient, '/ecgs/BSM_(nijmegen_64).elec'));
            STD_ref                     = loadmat(append(DATA_PATH, 'ECGsim_data/', patient, '/ecgs/standard_12.refECG'));
            obj.fig.UserData.STD_ref    = STD_ref(:, offset:end);
            obj.fig.UserData.STD_elc    = loadmat(append(DATA_PATH, 'ECGsim_data/', patient, '/ecgs/standard_12.elec'));
            obj.fig.UserData.dep        = loadmat(append(DATA_PATH, 'ECGsim_data/', patient, '/ventricular_beats/beat1/user.dep'));
            obj.fig.UserData.rep        = loadmat(append(DATA_PATH, 'ECGsim_data/', patient, '/ventricular_beats/beat1/user.rep'));
            obj.fig.UserData.L          = size(obj.fig.UserData.BSM_ref, 2);
            obj.fig.UserData.t_sim          = 0;
            obj.fig.UserData.t_ref          = 0;
            obj.fig.UserData.b              = 0.5;
            obj.fig.UserData.c              = 0.5;
            obj.fig.UserData.d              = 0.5;
            obj.fig.UserData.number_to_plot = 20;
            obj.fig.UserData.qtrtime        = 300;
            obj.fig.UserData.lead_number    = 8;
            obj.fig.UserData.showTOnly      = false;

        end

        % UI CREATOR
        function createUI(obj)
            % ADD PLOTS
            ax1 = uiaxes(obj.fig,'Position',[50 170 500 200]);
            ax1.NextPlot = 'add';
            ylim(ax1, [-0.2, 1.2]);
            
            ax2 = uiaxes(obj.fig,'Position',[50 420 500 250]);
            hold(ax2,'on')
            
            ax3 = uiaxes(obj.fig,'Position',[50 720 500 250]);
            hold(ax3,'on')

            obj.fig.UserData.ax2 = ax2;
            obj.fig.UserData.ax3 = ax3;

            obj.fig.UserData.ax4 = uiaxes( ...
                obj.fig, ...
                'Position', [310 720 240 250], ...
                'Visible', 'off', ...
                'Colormap', obj.getColorMap(), ...
                'XTick', [], 'YTick', [] ...
                );
            obj.fig.UserData.ax5 = uiaxes( ...
                obj.fig, ...
                'Position', [310 420 240 250], ...
                'Visible', 'off', ...
                'Colormap', obj.getColorMap(), ...
                'XTick', [], 'YTick', [] ...
                );

            obj.fig.UserData.tmp     = plot(ax1,nan,nan,'k');
            hold(ax2,'on')
            obj.fig.UserData.sig_ref = plot(ax2,nan,nan,'b','LineWidth',1);
            obj.fig.UserData.sig_sim = plot(ax2,nan,nan,'k','LineWidth',1.5);
            hold(ax2,'off')
            hold(ax2,'on')
            obj.fig.UserData.ecg_ref = plot(ax3,nan,nan,'r','LineWidth',1);
            obj.fig.UserData.ecg_sim = plot(ax3,nan,nan,'k','LineWidth',1.5);
            hold(ax2,'off')

            % ADD SLIDERS
            sld1 = uislider(obj.fig,...
                'Position',[100 700 400 3],...
                'Limits',[1 size(obj.fig.UserData.A64, 1)],...
                'Value',obj.fig.UserData.number_to_plot);
            
                % 'MajorTicks',0:20:round(mod(length(dst), 20) * 20),...
                % 'MinorTicks',[],...
            sld2 = uislider(obj.fig,...
                'Position',[100 150 400 3],...
                'Limits',[0 1],...
                'Value',obj.fig.UserData.b);
            sld3 = uislider(obj.fig,...
                'Position',[100 100 400 3],...
                'Limits',[0 1],...
                'Value',obj.fig.UserData.c);
            sld4 = uislider(obj.fig,...
                'Position',[100 50 400 3],...
                'Limits',[0 1],...
                'Value',obj.fig.UserData.d);
            
            sld5 = uislider(obj.fig,...
                'Position',[100 400 400 3],...
                'Limits',[1 obj.fig.UserData.L],...
                'Value',obj.fig.UserData.qtrtime);
            
            obj.addSlider(sld1,'number_to_plot','BSM(64) idx','%d',  true);
            obj.addSlider(sld2,'b','p0 (dep incl)','%.3f',          false);
            obj.addSlider(sld3,'c','p1 (st  decl)','%.3f',          false);
            obj.addSlider(sld4,'d','p2 (rep decl)','%.3f',          false);
            obj.addSlider(sld5,'qtrtime','time to model','%d',       true);

            % ADD L12_DROPDOWN
            dd = uidropdown(obj.fig, ...
                'Position',[70 970 80 25], ...
                'Items', obj.L12_names, ...
                'Value', 'V2', ...
                'ValueChangedFcn', @(btn,event) obj.dropdownCallback(btn));
            
            % ADD TOGGLE BUTTON (T-wave only)
            tgl = uibutton(obj.fig,'state', ...
                'Text','T-wave only (90:end)', ...
                'Position',[200 970 160 25], ...
                'Value',false, ...
                'ValueChangedFcn',@(btn,~) obj.toggleTCallback(btn));
            
            % ADD TOGGLE BUTTON (correlation map display)
            tglCorr = uibutton(obj.fig,'state', ...
                'Text','Extended correlation', ...
                'Position',[370 970 160 25], ... 
                'Value',false, ...
                'ValueChangedFcn',@(btn,~) obj.toggleCorrViewCallback(btn));

            % ADD ISOLINES
            obj.fig.UserData.iso2 = yline(ax2, 0, '--', ...
                'Color', [0.5 0.5 0.5], ...
                'LineWidth', 1);
            
            obj.fig.UserData.iso3 = yline(ax3, 0, '--', ...
                'Color', [0.5 0.5 0.5], ...
                'LineWidth', 1);

            % ADD TIME_MARKER
            obj.fig.UserData.time_marker = xline(ax1, 300, 'r');

            % ADD CORRELATION_DISPLAY
            obj.fig.UserData.corrText2 = text(ax2, ...
                0.98, 0.02, '', ...
                'Units','normalized', ...
                'FontWeight','bold', ...
                'VerticalAlignment','bottom', ...
                'HorizontalAlignment','right');
            
            obj.fig.UserData.corrText3 = text(ax3, ...
                0.98, 0.02, '', ...
                'Units','normalized', ...
                'FontWeight','bold', ...
                'VerticalAlignment','bottom', ...
                'HorizontalAlignment','right');

            % % ADD QTRIPLOT_BUTTON
            % btn = uibutton(fig, 'push', ...
            %     'Text', 'QTriplot', ...
            %     'Position', [450 900 90 25], ...
            %     'ButtonPushedFcn', @(src,~) btnCallback(fig, qWakeHandle));
            
            % INITIAL DRAW
            obj.sliderCallback( ...
                'number_to_plot', ...
                obj.fig.UserData.number_to_plot, ...
                [], ...          % label not needed on init
                '%d', ...
                true);
        end

        function addSlider(obj, sld, param, name, format, isInteger)

            % Name label
            uilabel(obj.fig,...
                'Text',name,...
                'Position',[20 sld.Position(2)-8 80 22]);
        
            % Value label
            lblVal = uilabel(obj.fig,...
                'Text',sprintf(format,sld.Value),...
                'Position',[470 sld.Position(2)-8 70 22],...
                'HorizontalAlignment','right');
        
            % Single callback
            sld.ValueChangingFcn = @(s,e) obj.sliderCallback(param,e.Value,lblVal,format,isInteger);
        end

        % COMPUTATION WRAPPER
        function runComputation(obj)        
            ud = obj.fig.UserData;

            % --- RE-ADDED: Get target indexes ---
            idx = ud.number_to_plot;
            lead = ud.lead_number;
        
            [tmp_all, sig_sim_all] = obj.tmpGeneratorHandle( ...
                [ud.b, ud.c, ud.d], ...
                ud.dep, ud.rep, ud.L, ud.A64);
        
            [~, std_sim_all] = obj.tmpGeneratorHandle( ...
                [ud.b, ud.c, ud.d], ...
                ud.dep, ud.rep, ud.L, ud.A12);

            % --- RE-ADDED: Slice the matrices for plotting ---
            if ud.showTOnly
                range = 90:ud.L;
            else
                range = 1:ud.L;
            end
            
            % BSM64
            tmp_sim = tmp_all(idx, range);
            sig_sim = sig_sim_all(idx, range);
            sig_ref = ud.BSM_ref(idx, range);
            
            % EKG12
            ecg12_sim = std_sim_all(:, range);
            ecg12_ref = ud.STD_ref(:, range);
            
            sel_ecg12_sim = ecg12_sim(lead, :);
            sel_ecg12_ref = ecg12_ref(lead, :);
        
            % --- plotting & correlation update ---
            set(ud.tmp, 'XData', 1:length(tmp_sim), 'YData', tmp_sim);
            
            set(ud.sig_sim, 'XData', range, 'YData', sig_sim);
            set(ud.sig_ref, 'XData', range, 'YData', sig_ref);

            set(ud.ecg_ref, 'XData', range, 'YData', sel_ecg12_ref);
            set(ud.ecg_sim, 'XData', range, 'YData', sel_ecg12_sim);
        
            ud.time_marker.Value = ud.qtrtime;
        
            % --- Update correlation coefficients ---
            corr_BSM        = obj.getSignalsCorr(ud.BSM_ref, sig_sim_all, false);
            corr2_global    = mean(corr_BSM);
            corr2_local     = obj.getSignalsCorr(sig_ref, sig_sim);
            corr3_global    = obj.getSignalsCorr(ecg12_ref, ecg12_sim);
            corr3_local     = obj.getSignalsCorr(sel_ecg12_ref, sel_ecg12_sim);

            set(ud.corrText2, 'String', sprintf('corr(full) = %.3f\ncorr(local) = %.3f', corr2_global, corr2_local));
            set(ud.corrText3, 'String', sprintf('corr(full) = %.3f\ncorr(local) = %.3f', corr3_global, corr3_local));        
            
            if strcmp(ud.ax4.Visible, 'on')
                imagesc(ud.ax4, obj.getECG12Map(ecg12_ref, ecg12_sim), [-1, 1]);
                axis(ud.ax4, 'tight');
                set(ud.ax4, 'XTick', [], 'YTick', []);

                colorMap = obj.colorMapGenerator.BSM64Corr2Img(corr_BSM);
                imagesc(ud.ax5, colorMap); 
                axis(ud.ax5, 'tight');
                axis(ud.ax5, 'off');
            end

            % qtriplot update
            [sim_dep, sim_rep] = obj.getTimesFromS(tmp_all);
            ud.t_sim = sim_rep; % sim_rep, tmp_all(:,ud.qtrtime)
            if(0 ~= size(ud.t_sim))
                obj.qValuesHandle(ud.t_sim);
            end
        
            % assign modified UserData to figure and draw
            obj.fig.UserData = ud;
            drawnow limitrate nocallbacks
        end

        % CALLBACKS
        function sliderCallback(obj, param,val,lbl,format,isInteger)

            if isInteger
                val = round(val);
            end
            if ~isempty(lbl)
                lbl.Text = sprintf(format,val);
            end
        
            obj.fig.UserData.(param) = val;

            % debouncing
            if strcmp(obj.computeTimer.Running,'on')
                stop(obj.computeTimer);
            end
            start(obj.computeTimer);
        end
        
        function dropdownCallback(obj, src)
        
            % Get selected lead index
            lead_idx = src.ValueIndex;
        
            % Store it
            obj.fig.UserData.lead_number = lead_idx;
        
            % Trigger redraw using existing logic
            obj.sliderCallback( ...
                'number_to_plot', ...
                obj.fig.UserData.number_to_plot, ...
                [], ...
                '%d', ...
                true);
        end
        
        function btnCallback(obj)
            obj.qValuesHandle(obj.fig.UserData.t_sim)
        end

        function toggleTCallback(obj, btn)

            obj.fig.UserData.showTOnly = btn.Value;
            obj.sliderCallback( ...
                'number_to_plot', ...
                obj.fig.UserData.number_to_plot, ...
                [], ...
                '%d', ...
                true);
        
        end

        function toggleCorrViewCallback(obj, btn)
            ud = obj.fig.UserData;
            isAdvView = btn.Value;
            
            if isAdvView
                % Zwężenie lewych wykresów o połowę
                ud.ax3.Position = [50 720 240 250];
                ud.ax2.Position = [50 420 240 250];
                
                % Pokazanie prawych wykresów
                ud.ax4.Visible = 'on';
                ud.ax5.Visible = 'on';
                
                % Pokazanie zawartości (jeśli coś już jest tam narysowane)
                if ~isempty(ud.ax4.Children)
                    set(ud.ax4.Children, 'Visible', 'on');
                end
                if ~isempty(ud.ax5.Children)
                    set(ud.ax5.Children, 'Visible', 'on');
                end
            else
                % Przywrócenie pierwotnych rozmiarów lewych wykresów
                ud.ax3.Position = [50 720 500 250];
                ud.ax2.Position = [50 420 500 250];
                
                % Ukrycie prawych wykresów
                ud.ax4.Visible = 'off';
                ud.ax5.Visible = 'off';
                
                % Ukrycie zawartości (aby obrazy/wykresy nie nachodziły)
                if ~isempty(ud.ax4.Children)
                    set(ud.ax4.Children, 'Visible', 'off');
                end
                if ~isempty(ud.ax5.Children)
                    set(ud.ax5.Children, 'Visible', 'off');
                end
            end
            
            % Zapisanie stanu i wymuszenie odświeżenia
            obj.fig.UserData = ud;
            
            % Opcjonalnie: wywołanie przeliczenia, aby narysować mapy, 
            % jeśli są generowane w czasie rzeczywistym
            obj.sliderCallback( ...
                'number_to_plot', ...
                ud.number_to_plot, ...
                [], '%d', true);
        end

        % HELPERS
        function obj = setL12Matrix(obj)

            obj.L12_idx.V1 = 1;
            obj.L12_idx.V2 = 2;
            obj.L12_idx.V3 = 3;
            obj.L12_idx.V4 = 4;
            obj.L12_idx.V5 = 5;
            obj.L12_idx.V6 = 6;
            obj.L12_idx.RA = 7;
            obj.L12_idx.LA = 8;
            obj.L12_idx.LL = 9;
        
            obj.L12_matrix = zeros(12,64);
        
            % Lead I
            obj.L12_matrix(1,obj.L12_idx.LA) =  1;
            obj.L12_matrix(1,obj.L12_idx.RA) = -1;
        
            % Lead II
            obj.L12_matrix(2,obj.L12_idx.LL) =  1;
            obj.L12_matrix(2,obj.L12_idx.RA) = -1;
        
            % Lead III
            obj.L12_matrix(3,obj.L12_idx.LL) =  1;
            obj.L12_matrix(3,obj.L12_idx.LA) = -1;
        
            % aVR
            obj.L12_matrix(4,obj.L12_idx.RA) =  1;
            obj.L12_matrix(4,obj.L12_idx.LA) = -0.5;
            obj.L12_matrix(4,obj.L12_idx.LL) = -0.5;
        
            % aVL
            obj.L12_matrix(5,obj.L12_idx.LA) =  1;
            obj.L12_matrix(5,obj.L12_idx.RA) = -0.5;
            obj.L12_matrix(5,obj.L12_idx.LL) = -0.5;
        
            % aVF
            obj.L12_matrix(6,obj.L12_idx.LL) =  1;
            obj.L12_matrix(6,obj.L12_idx.RA) = -0.5;
            obj.L12_matrix(6,obj.L12_idx.LA) = -0.5;
        
            % V1–V6
            for k = 1:6
                obj.L12_matrix(6+k, obj.L12_idx.(['V' num2str(k)])) = 1;
                obj.L12_matrix(6+k, obj.L12_idx.RA) = -1/3;
                obj.L12_matrix(6+k, obj.L12_idx.LA) = -1/3;
                obj.L12_matrix(6+k, obj.L12_idx.LL) = -1/3;
            end
        
            % names
            obj.L12_names = ["I", "II", "III", "aVR", "aVL", "aVF", "V1", "V2", "V3", "V4", "V5", "V6"];
        
        end
        
        function [dep, rep] = getTimesFromS(obj, S)
            dep = first_crossing_thr(S, 0.5, true);
            rep = first_crossing_thr(S, 0.5, false);
        end

        function corr = getSignalsCorr(obj, sig1, sig2, common)
            arguments
                obj
                sig1
                sig2
                common logical = true  % Domyślnie ustawione na true
            end
            len = min(size(sig1,2), size(sig2,2));
            if common
                tmp = corrcoef(sig1(:,1:len), sig2(:,1:len));
                corr = tmp(1,2);
            else
                numSignals = size(sig1, 1);
                corr = zeros(1, numSignals);
                for i = 1:numSignals
                    R = corrcoef(sig1(i,1:len), sig2(i,1:len));
                    corr(i) = R(1,2);
                end
            end
        end

        function map = getECG12Map(obj, sig1, sig2)
            map_cols = size(sig1,1)/4; % divisor is number of rows
            map_rows = size(sig1,1)/map_cols;
            map = zeros(map_rows, map_cols);
            for r = 1:map_rows
                for c = 1:map_cols
                    signal_num = (r-1) * map_cols + c;
                    map(r,c) = obj.getSignalsCorr(sig1(signal_num,:), sig2(signal_num,:));
                end
            end
        end

        function rgb = getColorMap(obj)
            n = 64; % rozdzielczość mapy kolorów (musi być parzysta)
            r = [linspace(1, 1, n/2), linspace(1, 0, n/2)]';
            g = [linspace(0, 1, n/2), linspace(1, 1, n/2)]';
            b = [linspace(0, 1, n/2), linspace(1, 0, n/2)]';
            rgb = [r, g, b];
        end
    end
end