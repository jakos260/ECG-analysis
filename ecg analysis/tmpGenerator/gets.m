function S=gets(T,dep,rep,p,mode)

    [nn, nt] = size(T);
    
    if size(dep,1)~=nn, dep=dep'; end
    if size(rep,1)~=nn, rep=rep'; end
    
    TDEP=T-dep*ones(1,nt);

    if mode==getsMode.Exponent
        S=1./(1+exp(-p(1)*TDEP));
    end
    
    if mode==getsMode.LineraST
        TDEP=T-dep*ones(1,nt);
        TREP=T-rep*ones(1,nt)-p(4);
        
        % --- upslope ---
        upslope =(1-p(2))./(1+exp(-p(1)*TDEP)); % apply the up-slope: logistic shape
        upslope = bsxfun(@rdivide,upslope,max(upslope,[],2)); % normalize

        % --- downslope ---
        downslope = 1 - (1./(1+exp(p(3)*TREP)));

        % --- st_section ---
        [n_rows, n_cols] = size(TREP);
        
        % st_section config
        if (length(p) == 6)
            val_start   = 1;
            val_end     = p(5);
            st_length   = p(6);
        elseif (length(p) == 5)
            val_start   = 1;
            [val_end, st_length] = get_dependent_xy_values(p(5));
        else
            error(sprintf('gets_v2() in mode = %d require 5 or 6 values in p, %d passed', mode, length(p)));
        end
        
        % generate linear st_section
        st_section = linspace(val_start, val_end, st_length)';
        S = upslope;

        % calculate index of first signal sample above treshold
        idx_up = first_crossing_thr(S, 0.999, true);
        
        for i = 1:n_rows
            % here ends upslope and st_section begins
            end_up_idx = idx_up(i);
            
            % here we start downslope, downslope samples difference on this
            % idx match st_section difference so transition will be smooth
            start_down_idx = first_samples_diff(downslope(i,:) * val_end, st_section(end-1)-st_section(end));

            % add st_section
            S(i,end_up_idx:end_up_idx+st_length-1) = S(i,end_up_idx:end_up_idx+st_length-1) .* st_section';
        
            % add downslope
            paste_length = min(n_cols - (end_up_idx+st_length), n_cols - start_down_idx);
            S(i,end_up_idx+st_length-1:end_up_idx+st_length-1 + paste_length) = (val_end / downslope(i,start_down_idx)) * downslope(i,start_down_idx:start_down_idx + paste_length);
            S(i,end_up_idx+st_length-1 + paste_length : end) = 0;
        end
    end

    if mode==getsMode.CosST
        % warning('experimental mode');

        TDEP=T-dep*ones(1,nt);
        TREP=T-rep*ones(1,nt)-p(4);
        
        % --- upslope ---
        upslope =(1-p(2))./(1+exp(-p(1)*TDEP)); % apply the up-slope: logistic shape
        upslope = bsxfun(@rdivide,upslope,max(upslope,[],2)); % normalize

        % --- downslope ---
        downslope = 1 - (1./(1+exp(p(3)*TREP)));

        % --- st_section ---
        
        % st_section config
        if (length(p) == 6)
            val_start   = 1;
            val_end     = p(5);
            st_length   = p(6);
        elseif (length(p) == 5)
            val_start   = 1;
            [val_end, st_length] = get_dependent_xy_values(p(5));
        else
            error(sprintf('gets_v2() in mode = %d require 5 or 6 values in p, %d passed', mode, length(p)));
        end
        
        % generate linear st_section
        st_section_x = linspace(0, pi/8, st_length)';
        st_corection = linspace(0, cos(pi/8), st_length)';
        st_section = (cos(st_section_x)-st_corection)*(val_start - val_end) + val_end;
        S = upslope;

        % calculate index of first signal sample above treshold
        [n_rows, n_cols] = size(TREP);
        idx_up = first_crossings_thr(upslope, 0.9999);
        
        for i = 1:n_rows
            % here ends upslope and st_section begins
            end_up_idx = idx_up(i);
            
            % here we start downslope, downslope samples difference on this
            % idx match st_section difference so transition will be smooth
            start_down_idx = first_samples_diff(downslope(i,:) * val_end, st_section(end-1)-st_section(end));

            % add st_section
            S(i,end_up_idx:end_up_idx+st_length-1) = S(i,end_up_idx:end_up_idx+st_length-1) .* st_section';
        
            % add downslope
            paste_length = min(n_cols - (end_up_idx+st_length), n_cols - start_down_idx);
            S(i,end_up_idx+st_length-1:end_up_idx+st_length-1 + paste_length) = (val_end / downslope(i,start_down_idx)) * downslope(i,start_down_idx:start_down_idx + paste_length);
            S(i,end_up_idx+st_length-1 + paste_length : end) = 0;
        end
    end

    if mode==getsMode.LinST_rept_dept_ang
        warning("uncontinous function");
        TDEP=T-dep*ones(1,nt);
        TREP=T-rep*ones(1,nt)-p(4);
        
        % --------------------------- upslope -----------------------------
        % apply the up-slope: logistic shape
        upslope =(1-p(2))./(1+exp(-p(1)*TDEP)); 
        upslope = bsxfun(@rdivide,upslope,max(upslope,[],2)); % normalize

        % --------------------------- downslope ---------------------------
        downslope = 1 - (1./(1+exp(p(3)*TREP)));

        % --------------------------- st_section --------------------------
        % get "first visible" highest point as ST begin
        st_begin_idx = first_crossing_up_thr(upslope, 0.999);

        % create ST section template
        st_section = ones(size(downslope));
        
        % calculate ST lenght based on provided angle
        st_x = ceil(cot(p(5)*pi/4) * 500); % 1y equals to 500x units

        % if ST section meets 0 in finite distance go there
        if(isreal(st_x) && isfinite(st_x))
            % pad pattern (stałe dla wszystkich wierszy)
            pad = linspace(1, 0, st_x);
            
            % pad length per row
            pad_len = min(st_x, nt - st_begin_idx);
            
            % compute absolute linear indices into st_section
            rows = (1:nn)';         % [nn x 1]
            maxPad = max(pad_len);  % maximum pad length across rows

            if maxPad > 0
                % generate column indices for each row
                cols = st_begin_idx + (0:maxPad-1);
            
                % mask for valid pad positions (within pad_len for each row)
                valid_mask = bsxfun(@le, (0:maxPad-1), (pad_len-1));
            
                % clip columns to avoid out-of-bounds indices
                cols_clipped = cols;
                cols_clipped(~valid_mask) = 1;
            
                % convert row/column indices to linear indices
                lin_idx = sub2ind([nn nt], repelem(rows, 1, maxPad), cols_clipped);
            
                % replicate pad pattern for each row
                pad_mat = repmat(pad(1:maxPad), nn, 1);   % nn x maxPad
            
                % select only valid pad elements
                pad_vals = pad_mat(valid_mask);
            
                % assign pad values to st_section
                st_section(lin_idx(valid_mask)) = pad_vals;
            else
                warning("padding failed, using flat ST section");
                st_section = ones(size(downslope));
            end


            % --------------- downslope position correction ---------------
            % original minimal diff index
            [~, stepest_downslope_idx] = min(diff(downslope, 1, 2), [], 2);

            % modified minimal diff index
            [~, stepest_mod_downslope_idx] = min(diff(downslope .* st_section, 1, 2), [], 2);

            % replace NaN indices (no step found) with safe default 1
            stepest_downslope_idx(~isfinite(stepest_downslope_idx)) = 1;
            stepest_mod_downslope_idx(~isfinite(stepest_mod_downslope_idx)) = 1;

            % compute shift
            correction = stepest_downslope_idx - stepest_mod_downslope_idx;

            % adjust column indices -> want matrix [nn x nt]
            % start with a vector of base columns (1:nt)
            base_cols = 1:nt;                    % 1 x nt

            % produce col_idx so that for row i, col j -> col = base_cols(j) - correction(i)
            % create matrix by subtracting correction per row
            col_idx = bsxfun(@minus, base_cols, correction(:));   % nn x nt

            % clip out-of-bounds
            col_idx(col_idx < 1) = 1;
            col_idx(col_idx > nt) = nt;

            % ensure integer indices
            col_idx = round(col_idx);

            % build row indices matrix [nn x nt]
            row_idx = repmat((1:nn).', 1, nt);

            % final safety: sizes must match
            if ~isequal(size(row_idx), size(col_idx))
                error('Internal indexing size mismatch: row_idx vs col_idx');
            end

            % safer sub2ind indexing
            lin_idx = sub2ind([nn nt], row_idx, col_idx);
            corrected_downslope = downslope(lin_idx);

            % combine all segments
            S = upslope .* corrected_downslope .* st_section;
        else
            % if ST section doesn't meet 0 there is no need to calculate it
            S = upslope .* downslope;
        end
        % normalize result S
        S = bsxfun(@rdivide,S,max(S,[],2));
        % S = (S-min(S))./(max(S)-min(S));
    end

    if mode==getsMode.LinST_rept_ang_dept
        TDEP=T-dep*ones(1,nt);
        TREP=T-rep*ones(1,nt)-p(4);
        
        % --------------------------- upslope -----------------------------
        % apply the up-slope: logistic shape
        upslope =(1-p(2))./(1+exp(-p(1)*TDEP)); 
        upslope = bsxfun(@rdivide,upslope,max(upslope,[],2)); % normalize

        % --------------------------- downslope ---------------------------
        downslope = 1 - (1./(1+exp(p(5)*TREP)));

        % --------------------------- st_section --------------------------
        % get "first visible" highest point as ST begin
        st_begin_idx = first_crossing_up_thr(upslope, 0.99);
        st_end_idx = first_crossing_down_thr(downslope, 0.99);

        % create ST section template
        st_section = ones(size(downslope));
        
        % calculate ST lenght
        for i = 1:nn
            idx_0 = min(st_begin_idx(i), nt);
            idx_1 = min(max(st_end_idx(i), idx_0+1), nt);
            if (p(3) ~= 0)
                st_pad_length = ceil((idx_1-idx_0)/p(3));
                st_pad_end = min(idx_0 + st_pad_length, size(st_section, 2));
                pad = linspace(1, 0, st_pad_length);
                st_section(i, idx_0:st_pad_end-1) = pad(1:st_pad_end-idx_0);
                try 
                    st_section(i, st_pad_end:end) = 0;
                catch e
                    sprintf e;
                end
            end
        end

        % --------------- downslope position correction -------------------

        % original minimal diff index
        [stepest_downslope, stepest_downslope_idx] = min(diff(downslope, 1, 2), [], 2);

        % modified minimal diff index
        [stepest_mod_downslope, stepest_mod_downslope_idx] = min(diff(downslope .* st_section, 1, 2), [], 2);

        % replace NaN indices (no step found) with safe default 1
        stepest_downslope_idx(~isfinite(stepest_downslope_idx)) = 1;
        stepest_mod_downslope_idx(~isfinite(stepest_mod_downslope_idx)) = 1;

        % compute shift
        correction = stepest_downslope_idx - stepest_mod_downslope_idx;

        % adjust column indices -> want matrix [nn x nt]
        % start with a vector of base columns (1:nt)
        base_cols = 1:nt;                    % 1 x nt

        % produce col_idx so that for row i, col j -> col = base_cols(j) - correction(i)
        % create matrix by subtracting correction per row
        col_idx = bsxfun(@minus, base_cols, correction(:));   % nn x nt

        % clip out-of-bounds
        col_idx(col_idx < 1) = 1;
        col_idx(col_idx > nt) = nt;

        % ensure integer indices
        col_idx = round(col_idx);

        % build row indices matrix [nn x nt]
        row_idx = repmat((1:nn).', 1, nt);

        % final safety: sizes must match
        if ~isequal(size(row_idx), size(col_idx))
            error(sprintf('Internal indexing size mismatch:\nrow_idx vs col_idx -> %s vs %s\n', mat2str(size(row_idx)), mat2str(size(col_idx))));
        end

        % safer sub2ind indexing
        lin_idx = sub2ind([nn nt], row_idx, col_idx);
        
        % get differences ratio modified_downslope/downslope
        amp_corr = stepest_mod_downslope./stepest_downslope;

        % apply correction for amplitude and time
        tmp_downslope = 1 - (1./(1+exp((p(5)*TREP)./amp_corr)));
        corrected_downslope = tmp_downslope(lin_idx);

        % combine all segments
        S = upslope .* corrected_downslope .* st_section;

        % normalize result S
        S = bsxfun(@rdivide,S,max(S,[],2));
        % S = (S-min(S))./(max(S)-min(S));
    end

    if mode==getsMode.NoReptCorrection
        p(3) = p(3)/10;
        p(5) = p(5) + 0.002;
        TDEP=T-dep*ones(1,nt);
        TREP=T-rep*ones(1,nt)-p(4);
        
        % --------------------------- upslope -----------------------------
        % apply the up-slope: logistic shape
        upslope =(1-p(2))./(1+exp(-p(1)*TDEP)); 
        upslope = bsxfun(@rdivide,upslope,max(upslope,[],2)); % normalize

        % --------------------------- downslope ---------------------------
        downslope = 1 - (1./(1+exp(p(5)*TREP)));

        % --------------------------- st_section --------------------------
        % get "first visible" highest point as ST begin
        st_begin_idx = first_crossing_up_thr(upslope, 0.999);
        st_end_idx = first_crossing_down_thr(downslope, 0.999);

        % create ST section template
        st_section = ones(size(downslope));

        % calculate ST lenght
        for i = 1:nn
            idx_0 = min(st_begin_idx(i), nt);
            idx_1 = min(max(st_end_idx(i), idx_0+1), nt);
            if (p(3) ~= 0)
                st_pad_length = ceil((idx_1-idx_0)/p(3));
                st_pad_end = min(idx_0 + st_pad_length, size(st_section, 2));
                pad = linspace(1, 0, st_pad_length);
                st_section(i, idx_0:st_pad_end-1) = pad(1:st_pad_end-idx_0);
                try 
                    st_section(i, st_pad_end:end) = 0;
                catch e
                    sprintf e;
                end
            end
        end

        % % vectorized calculation of st section
        % idx_0 = min(st_begin_idx(:), nt);
        % idx_1 = min(max(st_end_idx(:), idx_0 + 1), nt);
        % st_pad_length = ceil((idx_1 - idx_0) / p(3));
        % st_pad_end = min(idx_0 + st_pad_length, nt);
        % col_indices = 1:nt;
        % 
        % mask_fade = (col_indices >= idx_0) & (col_indices < st_pad_end);
        % mask_zero = (col_indices >= st_pad_end);
        % 
        % fade_values = 1 - (col_indices - idx_0) ./ st_pad_length;
        % st_section = ones(nn, nt);
        % st_section(mask_fade) = fade_values(mask_fade);
        % st_section(mask_zero) = 0;
        
        % get st_downslope as min sample
        st_downslope = min(downslope, st_section);
        win = hamming(50)';
        win = win / sum(win);
        st_downslope = conv2(st_downslope, win, 'same');

        % combine all segments
        S = upslope .* st_downslope;

        % normalize result S
        S = bsxfun(@rdivide,S,max(S,[],2));
        % S = (S-min(S))./(max(S)-min(S));
    end

    if mode==getsMode.Exp_Spline
        TDEP=T-dep*ones(1,nt);

        % apply the up-slope: logistic shape
        upslope = (1-p(2))./(1+exp(-p(1)*TDEP)); 
        upslope = bsxfun(@rdivide,upslope,max(upslope,[],2)); % normalize

        % --------------------------- st_section --------------------------
        % get "first visible" highest point as ST begin
        st_begin_idx = first_crossing_up_thr(upslope, 0.999);

        % create ST section template 
        t = 1:nt; % * p(3);
        delta_t = max(0, t - st_begin_idx(:));
        st_section = max(0, 1 + (p(3) - 1) * delta_t / (nt));
        % st_section = max(0, 1 - max(0, (t - st_begin_idx(:))/(nt * 5 * p(3))));

        % spline downslope
        x_stretch = ceil(50 + 200*p(5));
        y_stretch = 0.1 * abs(1 - p(5));

        rep_int = ceil(rep);
        x_start_base = ceil(rep_int - 0.5 * x_stretch);
        x_start_actual = round((1 - p(4)) * x_start_base + p(4) * st_begin_idx(:));
        MOD_rest = [0, 0.5, 1];
        x_rep = [x_start_actual, ceil(rep_int + x_stretch * MOD_rest)];

        st_t_joint_idx = sub2ind(size(st_section), (1:nn)', (x_rep(:, 1)));
        st_t_joint_val = st_section(st_t_joint_idx);
        y_rep = [st_t_joint_val, repmat([0.5, 0.15-y_stretch, 0], nn,1)];

        idx = sub2ind(size(st_section), (1:size(st_section,1))', rep_int + x_stretch);
        idx_prev = sub2ind(size(st_section), (1:size(st_section,1))', rep_int + x_stretch - 1);
        dv_start = st_section(idx) - st_section(idx_prev);
        dv_end = zeros(size(rep));

        y_clamped = [dv_start, y_rep, dv_end];

        st_downslope = ones(nn, nt);
        for i = 1:nn      
            idx_start = x_rep(i,1);
            idx_end   = x_rep(i,end);
            N = idx_end - idx_start + 1;

            pp = spline(x_rep(i,:), y_clamped(i,:));
            st_downslope(i, 1:idx_start-1) = st_section(i, 1:idx_start-1);
            st_downslope(i, idx_start:idx_end) = ppval(pp, linspace(idx_start, idx_end, N));
            st_downslope(i, idx_end+1:end) = 0;
        end

        % add spike
        st_spike = ones(nn, nt);
        fall_time = 20;
        fall_value = 0.9;
        for i = 1:nn
            start_idx = st_begin_idx(i);
            st_spike(i, start_idx : start_idx + fall_time - 1) = linspace(1, fall_value, fall_time);
            st_spike(i, start_idx + fall_time : end) = fall_value;
        end

        % combine all segments
        S = upslope .* st_downslope .* st_spike;

        % normalize result S
        S = bsxfun(@rdivide,S,max(S,[],2));
        % S = (S-min(S))./(max(S)-min(S));

    end

    if mode==getsMode.UpslopeDownslope
        TREP=T-rep*ones(1,nt)-p(5);
        % compute (-1*) derivative of the downward slope of the TMP;
    
        Y=(p(2)+1./(1+exp(p(3)*TREP)))./(1+exp(p(4)*TREP));
    
        % apply the up-slope: logistic shape;
        S=(1-Y)./(1+exp(-p(1)*TDEP));
    
        % force a unit upstroke
        S=bsxfun(@rdivide,S,max(S,[],2));   
    end

end

% -------------------------------- HELPERS --------------------------------

function idx = first_crossing_up_thr(A, pct)
    if nargin < 2, pct = 0.99; end
    idx = first_crossing_thr(A, pct, true);
end

function idx = first_crossing_down_thr(A, pct)
    if nargin < 2, pct = 0.99; end
    idx = first_crossing_thr(A, pct, false);
end

function idx = first_samples_diff(A, diff_thr)
    if (isvector(A))
        dA = diff(A);
        for idx = 1:length(dA)
            if(abs(dA(idx)) > diff_thr)
                break
            end
        end
        return
    end

    dA = diff(A, 1, 2);
    mask = abs(dA) > diff_thr;
    [~, idx] = max(mask, [], 2);
    idx(~any(mask, 2)) = NaN;
end


function [p5, p6] = get_dependent_xy_values(x)
    l = 250;
    p5 = sin(pi/4 + x*pi/4);       % Y
    p6 = round(cos(x*pi/4) * l);   % X
end