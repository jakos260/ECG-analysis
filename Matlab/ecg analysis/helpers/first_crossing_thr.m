function idx = first_crossing_thr(A, pct, direction_up)

    mx  = max(A, [], 2);       % max per row
    thr = pct .* mx;           % threshold per row

    % Shifted matrices for edge detection
    A1 = A(:,1:end-1);
    A2 = A(:,2:end);

    if direction_up
        cross = (A1 < thr) & (A2 >= thr);
    else
        cross = (A1 > thr) & (A2 <= thr);
    end

    % Find first crossing
    any_cross = any(cross, 2);
    idx = NaN(size(thr));

    idx(any_cross) = sum(cumprod(~cross(any_cross,:), 2), 2) + 1;

    % Shift by +1 because crossing occurs between samples
    idx(any_cross) = idx(any_cross) + 1;
end
