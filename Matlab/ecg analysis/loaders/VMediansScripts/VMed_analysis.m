close all

steepest_tdom = [0 -1];
flattest_tdom = [-1 -1];

loader = TwaveLoader();
figure(1);
for i = 1:loader.data_cnt
    [twave, tdom] = loader.getVMediansTwave(i);

    subplot(1, 2, 1);
    y = twave;
    % y = conv(y, ones(1, 25), "same");
    y = baselinecor(y);
    y = y./max(y);
    [~, max_idx] = max(y);
    x1 = (-max_idx:length(y)-max_idx-1);
    plot( ...
        x1, y, ...
        'HandleVisibility', 'off', ...
        'Color', [0 0 0 0.2] ...
        );
    hold on

    subplot(1, 2, 2)
    [min_val, min_idx] = min(diff(tdom));
    x2 = (-min_idx:length(tdom)-min_idx-1);
    plot( ...
        x2, tdom, ...
        'HandleVisibility', 'off', ...
        'Color', [0 0 0 0.2] ...
        );
    hold on

    if(min_val < steepest_tdom(1))
        steepest_tdom = [min_val, i];
    end

    if(min_val > flattest_tdom(1))
        flattest_tdom = [min_val, i];
    end

    % if(i == 19)
    %     plot( ...
    %     x2, tdom, ...
    %     'Color', 'r', ...
    %     'LineWidth', 3, ...
    %     'DisplayName',sprintf("%d", i) ...
    %     );
    % hold on
    % end
end

subplot(1, 2, 1);
title('RMS of Twave');
grid on
subplot(1, 2, 2);
title('Tdom');
grid on

%% the steepest downslope
if(steepest_tdom(2) ~= -1)
    [twave, tdom] = loader.getVMediansTwave(steepest_tdom(2));
    [name, twave_begin, twave_end] = loader.getVMediansPathAndRanges(steepest_tdom(2));
    
    subplot(1, 2, 1);
    y = twave;
    % y = conv(y, ones(1, 25), "same");
    y = baselinecor(y);
    y = y./max(y);
    [~, max_idx] = max(y);
    x1 = (-max_idx:length(y)-max_idx-1);
    plot( ...
        x1, y, ...
        'Color', 'r', ...
        'LineWidth', 3, ...
        'DisplayName', sprintf('The steepest Tdom\n%s [%d %d]', name, twave_begin, twave_end) ...
        );
    hold on
    legend()
    subplot(1, 2, 2)
    [min_val, min_idx] = min(diff(tdom));
    x2 = (-min_idx:length(tdom)-min_idx-1);
    plot( ...
        x2, tdom, ...
        'Color', 'r', ...
        'LineWidth', 3, ...
        'DisplayName', sprintf('The steepest Tdom\n%s [%d %d]', name, twave_begin, twave_end) ...
        );
    hold on
    legend()
end

%% the flattest downslope
if(flattest_tdom(2) ~= -1)
    [twave, tdom] = loader.getVMediansTwave(flattest_tdom(2));
    
    subplot(1, 2, 1);
    y = twave;
    % y = conv(y, ones(1, 25), "same");
    y = baselinecor(y);
    y = y./max(y);
    [~, max_idx] = max(y);
    x1 = (-max_idx:length(y)-max_idx-1);
    plot( ...
        x1, y, ...
        'Color', 'g', ...
        'LineWidth', 3, ...
        'DisplayName', sprintf('The flattest Tdom - %d', flattest_tdom(2)) ...
        );
    hold on
    legend()
    subplot(1, 2, 2)
    [min_val, min_idx] = min(diff(tdom));
    x2 = (-min_idx:length(tdom)-min_idx-1);
    plot( ...
        x2, tdom, ...
        'Color', 'g', ...
        'LineWidth', 3, ...
        'DisplayName', sprintf('The flattest Tdom - %d', flattest_tdom(2)) ...
        );
    hold on
    legend()
end





