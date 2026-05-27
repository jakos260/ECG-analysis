close all

loader = TwaveLoader();
fitted_theta = zeros(loader.data_cnt, 3);
fitted_stats = zeros(loader.data_cnt, 3);
figure(1);

for i = 1:loader.data_cnt
    [name, twave_begin, twave_end] = loader.getVMediansPathAndRanges(i);
    [twave, tdom] = loader.getVMediansTwave(i);
    printProgress(i, loader.data_cnt, name);

    ecg = loadmat(name);
    signal = rms([ecg(4:6,:)/1.5; ecg(7:end,:)]);

    [~, l] = size(tdom);
    y = zeros(1, 500);
    y(1:l) = tdom;

    [theta_best, y_best, err, iter] = twave_optimize_quamin(y);
    fitted_theta(i,1:3) = theta_best;
    fitted_stats(i,1) = err;
    fitted_stats(i,2) = iter;
    % fprintf("Best theta = [%.3f, %.3f, %.3f] \nError = %f\nIterations proceed = %d\n", theta_best(1), theta_best(2), theta_best(3), err, iter);
    
    subplot(1, 2, 1);
    y = twave;
    % y = conv(y, ones(1, 25), "same");
    y = baselinecor(y);
    % y = y./max(y);
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
end
theta_opt = [
    mean(fitted_theta(:, 1))
    mean(fitted_theta(:, 2))
    mean(fitted_theta(:, 3))
    ];
y_fit = twave_generator(theta_opt, (1:500));
[~, max_idx] = max(abs(diff(y_fit)));
x1 = (-max_idx:length(y_fit)-max_idx-1);
plot( ...
    x1, y_fit, ...
    'Color','r', ...
    'LineWidth',3, ...
    'DisplayName',sprintf("fitted err=%.2f", mean(fitted_stats(:, 1))) ...
    );

hold on


subplot(1, 2, 1);
title('RMS of Twave');
grid on
legend()

y2 = s2t(y_fit);
[max_val, max_idx] = max(y2);
x2 = (-max_idx:length(y2)-max_idx-1);
plot( ...
    x2, y2, ...
    'Color','r', ...
    'LineWidth',3, ...
    'DisplayName',sprintf("fitted Twave", mean(fitted_stats(:, 1))) ...
    );

plot(x2, y2*20, 'Color','b', 'LineWidth',3,'DisplayName','fitted Twave x20');
subplot(1, 2, 2);
title('Tdom');
grid on
legend()

disp(sprintf("mean\n" + ...
    "theta 1  = %f\n" + ...
    "theta 2  = %f\n" + ...
    "theta 3  = %f\n" + ...
    "mean err = %f\n", ...
    mean(fitted_theta(:, 1)), ...
    mean(fitted_theta(:, 2)), ...
    mean(fitted_theta(:, 3)), ...
    mean(fitted_stats(:, 1)) ...
    ));

% 
% %%
% figure
% hold on
% grid on
% 
% histogram(fitted_theta(:, 1), 15, 'FaceColor','r', 'FaceAlpha',0.5);
% histogram(fitted_theta(:, 2), 15, 'FaceColor','g', 'FaceAlpha',0.5);
% histogram(fitted_theta(:, 3), 15, 'FaceColor','b', 'FaceAlpha',0.5);
% % histogram(fitted_stats(:, 1), 'FaceColor','k', 'FaceAlpha',0.5);
% legend("theta1", "theta2", "theta3", "iter");
% hold off

%%
close all

loader = TwaveLoader();
fitted_theta = zeros(loader.data_cnt, 3);
fitted_stats = zeros(loader.data_cnt, 3);
% figure(1);

for i = 1:loader.data_cnt
    [name, twave_begin, twave_end] = loader.getVMediansPathAndRanges(i);
    [twave, tdom] = loader.getVMediansTwave(i);
    signal = loader.getVMediansOriginalSignalsRms(i);
    printProgress(i, loader.data_cnt, name);
    figure('Name', name, 'NumberTitle', 'off');

    [~, l] = size(tdom);
    y = zeros(1, 500);
    y(1:l) = tdom;

    [theta_best, y_best, err, iter] = twave_optimize_quamin(y);
    fitted_theta(i,1:3) = theta_best;
    fitted_stats(i,1) = err;
    fitted_stats(i,2) = iter;
    % fprintf("Best theta = [%.3f, %.3f, %.3f] \nError = %f\nIterations proceed = %d\n", theta_best(1), theta_best(2), theta_best(3), err, iter);
    
    subplot(1, 2, 1);
    y = twave;
    % y = baselinecor(y);
    [max_val_y, max_idx_y] = max(y);
    [min_val_y, min_idx_y] = min(y);
    x1 = (-max_idx_y:length(y)-max_idx_y-1);
    x2 = (-max_idx_y-twave_begin:length(signal)-max_idx_y-twave_begin-1);
    plot(x2, signal, 'Color', [0.5 0.5 0.5 1], 'HandleVisibility', 'off');
    hold on
    plot( ...
        x1, y, ...
        'HandleVisibility', 'off', ...
        'LineWidth',3, ...
        'Color', 'k' ...
        );
    hold on
    y2 = s2t(y_fit);
    [max_val, max_idx] = max(y2);
    [min_val, min_idx] = min(y2);
    x2 = (-max_idx:length(y2)-max_idx-1);
    plot( ...
        x2, y2, ...
        'Color','r', ...
        'LineWidth',2, ...
        'DisplayName',sprintf("fitted Twave") ...
        );
    y2_norm = (y2-min_val)/(max_val-min_val)*((max_val_y-min_val_y)+min_val_y);
    % err = 
    plot(x2, y2_norm, 'Color','g', 'LineWidth',2,'DisplayName','fitted Twave norm');
    title('RMS of Twave');
    grid on
    legend()

    subplot(1, 2, 2)
    [min_val, min_idx] = min(diff(tdom));
    x2 = (-min_idx:length(tdom)-min_idx-1);
    plot( ...
        x2, tdom, ...
        'HandleVisibility', 'off', ...
        'Color', 'k' ...
        );
    hold on
    y_fit = twave_generator(theta_best, (1:500));
    [~, max_idx] = max(abs(diff(y_fit)));
    x1 = (-max_idx:length(y_fit)-max_idx-1);
    plot( ...
        x1, y_fit, ...
        'Color','r', ...
        'LineWidth',3, ...
        'DisplayName',sprintf("fitted err=%.2f\nt1 = %.2f\nt2 = %.2f", err, theta_best(1), theta_best(2)) ...
        );
    
    hold on

    title('Tdom');
    grid on
    legend()
end




figure
hold on
grid on

histogram(fitted_theta(:, 1), 15, 'FaceColor','r', 'FaceAlpha',0.5);
histogram(fitted_theta(:, 2), 15, 'FaceColor','g', 'FaceAlpha',0.5);
histogram(fitted_theta(:, 3), 15, 'FaceColor','b', 'FaceAlpha',0.5);
% histogram(fitted_stats(:, 1), 'FaceColor','k', 'FaceAlpha',0.5);
legend("theta1", "theta2", "theta3", "iter");
hold off