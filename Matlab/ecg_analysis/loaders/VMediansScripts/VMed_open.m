
dirName = 'data/VMedians/';
fileType = '*.vmedianecg';
files = dir(append(dirName, fileType));
usedPrefixes = {};

for i=1:size(files, 1)
    fileName = files(i, 1).name;
    prefix = fileName(1:4);

    printProgress(i, size(files, 1), sprintf("%s %d uniqe ecg add", prefix, size(usedPrefixes, 2)));
    if ismember(prefix, usedPrefixes)
        continue
    end
    usedPrefixes{end+1} = prefix;

    ecg = loadmat(append(dirName, fileName));

    figure('Name', fileName, 'NumberTitle', 'off');
    leadv12(ecg);
    hold on
    % if (i == 20); break; end
end

%%
close all
data_list = dir("data/VMedians/*.vmedianecg");

% ?? 13, 15, 32, 64
% N = N+1;
N = 44;
name = append('data/VMedians/', data_list(N).name);
ecg = loadmat(name);
plot(ecg(:, :)');
hold on
plot(rms(ecg(:, :)), "LineWidth",2,"Color","k");
hold on
grid("on");
title(sprintf("%d - %s\n", N, name));
fprintf("%d - %s\n", N, name);

%%
close all

number = 44;
loader = TwaveLoader();
[name, ~, ~] = loader.getVMediansPathAndRanges(number);

ecg = loadmat(name);
plot(ecg(:, :)');
hold on
plot(rms(ecg(:, :)), "LineWidth",2,"Color","k");
hold on
grid("on");
title(sprintf("%d - %s\n", number, name));
fprintf("%d - %s\n", number, name);

%%
close all

number = 19;
loader = TwaveLoader();

[name, twave_begin, twave_end] = loader.getVMediansPathAndRanges(number);
[twave, tdom] = loader.getVMediansTwave(number);
ecg = loadmat(name);

plot(rms(ecg), 'k');
hold on
plot(tdom, 'r');
title(sprintf("%d - %s", number, name));
%%
close all

loader = TwaveLoader();
figure(1);
for i = 1:loader.data_cnt
    [twave, tdom] = loader.getVMediansTwave(i);

    subplot(1, 2, 1);
    y = conv(rms(twave), ones(1, 25), "same");
    y = y./max(y);
    [~, max_idx] = max(y);
    x1 = (-max_idx:length(y)-max_idx-1);
    plot( ...
        x1, y, ...
        'Color', [0 0 0 0.2], ...
        'DisplayName',sprintf("%d", i) ...
        );
    hold on

    subplot(1, 2, 2)
    [~, min_idx] = min(diff(tdom));
    x2 = (-min_idx:length(tdom)-min_idx-1);
    plot( ...
        x2, tdom, ...
        'Color', [0 0 0 0.2], ...
        'DisplayName',sprintf("%d", i) ...
        );
    hold on

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

%%
close all


loader = TwaveLoader();
for i = 1:loader.data_cnt
    % if ~ismember(i, [3 15 19 25 34 43])
    %     continue
    % end

    [name, twave_begin, twave_end] = loader.getVMediansPathAndRanges(i);
    [twave, tdom] = loader.getVMediansTwave(i);

    ecg = loadmat(name);
    signal = rms([ecg(4:6,:)/1.5; ecg(7:end,:)]);

    figure('Name', sprintf("%d - %s", i, name), 'NumberTitle', 'off');
    plot(signal, "k");
    hold on
    x = (twave_begin : twave_end);
    plot( ...
        x, signal(1, twave_begin : twave_end), ...
        "LineWidth", 3, ...
        "Color", "r", ...
        "DisplayName", sprintf("twave range = (%d : %d)", twave_begin, twave_end) ...
        );
    hold on;
    plot( ...
        x, baselinecor(twave./max(twave)), ...
        "LineWidth", 1, ...
        "Color", "b", ...
        "DisplayName", "baselinecor(twave/max(twave))" ...
        );
    legend;
    grid;

end