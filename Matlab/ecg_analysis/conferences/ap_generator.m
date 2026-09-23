clf;
CT = 1;

type_list = {'Epicardium', 'Mid-myocardial', 'Endocardium'};
type = type_list{CT};
x_labels = {'Phase 1', 'Phase 2', 'Phase 3'};
params = [1,1,1];


HT = 0.02;
STOPTIME = 380;
LIM = 50;

figure(23);

t = linspace(HT,STOPTIME,STOPTIME/HT);
hold on;

[t_tmpl, V_tmpl] = wrapper_TenTusscher2mod(HT, STOPTIME, CT, params);
V_min = min(V_tmpl);
V_max = max(V_tmpl);
x_size = @(x) (t_tmpl(13297) - t_tmpl(101)) * x;
y_size = @(y) (V_max - (V_max - V_min) * y);

plot([t_tmpl(101), t_tmpl(13297)],[-30,-30], 'g', 'LineWidth', 2.5) % APD
text(x_size(0.5), y_size(0.5), 'APD', ...
    'HorizontalAlignment', 'center', ...
    'VerticalAlignment', 'bottom', ...
    'Color', 'g', ...
    'FontWeight', 'bold', ...
    'FontSize', 12); % APD
plot([t_tmpl(101), t_tmpl(101)],[-30,-90], 'b--', 'LineWidth', 1) % DEP
text(x_size(0.1), y_size(1.05), 'DEP time', ...
    'HorizontalAlignment', 'center', ...
    'VerticalAlignment', 'bottom', 'Color', 'b', 'FontWeight', 'bold', 'FontSize', 12); % DEP

plot([t_tmpl(13297), t_tmpl(13297)],[-30,-90], 'r--', 'LineWidth', 1) % REP
text(x_size(1.10), y_size(1.05), 'REP time', ...
    'HorizontalAlignment', 'center', ...
    'VerticalAlignment', 'bottom', 'Color', 'r', 'FontWeight', 'bold', 'FontSize', 12); % REP

plot(t_tmpl, V_tmpl, 'k', 'LineWidth', 2.5);
ylabel('Amplitude [mV]', 'FontSize', 15);
xlabel('Time [ms]', 'FontSize', 15);

title('AP phases', 'FontWeight', 'bold', 'FontSize', 20);
xlim([-LIM, STOPTIME+LIM]);

text(x_size(-0.08), y_size(0.60), sprintf('phase 0'), 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', 'FontWeight', 'bold', 'FontSize', 15);
text(x_size(0.10), y_size(0.25), sprintf('phase 1\nI_{to}'), 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', 'FontWeight', 'bold', 'FontSize', 15);
text(x_size(0.50), y_size(0.08), sprintf('phase 2\nI_{Ca}'), 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', 'FontWeight', 'bold', 'FontSize', 15);
text(x_size(1.15), y_size(0.60), sprintf('phase 3\nI_{Kr} + I_{Ks}'), 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', 'FontWeight', 'bold', 'FontSize', 15);
text(x_size(1.30), y_size(0.95), sprintf('phase 4'), 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', 'FontWeight', 'bold', 'FontSize', 15);