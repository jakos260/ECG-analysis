clf;
CT = 1;

type_list = {'Epicardium', 'Mid-myocardial', 'Endocardium'};
type = type_list{CT};
x_labels = {'Phase 1', 'Phase 2', 'Phase 3'};
params = [1,1,1];


HT = 0.02;
STOPTIME = 500;
figure(23);

t = linspace(HT,STOPTIME,STOPTIME/HT);
hold on;

[t_tmpl, V_tmpl] = wrapper_TenTusscher2mod(HT, STOPTIME, CT, params);
plot([t_tmpl(101), t_tmpl(13297)],[-30,-30], 'r', 'LineWidth', 2.5) % APD
text(0.28, 0.52, 'APD', 'Units', 'normalized', 'Color', 'r', 'FontWeight', 'bold', 'FontSize', 12); % APD
plot(t_tmpl, V_tmpl, 'k', 'LineWidth', 2.5);
ylabel('Amplitude [mV]', 'FontSize', 15);
xlabel('Time [ms]', 'FontSize', 15);

title('AP phases', 'FontWeight', 'bold', 'FontSize', 20);
xlim([-50, 550]);

text(0.04, 0.40, sprintf('phase 0'), 'Units', 'normalized', 'HorizontalAlignment', 'center', 'FontWeight', 'bold', 'FontSize', 15);
text(0.13, 0.78, sprintf('phase 1\nI_{to}'), 'Units', 'normalized', 'HorizontalAlignment', 'center', 'FontWeight', 'bold', 'FontSize', 15);
text(0.30, 0.90, sprintf('phase 2\nI_{Ca}'), 'Units', 'normalized', 'HorizontalAlignment', 'center', 'FontWeight', 'bold', 'FontSize', 15);
text(0.60, 0.40, sprintf('phase 3\nI_{Kr} + I_{Ks}'), 'Units', 'normalized', 'HorizontalAlignment', 'center', 'FontWeight', 'bold', 'FontSize', 15);
text(0.80, 0.13, sprintf('phase 4'), 'Units', 'normalized', 'HorizontalAlignment', 'center', 'FontWeight', 'bold', 'FontSize', 15);