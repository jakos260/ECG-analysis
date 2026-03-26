close all

% get data to BSM simulation
V = loadmat("data/ExportData/normal_male/model/ventricles2Thorax.mat");
dst = loadmat('data/ExportData/normal_male/model/ventricle.dst3d');      % distances (mm)
A64 = loadmat('data/ExportData/normal_male/model/ventricles2BSM_(nijmegen_64).mat');

% get reference BSM (nijmegen_64)
BSM_ref = loadmat('data/ExportData/normal_male/ecgs/BSM_(nijmegen_64).refECG');


v = 1.5;          % m/s
dst_mm = dst;     % ECGSIM stores mm
time_ms = (dst_mm / 1000) / v * 1000;
root_nodes = [20];

% dep compute (Dijkstra)
time_ms(time_ms == 0) = Inf;

G = graph(time_ms,'upper');  % undirected graph

dep = distances(G, root_nodes);
dep = min(dep,[],1)';        % earliest arrival
dep = dep - min(dep);        % start at 0 ms

disp([min(dep), max(dep)]);

rep = 350 + mean(dep) - dep*0.6;

% time vector setup
start = 700;
ending = 1800;
maxt = ending-start;
maxt = 599;
T = ones(length(dep),1)*(0:maxt);
p = [-0.025, 0.2];

%%
S = twave_generator_full(T, dep, rep, p);
potentials = V * S; % potentials map [torso_nodes x signals]


% Forward projection
BSM_sim = A64 * S;          % [64 x T] true BSM

BSM_sim = BSM_sim - mean(BSM_sim,2);
BSM_ref = BSM_ref - mean(BSM_ref,2);

corr(BSM_sim, BSM_ref);

%%
close all

num_to_plot = 23;
figure(2)
sgtitle(sprintf('BSM(%d)', num_to_plot));

subplot(2,1,1);
for i = 1:64
    if(i ~= num_to_plot)
        plot(S(i,:), 'Color', [0 0 0 0.2]);
    end
    hold on
end

plot(S(i,:), 'r');
hold on
ylabel('mV');
grid on

subplot(2,1,2);
% plot(potentials(num_to_plot,:), 'r', 'DisplayName', 'potentials');
% hold on
plot(BSM_ref(num_to_plot,:), 'b', 'DisplayName', 'BSM\_ref');
hold on
plot(BSM_sim(num_to_plot,:),'k', 'DisplayName', 'BSM\_sim');
hold on
ylabel('mV');
xlabel('Time (s)')
grid on
legend

% plot(t,V2*1000)
% title('V2')
% ylabel('mV')
% xlabel('Time (s)')

%% qtriplot

[LCAV_V, LCAV_T] = read_ecgsim_geometry('data/ExportData/normal_male/model/lcav.tri');
[RCAV_V, RCAV_T] = read_ecgsim_geometry('data/ExportData/normal_male/model/rcav.tri');
[LLUNG_V, LLUNG_T] = read_ecgsim_geometry('data/ExportData/normal_male/model/llung.tri');
[RLUNG_V, RLUNG_T] = read_ecgsim_geometry('data/ExportData/normal_male/model/rlung.tri');
[THORAX_V, THORAX_T] = read_ecgsim_geometry('data/ExportData/normal_male/model/thorax.tri');

q = initQtripy();
q.surface(LCAV_V, LCAV_T)
q.surface(RCAV_V, RCAV_T)
q.surface(LLUNG_V, LLUNG_T)
q.transparency(0.5)
q.surface(RLUNG_V, RLUNG_T)
q.transparency(0.5)
q.surface(THORAX_V, THORAX_T)
q.transparency(0.7)
q.edge('w')
q.property_on_mouse_click('ver')


%%
figure
scatter3(1:length(dep), dep, dep*0, 40, dep, 'filled')
colorbar
xlabel('Node index')
ylabel('Activation time (ms)')
title('Ventricular activation times')

%%
figure
plot(dep, (1:length(dep)), '.')
ylabel('Ventricular node')
xlabel('Activation time (ms)')
title('Depolarization times')


