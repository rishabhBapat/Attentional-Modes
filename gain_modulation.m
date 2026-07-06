% Effect of multiplicative gain applied to either hemisphere on SCO frequency

% Load dependencies

% Human
%C = load('./dependencies/sc_schaefer.mat').sc_schaefer; % Load connectivity
%LUT = readtable('./dependencies/LUT_schaefer.csv').ROIName; % Load lookup table

% Macaque
C = readmatrix('./dependencies/sc_shen.txt'); 
LUT = readtable('./dependencies/LUT_shen.txt');

C_norm = C/max(C(:)); % Normalize connectivity

% Find indices of each network

% Human
%tp_ind  = find(contains(LUT, 'TempPar'));
%tp_lh   = find(contains(LUT, 'LH_TempPar'));
%tp_rh   = find(contains(LUT, 'RH_TempPar'));

% Macaque
temporoparietalRegions = {'secondary auditory cortex', 'inferior parietal cortex', 'superior temporal cortex', 'central temporal cortex'};
tp_ind = find(ismember(LUT.name, temporoparietalRegions));
tp_rh = find(ismember(LUT.name, temporoparietalRegions) & strcmp(LUT.hem, 'right'));
tp_lh = find(ismember(LUT.name, temporoparietalRegions) & strcmp(LUT.hem, 'left'));

target = tp_rh ; % Network to scale gain within

% Intrinsic frequency range

%high = 12; low = 8; % Human
high = 11; low = 9; % Macaque

S_n = sum(C_norm, 2);
Omega = 2 * pi * (high - (high - low) * (((S_n - min(S_n))/(min(S_n) - max(S_n))).^2));

N = size(C_norm, 1);
t = 1000; % Duration of simulation in seconds
dt = 0.001; % Timestep in seconds
nSteps = t / dt;

gains = linspace(1, 5, 10);
gains_repl = repelem(gains, 5);
n_gains = length(gains_repl);

if isempty(gcp('nocreate')); parpool('local'); end % initiate parpool

% CONDITION 1: delayed coupling

% G = 2; delay = 4; % Human
G = 2.6; delay = 5; % Macaque

buf_size = delay + 2;

gain_coh_peaks = zeros(1, n_gains);
parfor i = 1:n_gains
    rng(i, 'twister');
    theta_init = mod(1000 * randn(1, N), 2*pi);

    Theta_buf = repmat(theta_init, buf_size, 1);

    g_conn = C_norm;
    g_conn(target, target) = g_conn(target, target) * gains_repl(i);

    r_tp = zeros(nSteps + 1, 1);
    r_tp(1) = abs(mean(exp(1i * theta_init(tp_ind))));

    for tt = delay : nSteps - 1
    cur_idx = mod(tt, buf_size) + 1;
    del_idx = mod(tt - delay, buf_size) + 1;
    nxt_idx = mod(tt + 1, buf_size) + 1;

    coupling = G * sum(g_conn .* sin(repmat(Theta_buf(del_idx,:)', [1 N]) - Theta_buf(cur_idx,:)), 1);
    Theta_buf(nxt_idx,:) = mod((Omega' + coupling) * dt + Theta_buf(cur_idx,:), 2*pi);

    r_tp(tt + 2) = abs(mean(exp(1i * Theta_buf(nxt_idx, tp_ind))));
    end

    pks = findpeaks(-1 * r_tp, 'MinPeakProminence', 0.4);
    gain_coh_peaks(i) = length(pks) / t;
end

gain_coh = gain_coh_peaks;
mean_gain_peaks = mean(reshape(gain_coh_peaks, 5, []));
fprintf('Condition 1 done.\n');

% CONDITION 2: no delay
% G = 1.2; delay = 0; % Human
G = 2; delay = 0; % Macaque

buf_size = delay + 2;

gain_coh_peaks = zeros(1, n_gains);
parfor i = 1:n_gains
    rng(i, 'twister');
    theta_init = mod(1000 * randn(1, N), 2*pi);

    Theta_buf = repmat(theta_init, buf_size, 1);

    g_conn = C_norm;
    g_conn(target, target) = g_conn(target, target) * gains_repl(i);

    r_tp = zeros(nSteps + 1, 1);
    r_tp(1) = abs(mean(exp(1i * theta_init(tp_ind))));

    for tt = delay : nSteps - 1
    cur_idx = mod(tt, buf_size) + 1;
    del_idx = mod(tt - delay, buf_size) + 1;
    nxt_idx = mod(tt + 1, buf_size) + 1;

    coupling = G * sum(g_conn .* sin(repmat(Theta_buf(del_idx,:)', [1 N]) - Theta_buf(cur_idx,:)), 1);
    Theta_buf(nxt_idx,:) = mod((Omega' + coupling) * dt + Theta_buf(cur_idx,:), 2*pi);

    r_tp(tt + 2) = abs(mean(exp(1i * Theta_buf(nxt_idx, tp_ind))));
    end

    pks = findpeaks(-1 * r_tp, 'MinPeakProminence', 0.4);
    gain_coh_peaks(i) = length(pks) / t;
end

gain_coh_nd = gain_coh_peaks;
mean_gain_peaks_nd = mean(reshape(gain_coh_peaks, 5, []));
fprintf('Condition 2 done.\n');

% Save results

% Human
%if target == tp_rh
    %save('./figure_data/gainmod_results_human_right.mat', 'mean_gain_peaks', 'mean_gain_peaks_nd', 'gain_coh', 'gain_coh_nd');
%else
    %save('./figure_data/gainmod_results_human_left.mat', 'mean_gain_peaks', 'mean_gain_peaks_nd', 'gain_coh', 'gain_coh_nd');
%end

% Macaque
if target == tp_rh
    save('./figure_data/gainmod_results_macaque_right.mat', 'mean_gain_peaks','mean_gain_peaks_nd', 'gain_coh', 'gain_coh_nd');
else
    save('./figure_data/gainmod_results_macaque_left.mat', 'mean_gain_peaks', 'mean_gain_peaks_nd', 'gain_coh', 'gain_coh_nd');
end