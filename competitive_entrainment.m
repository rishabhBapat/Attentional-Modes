% Time-frequency analysis of node alternating between stimulus entrainment and network

% Load dependencies

% Human
%Theta_ds = load('./dependencies/theta_human.mat').Theta_Down; % Load simulated phase timeseries
%LUT = readtable('./dependencies/LUT_schaefer.csv').ROIName; % Load lookup table
%tp_ind  = find(contains(LUT, 'TempPar'));

% Macaque
Theta = load('./dependencies/theta_macaque.mat').Theta;
Theta_ds = downsample(Theta,10);
LUT = readtable('./dependencies/LUT_shen.txt');
temporoparietalRegions = {'secondary auditory cortex', 'inferior parietal cortex', 'superior temporal cortex', 'central temporal cortex'};
tp_ind = find(ismember(LUT.name, temporoparietalRegions));

tspan_Down = 0:0.01:100;

f_0 = 5*2*pi; % Sensory node
f_s = 3*2*pi;  % Stimulus frequency
k_s = 100; % Stimulus coupling
k_n = 200; % Network coupling

% Run simulation

dt = .01;
n = length(tp_ind);
stim = f_s*tspan_Down;
p = randn(length(tspan_Down),1);

for t = 1:length(tspan_Down)-1
    p(t+1) = (f_0 + (k_n/n)* sum(sin(Theta_ds(t,tp_ind) - p(t))) + k_s*sin(stim(t)-p(t)))*dt + p(t);
end

[wt,f] = cwt(cos(p),1/dt);

% Save results

% Human
%save("./figure_data/entrainment_results_human.mat", 'wt', 'f')

% Macaque
%save("./figure_data/entrainment_results_macaque.mat", 'wt', 'f')