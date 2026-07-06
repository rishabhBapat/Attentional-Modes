% Randomised human and macaque connectomes

% Load dependencies

C = load('./dependencies/sc_schaefer.mat').sc_schaefer; % Human
%C = load('./dependencies/sc_shen.txt'); % Macaque

C_norm = C/max(C(:)); % Normalize connectivity
N = size(C_norm,1);
C_norm = C_norm(reshape(randperm(numel(C_norm)),[N N])); % Randomize connectivity
C_norm = C_norm - diag(diag(C_norm));

% Intrinsic frequency range

high = 12; low = 8; % Human
%high = 11; low = 9; % Macaque

S_n = sum(C_norm, 2);
Omega = 2 * pi * (high - (high - low) * (((S_n - min(S_n))/(min(S_n) - max(S_n))).^2));

t = 100; % Duration of simulation in seconds
dt = 0.001; % Simulation step in seconds
tspan = 0:dt:t; % Time steps
t_down = length(downsample(tspan,100));

G = 0:.1:8;

delay = 0:1:15;
mtsb = zeros(length(G),length(delay));
ord_param = zeros(length(G),length(delay),t_down);

if isempty(gcp('nocreate')); parpool('local'); end % initialise parpool

parfor g = 1:length(G)
    temp_mtsb = zeros(length(delay),1);
    r_down = zeros(length(delay),t_down);
    for d = 1:length(delay)
        Theta = 1000*randn(length(tspan),N);

        % Phase timeseries for current value of delay and coupling
        for tt = delay(d)+1:length(tspan)-1
            Theta(tt+1,:) = (Omega' + G(g)*sum(C_norm.*sin(repmat(Theta(tt-delay(d),:)',[1 N]) - Theta(tt,:)),1))*dt + Theta(tt,:);
        end

        r = abs(mean(exp(1i*Theta),2));
        temp_mtsb(d) = std(r(tspan>10));
        r_down(d,:) = downsample(r,100);
    end
    mtsb(g,:) = temp_mtsb;
    ord_param(g,:,:) = r_down;

end

% Save results

save('./figure_data/randconn_results_human.mat', "mtsb", "G", "ord_param", "delay") % Human
%save('./figure_data/randconn_results_macaque.mat', "mtsb", "G", "ord_param", "delay") % Macaque