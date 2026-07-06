% Simulated phase timeseries across the coupling-delay parameter space

%Load dependencies

% AAL
C = load('./dependencies/sc_AAL.mat').C; 
LUT = load('./dependencies/LUT_AAL.mat').label90;

conn = C/max(C(:)); % Normalize connectivity
N = size(conn,1);

temporoparietalRegions = {'L Angular','R Angular', 'L SupraMarginal','R SupraMarginal', 'L Parietal Inf','R Parietal Inf', 'L Temporal Sup','R Temporal Sup'};
DMNRegions = {'L Front Sup Med','R Front Sup Med', 'L Cingulum Post','R Cingulum Post', 'L Precuneus','R Precuneus', 'L Angular','R Angular'};

tp_ind = find(ismember(LUT, temporoparietalRegions));
dmn_ind = find(ismember(LUT, DMNRegions));
all_ind = 1:N;
network_idx = {tp_ind;dmn_ind;all_ind};

t = 100; % Duration of simulation in seconds
dt = 0.001; % Timestep in seconds
tspan = 0:dt:t;
t_down = length(downsample(tspan,100));

% Intrinsic frequency range

% Human
high = 12; low = 8;
low = 12; % Use for homogenous intrinsic frequencies 

S_n = sum(conn, 2);
Omega = 2 * pi * (high - (high - low) * (((S_n - min(S_n))/(min(S_n) - max(S_n))).^2));

% Single run

G = 2;
delay = 4;

Theta = 1000*randn(length(tspan),N);
for tt = delay+1:length(tspan)-1
    Theta(tt+1,:) = (Omega' + G*sum(conn.*sin(repmat(Theta(tt-delay,:)',[1 N]) - Theta(tt,:)),1))*dt + Theta(tt,:);
end
Theta_ds = downsample(Theta,10);

%save('./dependencies/theta_AAL.mat', 'Theta_ds');

% Parameter Sweep
G = 0:.1:8; % Range of coupling values to sweep over
delay = 0:1:15; % Range of delay values to sweep over

mtsb = zeros(length(G),length(delay));
ord_param = zeros(length(G),length(delay),length(network_idx),t_down);

%parpool('local') % Initialise parpool 

parfor g = 1:length(G)
    temp_mtsb = zeros(length(delay),1);
    r_down = zeros(length(delay),length(network_idx),t_down);
    for d = 1:length(delay)
        Theta = 1000*randn(length(tspan),N);
        
        % Phase timeseries for current value of delay and coupling
        for tt = delay(d)+1:length(tspan)-1
            Theta(tt+1,:) = (Omega' + G(g)*sum(conn.*sin(repmat(Theta(tt-delay(d),:)',[1 N]) - Theta(tt,:)),1))*dt + Theta(tt,:);
        end
        
        % Coherence timeseries for each subnetwork for current value of delay and coupling
        for n = 1:length(network_idx)
            r_down(d,n,:) = downsample(abs(mean(exp(1i*Theta(:,network_idx{n})),2)),100);
        end

        r = abs(mean(exp(1i*Theta),2));
        temp_mtsb(d) = std(r(tspan>10));

    end
    mtsb(g,:) = temp_mtsb; % Log std(r) at current coupling for all delays
    ord_param(g,:,:,:) = r_down; % Log subnetwork SCOs at current coupling for all delays
end

% Save results

% Human
%save('./figure_data/sweep_hetero_AAL.mat','G','delay','mtsb','ord_param');
save('./figure_data/sweep_homo_AAL.mat','G','delay','mtsb','ord_param'); % Use for homogenous frequencies