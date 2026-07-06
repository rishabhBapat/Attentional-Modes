% Model with all to all connectivity and Lorentzian frequency distribution
if isempty(gcp('nocreate')); parpool('local'); end

N = 1000;
conn = ones(N,N);

t = 100; % Duration of simulation in seconds
dt = 0.001; % Simulation step in seconds
tspan = 0:dt:t; % Time steps
t_down = length(downsample(tspan,100));

delta = .02;
mu = 10; % Location parameter

Omega = mu + delta*tan(pi*(rand(N,1)-1/2)); % Cauchy distributed natural frequencies
Omega_ = Omega*2*pi;

G = 0:.000012:.001;

delay = 0:1:15;
mtsb = zeros(length(G),length(delay));
ord_param = zeros(length(G),length(delay),t_down);

parfor g = 1:length(G)
    temp_mtsb = zeros(length(delay),1);
    r_down = zeros(length(delay),t_down);
    for d = 1:length(delay)
        Theta = 1000*randn(length(tspan),N);

        % Phase timeseries for current value of delay and coupling
        for tt = delay(d)+1:length(tspan)-1
            Theta(tt+1,:) = (Omega_' + G(g)*sum(conn.*sin(repmat(Theta(tt-delay(d),:)',[1 N]) - Theta(tt,:)),1))*dt + Theta(tt,:);
        end
        
        r = abs(mean(exp(1i*Theta),2));
        temp_mtsb(d) = std(r(tspan>10));
        r_down(d,:) = downsample(r,100);
    end
    mtsb(g,:) = temp_mtsb;
    ord_param(g,:,:) = r_down;
end

save('./figure_data/lorentz_results.mat', "mtsb", "G", "ord_param", "delay")