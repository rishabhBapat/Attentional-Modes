% Simulated SCO frequencies for a toy model acoss local coupling

runtime = 10; % Duration of simulation in seconds
dt = 0.001; % Timestep in seconds
tsteps = 0:dt:runtime;
N = 1000;
Omega = zeros(N,1);

% Intrinsic frequenices
Omega(1:N/2) = 10;
Omega(N/2 + 1:N) = 10;
Omega = Omega * 2*pi;

delay = 5;

% Local and global coupling
k_A = .2:.01:1; % Range of local coupling values for module A
k_local_B = .1;
k_global = 0;
G = 1;
k_diff = k_A - k_local_B;

% Run simulation
r_a = zeros(length(k_A),length(tsteps));
r_b = zeros(length(k_A),length(tsteps));
r = zeros(length(k_A),length(tsteps));

if isempty(gcp('nocreate')); parpool('local'); end % Initialise parpool

parfor k=1:length(k_A)

    k_local_A = k_A(k);
    Conn = zeros(N,N);
    Conn(1:N/2,1:N/2) = k_local_A*ones(N/2,N/2);
    Conn(1:N/2,N/2+1:N) = k_global*ones(N/2,N/2);
    Conn(N/2+1:N,1:N/2) = k_global*ones(N/2,N/2);
    Conn(N/2+1:N,N/2+1:N) = k_local_B*ones(N/2,N/2);
    Conn = Conn - diag(diag(Conn)); % Diagonal elements are zeros
    
    theta = 1000*randn(length(tsteps),N);
    
    for tt = delay+2:length(tsteps)
            theta(tt,:) = (Omega' + G*sum(Conn.*sin(repmat(theta(tt-1-delay,:)',[1 N]) - theta(tt-1,:)),1))*dt + theta(tt-1,:);
    end
    
    % Log the coherence timeseries for the current delay
    r_a(k,:) = abs(mean(exp(1i*(theta(:,1:N/2))),2));
    r_b(k,:) = abs(mean(exp(1i*(theta(:,N/2+1:N))),2));
    r(k,:) =  abs(mean(exp(1i*theta),2));

end

save('./figure_data/toy_hcoupling_results.mat','k_A','r','k_diff','Omega','k_local_B')  