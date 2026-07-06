% Simulated SCO frequencies for a toy model over a range of instrinic frequencies

runtime = 10; % Duration of simulation in seconds
dt = 0.001; % Timestep in seconds
tsteps = 0:dt:runtime;
N = 1000;

% Intrinsic frequenices
f_A = 10.5:.05:15; % Range of intrinsic frequencies for module A
f_B = 10;
f_diff = f_A - f_B;

% Local and global coupling
k_local_A = .1;
k_local_B = .1;
k_global = 0;
G = 1;

% Run simulation
r_a = zeros(length(f_A),length(tsteps));
r_b = zeros(length(f_A),length(tsteps));
r = zeros(length(f_A),length(tsteps));

if isempty(gcp('nocreate')); parpool('local'); end % Initialise parpool

parfor f=1:length(f_A)

    Omega = zeros(N,1);
    Omega(1:N/2) = f_A(f);
    Omega(N/2 + 1:N) = f_B;
    Omega = Omega * 2*pi;

  
    Conn = zeros(N,N);
    Conn(1:N/2,1:N/2) = k_local_A*ones(N/2,N/2);
    Conn(1:N/2,N/2+1:N) = k_global*ones(N/2,N/2);
    Conn(N/2+1:N,1:N/2) = k_global*ones(N/2,N/2);
    Conn(N/2+1:N,N/2+1:N) = k_local_B*ones(N/2,N/2);
    Conn = Conn - diag(diag(Conn));  % Diagonal elements are zeros
    
    theta = 1000*randn(length(tsteps),N);
    
    delay = 5;
    
    for tt = delay+2:length(tsteps)
        
            theta(tt,:) = (Omega' + G*sum(Conn.*sin(repmat(theta(tt-1-delay,:)',[1 N]) - theta(tt-1,:)),1))*dt + theta(tt-1,:);
        
    end
    
    % Log the coherence timeseries for the current intrinsic frequency
    r_a(f,:) = abs(mean(exp(1i*(theta(:,1:N/2))),2));
    r_b(f,:) = abs(mean(exp(1i*(theta(:,N/2+1:N))),2));
    r(f,:) =  abs(mean(exp(1i*theta),2));

end

save('./figure_data/toy_hfreq_results.mat','f_A','r','f_diff','f_B','k_local_B','k_local_A','N')  