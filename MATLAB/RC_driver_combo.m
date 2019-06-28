clear; close all; clc;

tic

fileName = datestr(datetime('now'));

num_sim = 10; % number of simulations #sim
thold_set = [0.2, 0.3, 0.4, 0.5]; % set of input hold times
Sm_base_set = [1e-4, 2e-4, 3e-4, 4e-4, 5e-4]; % set of base influx rates

% Pre-allocation
std_NRMSE_A = zeros(length(Sm_base_set), length(thold_set));
avg_NRMSE_A = zeros(length(Sm_base_set), length(thold_set));
std_NRMSE_B = zeros(length(Sm_base_set), length(thold_set));
avg_NRMSE_B = zeros(length(Sm_base_set), length(thold_set));

% Run experiments
for smbaseset = 1 : length(Sm_base_set)
    
    for tholdset = 1 : length(thold_set)
        
        thold = thold_set(tholdset);
        Sm_base = Sm_base_set(smbaseset);
        
        % Pre-allocation
        NRMSE_A_1set = zeros(num_sim, 1);
        NRMSE_B_1set = zeros(num_sim, 1);
        
        % Run #sim experiments
        for numsim = 1 : num_sim
            
            % Task A
            TaskA_Learning;
            NRMSE_A_1set(numsim) = NRMSE_A;
            std_NRMSE_A(smbaseset,tholdset) = std(NRMSE_A_1set);
            avg_NRMSE_A(smbaseset,tholdset) = mean(NRMSE_A_1set);
            
            % Task B
            TaskB_Learning;
            NRMSE_B_1set(numsim) = NRMSE_B;
            std_NRMSE_B(smbaseset,tholdset) = std(NRMSE_B_1set);
            avg_NRMSE_B(smbaseset,tholdset) = mean(NRMSE_B_1set);
            
        end
        
    end
    
end

% Save useful data
save(fileName, 'std_NRMSE_A', 'avg_NRMSE_A', 'std_NRMSE_B', 'avg_NRMSE_B');

toc