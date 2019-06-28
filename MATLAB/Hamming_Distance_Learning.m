clear; clc; close all;

tic
fileName = datestr(datetime('now'));
Sm_base = 0.0003;
num_sim = 100;
thold_set = 0.1:0.1:0.4;

std_NRMSE_H = zeros(length(thold_set), 1);
avg_NRMSE_H = zeros(length(thold_set), 1);

for tholdset = 1 : length(thold_set)
    
    thold = thold_set(tholdset);
    NRMSE_H_1set = zeros(num_sim, 1);
    
    for numsim = 1 : num_sim    
        
        % Run chemistry\
        [conS, S, freq, R_total, iniCon_ss, iniCon_ds, orderTable, Influx, ...
            timeVec, tmax, tPerturb, numPerturb, thold, perturbTimes, index1, index2, kI_base, nS] ...
            = random_DNA_Strand_circuit_ver5(Sm_base, thold);

        bit = numPerturb;

        % original inputs
        original_input1 = Influx(:, index1);
        original_input2 = Influx(:, index2);

        % rounded and normalized inputs
        input1 = round(original_input1 / (max(original_input1) + min(original_input1)));
        input2 = round(original_input2 / (max(original_input2) + min(original_input2)));

        % Calculate the target vector and the Hamming distance between the 2 inputs
        target = input1 ~= input2; 
        Hdist = sum(target);

        % Concentration dataset
        conS_2 = cell2mat(conS(2:end,:)); % Convert concentration dataset from cell type to matrix type
        conS_3 = zeros(size(conS_2, 1), nS); % Pre-allocation for the full concentration dataset
        conS_4 = zeros(bit, nS); % Pre-allocation for the scaled concentration dataset (to #bits)
        
        for iCon = 1 : nS
            conS_3(:,iCon) = conS_2(:,iCon) / (max(conS_2(:,1)) + min(conS_2(:,1))); % Normalized dataset
            conS_4(:,iCon) = conS_3(randi(size(conS_2, 1), [bit, 1])); % pick a random region of the dataset that has the same length as the number of bits
        end
        
        bias = 1; % Bias
        wBias = rand; wSpecies = rand(1,nS); % Weights
        coeff = 0.01; % Learning coefficient
        iterations = 10000; % Number of iterations

        % Preallocation
        actualOut_matrix = zeros(bit, iterations);
        error_matrix = zeros(bit, iterations);
        weights = [wSpecies, wBias];

        for i = 1 : iterations
            for j = 1 : bit
                net = sum(conS_4(j,:) .* wSpecies) + bias * wBias;    
                actualOut_matrix(j,i) = 1./(1+exp(-net));
                error_matrix(j,i) = target(j) - actualOut_matrix(j,i);

                wSpecies = wSpecies + coeff * conS_4(j,:) .* error_matrix(j,i); % species weight update
                wBias = wBias + coeff * bias * error_matrix(j,i); % bias weight update
                weights = [weights; wSpecies, wBias];
            end
        end
        error_matrix = abs(error_matrix);
        error = error_matrix(:,end);
        actualOut = actualOut_matrix(:,end);
        actual_Hdist = sum(actualOut);

        NRMSE_H = (1 / (max(actualOut) - min(actualOut))) * sqrt(sum(error.^2) / bit);
        
        NRMSE_H_1set(numsim) = NRMSE_H;
        std_NRMSE_H(tholdset) = std(NRMSE_H_1set);
        avg_NRMSE_H(tholdset) = mean(NRMSE_H_1set);
        
    end
    
end

% Save useful data
save(fileName, 'std_NRMSE_H', 'avg_NRMSE_H');

toc
    
% fprintf('\nResevoir computing: \n');
% fprintf(' Setup: \n');
% fprintf('  Number of iteration = %d \n', iterations);
% fprintf('  Learning rate = %g \n',  coeff);
% fprintf('  Number of bits = %d \n', bit);
% fprintf('  Bias = %d \n', bias);
% fprintf(' Results: \n');
% fprintf('  Target = \n');
% disp(target);
% fprintf('  Actual output = \n');
% disp(actualOut);
% fprintf('  Calculated Hamming distance = %d \n', Hdist);
% fprintf('  Learned Hamming distance = %f \n', actual_Hdist);
% fprintf('  Error = \n');
% disp(error);
% fprintf(' Comment: \n');
% if (error <= 0.1)
%     fprintf('  The readout layer SUCCESSFULLY learned the Hamming distance using the concentration data from the reservoir. \n');
% else
%     fprintf('  The readout layer UNSUCCESSFULLY learned the Hamming distance using the concentration data from the reservoir. \n');
% end
    