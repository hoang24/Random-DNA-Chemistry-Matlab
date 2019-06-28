Plot_Influx_TaskA_TaskB;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Task A Computation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

t1st_TaskA = time_vector - 0.001;
t2nd_TaskA = time_vector - 0.002;

t_A = time_vector(find(t2nd_TaskA == 0) : end);

% Task A for input 1
Y_A1 = zeros(length(time_vector), 1);
for i = find(t2nd_TaskA == 0) : length(time_vector)
    Y_A1(i) = Influx_vector1(abs(time_vector - t1st_TaskA(i)) <= 1e-10) + ...
        2*Influx_vector1(abs(time_vector - t2nd_TaskA(i)) <= 1e-10);
end

% % Task A for input 2
% Y_A2 = zeros(length(time_vector), 1);
% for i = find(t2nd_TaskA == 0) : length(time_vector)
%     Y_A2(i) = Influx_vector2(abs(time_vector - t1st_TaskA(i)) <= 1e-10) + ...
%         2*Influx_vector2(abs(time_vector - t2nd_TaskA(i)) <= 1e-10);
% end

% % Plot Task A output
% figure;
% hold on;
% grid on;
% plot(t_A, Y_A1(find(t2nd_TaskA == 0) : end), ...
%     t_A, Y_A2(find(t2nd_TaskA == 0) : end), 'LineWidth', 2);
% titleName = sprintf('Task A Computation');
% title(titleName);
% xlabel('Time');
% ylabel('Y_{A}(t)');
% legendName1 = sprintf('Input 1: Influx of %s', S{index1});
% legendName2 = sprintf('Input 2: Influx of %s', S{index2});
% legend(legendName1, legendName2);
% hold off;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Task A Learning %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

bias = 1; % Bias
coeff = 0.01; % Learning coefficient
iterations = 1000; % Number of iterations

target_A = Y_A1(Y_A1 ~= 0); % get the non-zero elements in output of Task A (elements that was filled in)
target_A_norm = target_A / (max(target_A) + min(target_A)); % normalize target

% Number of bits
bit = 10; 
num20bits = floor(length(target_A_norm) / bit); % number of sections that have #bit bits
lastbit = length(target_A_norm) - num20bits * bit; % number of bits in the last section

target = zeros(bit, num20bits);
error_20 = zeros(bit, num20bits);
actualOut_20 = zeros(bit, num20bits);

for i20 = 1 : num20bits
    target(:,i20) = target_A_norm(1 + bit*(i20-1) : bit*i20); % sections of #bit targets 
    
    % Concentration dataset
    conS_2 = cell2mat(conS(2:end,:)); % Convert concentration dataset from cell type to matrix type
    conS_3 = zeros(size(conS_2, 1), nS); % Pre-allocation for the full concentration dataset
    conS_4 = zeros(bit, nS); % Pre-allocation for the scaled concentration dataset (to #bits)
    for iCon = 1 : nS
        conS_3(:,iCon) = conS_2(:,iCon) / (max(conS_2(:,1)) + min(conS_2(:,1))); % Normalized dataset
        conS_4(:,iCon) = conS_3(randi(size(conS_2, 1), [bit, 1])); % pick a random region of the dataset that has the same length as the number of bits
    end

    % Weights
    wBias = rand;
    wSpecies = rand(1,nS);

    % Preallocation
    actualOut_matrix = zeros(bit, iterations);
    error_matrix = zeros(bit, iterations);
    weights = [wSpecies, wBias];

    for i = 1 : iterations
        for j = 1 : bit
            net = sum(conS_4(j,:) .* wSpecies) + bias * wBias;    
            actualOut_matrix(j,i) = 1./(1+exp(-0.05*net));
            error_matrix(j,i) = target(j, i20) - actualOut_matrix(j,i);

            wSpecies = wSpecies + coeff * conS_4(j,:) .* error_matrix(j,i); % species weight update
            wBias = wBias + coeff * bias * error_matrix(j,i); % bias weight update

            weights = [weights; wSpecies, wBias];
        end
    end
    error_matrix = abs(error_matrix);
    error_20(:, i20) = error_matrix(:,end);
    actualOut_20(:, i20) = actualOut_matrix(:,end);    
end

target_last = target_A_norm(length(target_A_norm)-lastbit+1 : length(target_A_norm)); % last section of target

% Concentration dataset
conS_2 = cell2mat(conS(2:end,:)); % Convert concentration dataset from cell type to matrix type
conS_3 = zeros(size(conS_2, 1), nS); % Pre-allocation for the full concentration dataset
conS_4 = zeros(lastbit, nS); % Pre-allocation for the scaled concentration dataset (to #bits)
for iCon = 1 : nS
    conS_3(:,iCon) = conS_2(:,iCon) / (max(conS_2(:,1)) + min(conS_2(:,1))); % Normalized dataset
    conS_4(:,iCon) = conS_3(randi(size(conS_2, 1), [lastbit, 1])); % pick a random region of the dataset that has the same length as the number of bits
end

% Weights
wBias = rand;
wSpecies = rand(1,nS);

% Preallocation
actualOut_matrix = zeros(lastbit, iterations);
error_matrix = zeros(lastbit, iterations);
weights = [wSpecies, wBias];

for i = 1 : iterations
    for j = 1 : lastbit
        net = sum(conS_4(j,:) .* wSpecies) + bias * wBias;    
        actualOut_matrix(j,i) = 1./(1+exp(-0.05*net));
        error_matrix(j,i) = target_last(j) - actualOut_matrix(j,i);

        wSpecies = wSpecies + coeff * conS_4(j,:) .* error_matrix(j,i); % species weight update
        wBias = wBias + coeff * bias * error_matrix(j,i); % bias weight update

        weights = [weights; wSpecies, wBias];
    end
end
error_matrix = abs(error_matrix);
error_last = error_matrix(:,end);
actualOut_last = actualOut_matrix(:,end); 

error = [reshape(error_20,numel(error_20),1) ; error_last];
actualOut_A_norm = [reshape(actualOut_20,numel(actualOut_20),1) ; actualOut_last];
actualOut_A = actualOut_A_norm * (max(target_A) + min(target_A));

% fprintf('\nResevoir computing: \n');
% fprintf(' Setup: \n');
% fprintf('  Number of iteration = %d \n', iterations);
% fprintf('  Learning rate = %g \n',  coeff);
% fprintf('  Number of bits = %d \n', bit);
% fprintf('  Bias = %d \n', bias);
% fprintf(' Results: \n');
% fprintf('  Target = \n');
% disp(target_A);
% fprintf('  Actual output = \n');
% disp(actualOut_A);
% fprintf('  Error = \n');
% disp(error);

% Plot Task A output
figure;
hold on;
grid on;
plot(t_A, target_A, t_A, actualOut_A, 'LineWidth', 2);
%titleName = sprintf('Task A');
%title(titleName);
xlabel('Time');
ylabel('Y(t)');
legendName1 = sprintf('Target');
legendName2 = sprintf('Actual Output');
legend(legendName1, legendName2);
hold off;

% Normalized Root Mean Square Error (NRMSE) of Task A
NRMSE_A = (1 / (max(actualOut_A) - min(actualOut_A))) * sqrt(sum((target_A - actualOut_A).^2) / length(t_A));
% fprintf('Task A NRMSE = %f\n', NRMSE_A);



% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Task A Learning %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% % target_A = Y_A1(Y_A1 ~= 0);
% % target_A_norm = round(target_A / (max(target_A) + min(target_A)));
% % 
% % % Number of bits
% % bit = length(target_A_norm);
% % 
% % % Concentration dataset
% % conS_2 = cell2mat(conS(2:end,:)); % Convert concentration dataset from cell type to matrix type
% % conS_3 = zeros(size(conS_2, 1), nS); % Pre-allocation for the full concentration dataset
% % conS_4 = zeros(bit, nS); % Pre-allocation for the scaled concentration dataset (to #bits)
% % for i = 1 : nS
% %     conS_3(:,i) = conS_2(:,i) / (max(conS_2(:,1)) + min(conS_2(:,1))); % Normalized dataset
% %     conS_4(:,i) = conS_3(randi(size(conS_2, 1), [bit, 1])); % pick a random region of the dataset that has the same length as the number of bits
% % end
% % 
% % % Bias
% % bias = 1;
% % 
% % % Weights
% % wBias = rand;
% % wSpecies = rand(1,nS);
% % 
% % % Learning coefficient
% % coeff = 0.7;
% % 
% % % Number of iterations
% % iterations = 1000;
% % 
% % % Reallocation
% % actualOut_matrix = zeros(bit, iterations);
% % error_matrix = zeros(bit, iterations);
% % weights = [wSpecies, wBias];
% % 
% % for i = 1 : iterations
% %     for j = 1 : bit
% %         net = sum(conS_4(j,:) .* wSpecies) + bias * wBias;    
% %         actualOut_matrix(j,i) = 1./(1+exp(-net));
% %         error_matrix(j,i) = target_A_norm(j) - actualOut_matrix(j,i);
% %         
% %         wSpecies = wSpecies + coeff * conS_4(j,:) .* error_matrix(j,i); % species weight update
% %         wBias = wBias + coeff * bias * error_matrix(j,i); % bias weight update
% % 
% %         weights = [weights; wSpecies, wBias];
% %     end
% % end
% % error_matrix = abs(error_matrix);
% % error = error_matrix(:,end);
% % actualOut = actualOut_matrix(:,end);
% % 
% % fprintf('\nResevoir computing: \n');
% % fprintf(' Setup: \n');
% % fprintf('  Number of iteration = %d \n', iterations);
% % fprintf('  Learning rate = %g \n',  coeff);
% % fprintf('  Number of bits = %d \n', bit);
% % fprintf('  Bias = %d \n', bias);
% % fprintf(' Results: \n');
% % fprintf('  Target = \n');
% % disp(target_A_norm);
% % fprintf('  Actual output = \n');
% % disp(actualOut);
% % fprintf('  Error = \n');
% % disp(error);
% % fprintf(' Comment: \n');
% % if (error <= 0.1)
% %     fprintf('  The readout layer SUCCESSFULLY learned the Hamming distance using the concentration data from the reservoir. \n');
% % else
% %     fprintf('  The readout layer UNSUCCESSFULLY learned the Hamming distance using the concentration data from the reservoir. \n');
% % end
