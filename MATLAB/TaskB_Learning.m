Plot_Influx_TaskA_TaskB;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Task B Computation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

t1st_TaskB = time_vector - thold;
t2nd_TaskB = time_vector - (3/2)*thold;

t_B = time_vector(find(abs(t2nd_TaskB) <= 1e-10) : end);

% Task B for input 1
Y_B1 = zeros(length(time_vector), 1);
for i = find(abs(t2nd_TaskB) <= 1e-10) : length(time_vector)
    Y_B1(i) = Influx_vector1(abs(time_vector - t1st_TaskB(i)) <= 1e-10) + ...
        (1/2)*Influx_vector1(abs(time_vector - t2nd_TaskB(i)) <= 1e-10);
end

% % Task B for input 2
% Y_B2 = zeros(length(time_vector), 1);
% for i = find(abs(t2nd_TaskB) <= 1e-10) : length(time_vector)
%     Y_B2(i) = Influx_vector2(abs(time_vector - t1st_TaskB(i)) <= 1e-10) + ...
%         (1/2)*Influx_vector2(abs(time_vector - t2nd_TaskB(i)) <= 1e-10);
% end

% % Plot Task B output
% figure;
% hold on;
% grid on;
% plot(t_B, Y_B1(find(abs(t2nd_TaskB) <= 1e-10) : end), ...
%     t_B, Y_B2(find(abs(t2nd_TaskB) <= 1e-10) : end), 'LineWidth', 2);
% titleName = sprintf('Task B Computation');
% title(titleName);
% xlabel('Time');
% ylabel('Y_{B}(t)');
% legendName1 = sprintf('Input 1: Influx of %s', S{index1});
% legendName2 = sprintf('Input 2: Influx of %s', S{index2});
% legend(legendName1, legendName2);
% hold off;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Task B Learning %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

bias = 1; % Bias
coeff = 0.01; % Learning coefficient
iterations = 1000; % Number of iterations

target_B = Y_B1(Y_B1 ~= 0); % get the non-zero elements in output of Task A (elements that was filled in)
target_B_norm = target_B / (max(target_B) + min(target_B)); % normalize target

% Number of bits
bit = 10; 
num20bits = floor(length(target_B_norm) / bit); % number of sections that have #bit bits
lastbit = length(target_B_norm) - num20bits * bit; % number of bits in the last section

target = zeros(bit, num20bits);
error_20 = zeros(bit, num20bits);
actualOut_20 = zeros(bit, num20bits);

for i20 = 1 : num20bits
    target(:,i20) = target_B_norm(1 + bit*(i20-1) : bit*i20); % sections of #bit targets 
    
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

target_last = target_B_norm(length(target_B_norm)-lastbit+1 : length(target_B_norm)); % last section of target

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
actualOut_B_norm = [reshape(actualOut_20,numel(actualOut_20),1) ; actualOut_last];
actualOut_B = actualOut_B_norm * (max(target_B) + min(target_B));

% fprintf('\nResevoir computing: \n');
% fprintf(' Setup: \n');
% fprintf('  Number of iteration = %d \n', iterations);
% fprintf('  Learning rate = %g \n',  coeff);
% fprintf('  Number of bits = %d \n', bit);
% fprintf('  Bias = %d \n', bias);
% fprintf(' Results: \n');
% fprintf('  Target = \n');
% disp(target_B);
% fprintf('  Actual output = \n');
% disp(actualOut_B);
% fprintf('  Error = \n');
% disp(error);

% Plot Task B output
figure;
hold on;
grid on;
plot(t_B, target_B, t_B, actualOut_B, 'LineWidth', 2);
%titleName = sprintf('Task B');
%title(titleName);
xlabel('Time');
ylabel('Y(t)');
legendName1 = sprintf('Target');
legendName2 = sprintf('Actual Output');
legend(legendName1, legendName2);
hold off;

% Normalized Root Mean Square Error (NRMSE) of Task B
NRMSE_B = (1 / (max(actualOut_B) - min(actualOut_B))) * sqrt(sum((target_B - actualOut_B).^2) / length(t_B));
% fprintf('Task B NRMSE = %f\n', NRMSE_B);

