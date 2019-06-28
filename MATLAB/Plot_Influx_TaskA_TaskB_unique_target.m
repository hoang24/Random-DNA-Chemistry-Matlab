%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Run chemistry %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

random_DNA_Strand_circuit_ver5;

%%%%%%%%%%%%%%%%%%%%%%%%%% Create vector of influx rates as a function of time %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% original inputs
original_input1 = Influx(:, index1);
original_input2 = Influx(:, index2);
time_vector = (0:0.001:1)'; % vector of time

% Input 1
Influx_vector1 = zeros(length(time_vector), 1);
for i0 = 1 : find(abs(time_vector - perturbTimes(1)) <= 1e-10)-1
   Influx_vector1(i0) = kI_base; 
end
for i1 = find(abs(time_vector - perturbTimes(1)) <= 1e-10) : find(abs(time_vector - perturbTimes(2)) <= 1e-10)-1 
   Influx_vector1(i1) = original_input1(1);
end
for i2 = find(abs(time_vector - perturbTimes(2)) <= 1e-10) : find(abs(time_vector - perturbTimes(3)) <= 1e-10)-1
   Influx_vector1(i2) = original_input1(2);
end
for i3 = find(abs(time_vector - perturbTimes(3)) <= 1e-10) : find(abs(time_vector - perturbTimes(4)) <= 1e-10)-1 
   Influx_vector1(i3) = original_input1(3);
end
for i4 = find(abs(time_vector - perturbTimes(4)) <= 1e-10) : find(abs(time_vector - perturbTimes(5)) <= 1e-10)-1 
   Influx_vector1(i4) = original_input1(4);
end
for i5 = find(abs(time_vector - perturbTimes(5)) <= 1e-10) : find(time_vector == tmax)
   Influx_vector1(i5) = original_input1(5);
end

% Input 2
Influx_vector2 = zeros(length(time_vector), 1);
for i0 = 1 : find(abs(time_vector - perturbTimes(1)) <= 1e-10)-1
   Influx_vector2(i0) = kI_base; 
end
for i1 = find(abs(time_vector - perturbTimes(1)) <= 1e-10) : find(abs(time_vector - perturbTimes(2)) <= 1e-10)-1 
   Influx_vector2(i1) = original_input2(1);
end
for i2 = find(abs(time_vector - perturbTimes(2)) <= 1e-10) : find(abs(time_vector - perturbTimes(3)) <= 1e-10)-1
   Influx_vector2(i2) = original_input2(2);
end
for i3 = find(abs(time_vector - perturbTimes(3)) <= 1e-10) : find(abs(time_vector - perturbTimes(4)) <= 1e-10)-1 
   Influx_vector2(i3) = original_input2(3);
end
for i4 = find(abs(time_vector - perturbTimes(4)) <= 1e-10) : find(abs(time_vector - perturbTimes(5)) <= 1e-10)-1 
   Influx_vector2(i4) = original_input2(4);
end
for i5 = find(abs(time_vector - perturbTimes(5)) <= 1e-10) : find(time_vector == tmax)
   Influx_vector2(i5) = original_input2(5);
end

% Table of time and influx rate as a function of time
Influx_matrix = [time_vector, Influx_vector1, Influx_vector2];

% Plot influx rate as a function of time
figure;
hold on;
grid on;
plot(time_vector, Influx_vector1, time_vector, Influx_vector2, 'LineWidth', 2);
titleName = sprintf('Influx rates (reservoir inputs) over simulation time');
title(titleName);
xlabel('Time');
ylabel('Influx rates');
legendName1 = sprintf('Input 1: Influx of %s', S{index1});
legendName2 = sprintf('Input 2: Influx of %s', S{index2});
legend(legendName1, legendName2);
hold off;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Task A Computation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

t1st_TaskA = time_vector - 0.001;
t2nd_TaskA = time_vector - 0.002;

% Task A for input 1
Y_A1 = zeros(length(time_vector), 1);
for i = find(t2nd_TaskA == 0) : length(time_vector)
    Y_A1(i) = Influx_vector1(abs(time_vector - t1st_TaskA(i)) <= 1e-10) + ...
        2*Influx_vector1(abs(time_vector - t2nd_TaskA(i)) <= 1e-10);
end

% Task A for input 2
Y_A2 = zeros(length(time_vector), 1);
for i = find(t2nd_TaskA == 0) : length(time_vector)
    Y_A2(i) = Influx_vector2(abs(time_vector - t1st_TaskA(i)) <= 1e-10) + ...
        2*Influx_vector2(abs(time_vector - t2nd_TaskA(i)) <= 1e-10);
end

% Plot Task A output
figure;
hold on;
grid on;
plot(time_vector(find(t2nd_TaskA == 0) : end), Y_A1(find(t2nd_TaskA == 0) : end), ...
    time_vector(find(t2nd_TaskA == 0) : end), Y_A2(find(t2nd_TaskA == 0) : end), 'LineWidth', 2);
titleName = sprintf('Task A Computation');
title(titleName);
xlabel('Time');
ylabel('Y_{A}(t)');
legendName1 = sprintf('Input 1: Influx of %s', S{index1});
legendName2 = sprintf('Input 2: Influx of %s', S{index2});
legend(legendName1, legendName2);
hold off;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Task B Computation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

t1st_TaskB = time_vector - thold;
t2nd_TaskB = time_vector - (3/2)*thold;

% Task B for input 1
Y_B1 = zeros(length(time_vector), 1);
for i = find(abs(t2nd_TaskB) <= 1e-10) : length(time_vector)
    Y_B1(i) = Influx_vector1(abs(time_vector - t1st_TaskB(i)) <= 1e-10) + ...
        (1/2)*Influx_vector1(abs(time_vector - t2nd_TaskB(i)) <= 1e-10);
end

% Task B for input 2
Y_B2 = zeros(length(time_vector), 1);
for i = find(abs(t2nd_TaskB) <= 1e-10) : length(time_vector)
    Y_B2(i) = Influx_vector2(abs(time_vector - t1st_TaskB(i)) <= 1e-10) + ...
        (1/2)*Influx_vector2(abs(time_vector - t2nd_TaskB(i)) <= 1e-10);
end

% Plot Task B output
figure;
hold on;
grid on;
plot(time_vector(find(abs(t2nd_TaskB) <= 1e-10) : end), Y_B1(find(abs(t2nd_TaskB) <= 1e-10) : end), ...
    time_vector(find(abs(t2nd_TaskB) <= 1e-10) : end), Y_B2(find(abs(t2nd_TaskB) <= 1e-10) : end), 'LineWidth', 2);
titleName = sprintf('Task B Computation');
title(titleName);
xlabel('Time');
ylabel('Y_{B}(t)');
legendName1 = sprintf('Input 1: Influx of %s', S{index1});
legendName2 = sprintf('Input 2: Influx of %s', S{index2});
legend(legendName1, legendName2);
hold off;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Task A Learning %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

target_A = Y_A1(Y_A1 ~= 0); % target vector 
target_A_unique = unique(target_A, 'stable'); % target vector with only unique values
target_A_norm = target_A_unique / (max(target_A) + min(target_A)); % normalize unique target vector

% Number of bits
bit = length(target_A_norm);

% Concentration dataset
conS_2 = cell2mat(conS(2:end,:)); % Convert concentration dataset from cell type to matrix type
conS_3 = zeros(size(conS_2, 1), nS); % Pre-allocation for the full concentration dataset
conS_4 = zeros(bit, nS); % Pre-allocation for the scaled concentration dataset (to #bits)
for i = 1 : nS
    conS_3(:,i) = conS_2(:,i) / (max(conS_2(:,1)) + min(conS_2(:,1))); % Normalized dataset
    conS_4(:,i) = conS_3(randi(size(conS_2, 1), [bit, 1])); % pick a random region of the dataset that has the same length as the number of bits
end

% Bias
bias = 1;

% Weights
wBias = rand;
wSpecies = rand(1,nS);

% Learning coefficient
coeff = 0.01;

% Number of iterations
iterations = 1000;

% Preallocation
actualOut_matrix = zeros(bit, iterations);
error_matrix = zeros(bit, iterations);
weights = [wSpecies, wBias];

for i = 1 : iterations
    for j = 1 : bit
        net = sum(conS_4(j,:) .* wSpecies) + bias * wBias;    
        actualOut_matrix(j,i) = 1./(1+exp(-0.5 * net));
        error_matrix(j,i) = target_A_norm(j) - actualOut_matrix(j,i);
        
        wSpecies = wSpecies + coeff * conS_4(j,:) .* error_matrix(j,i); % species weight update
        wBias = wBias + coeff * bias * error_matrix(j,i); % bias weight update

        weights = [weights; wSpecies, wBias];
    end
end
error_matrix = abs(error_matrix);
error = error_matrix(:,end);
actualOut_A_norm = actualOut_matrix(:,end);
actualOut_A_unique = actualOut_A_norm * (max(target_A) + min(target_A));

actualOut_A(1:9) = actualOut_A_unique(1);
actualOut_A(10) = actualOut_A_unique(2);
actualOut_A(11:209) = actualOut_A_unique(3);
actualOut_A(210) = actualOut_A_unique(4);
actualOut_A(211:409) = actualOut_A_unique(5);
actualOut_A(410) = actualOut_A_unique(6);
actualOut_A(411:609) = actualOut_A_unique(7);
actualOut_A(610) = actualOut_A_unique(8);
actualOut_A(611:809) = actualOut_A_unique(9);
actualOut_A(810) = actualOut_A_unique(10);
actualOut_A(811:999) = actualOut_A_unique(11);

% for i0 = 1 : 9
%     actualOut_A(
% end

fprintf('\nResevoir computing: \n');
fprintf(' Setup: \n');
fprintf('  Number of iteration = %d \n', iterations);
fprintf('  Learning rate = %g \n',  coeff);
fprintf('  Number of bits = %d \n', bit);
fprintf('  Bias = %d \n', bias);
fprintf(' Results: \n');
fprintf('  Target = \n');
disp(target_A);
fprintf('  Actual output = \n');
disp(actualOut_A);
fprintf('  Error = \n');
disp(error);
fprintf(' Comment: \n');
% if (error <= 0.1)
%     fprintf('  The readout layer SUCCESSFULLY learned the Hamming distance using the concentration data from the reservoir. \n');
% else
%     fprintf('  The readout layer UNSUCCESSFULLY learned the Hamming distance using the concentration data from the reservoir. \n');
% end

% Plot Task A output
figure;
hold on;
grid on;
plot(time_vector(find(t2nd_TaskA == 0) : end), target_A, ...
    time_vector(find(t2nd_TaskA == 0) : end), actualOut_A, 'LineWidth', 2);
titleName = sprintf('Task A');
title(titleName);
xlabel('Time');
ylabel('Y_{A}(t)');
legendName1 = sprintf('Computation');
legendName2 = sprintf('Learning');
legend(legendName1, legendName2);
hold off;

