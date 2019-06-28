%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Run chemistry %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[conS, S, freq, R_total, iniCon_ss, iniCon_ds, orderTable, Influx, ...
    timeVec, tmax, tPerturb, numPerturb, thold, perturbTimes, index1, index2, kI_base, nS] = random_DNA_Strand_circuit_ver5(Sm_base, thold);

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
for pT = 1 : length(perturbTimes)-1
    for imid = find(abs(time_vector - perturbTimes(pT)) <= 1e-10) : find(abs(time_vector - perturbTimes(pT+1)) <= 1e-10)-1 
       Influx_vector1(imid) = original_input1(pT);
    end
end
for ilast = find(abs(time_vector - perturbTimes(end)) <= 1e-10) : find(time_vector == tmax)
   Influx_vector1(ilast) = original_input1(end);
end

% % Input 2
% Influx_vector2 = zeros(length(time_vector), 1);
% for i0 = 1 : find(abs(time_vector - perturbTimes(1)) <= 1e-10)-1
%    Influx_vector2(i0) = kI_base; 
% end
% for i1 = find(abs(time_vector - perturbTimes(1)) <= 1e-10) : find(abs(time_vector - perturbTimes(2)) <= 1e-10)-1 
%    Influx_vector2(i1) = original_input2(1);
% end
% for i2 = find(abs(time_vector - perturbTimes(2)) <= 1e-10) : find(abs(time_vector - perturbTimes(3)) <= 1e-10)-1
%    Influx_vector2(i2) = original_input2(2);
% end
% for i3 = find(abs(time_vector - perturbTimes(3)) <= 1e-10) : find(abs(time_vector - perturbTimes(4)) <= 1e-10)-1 
%    Influx_vector2(i3) = original_input2(3);
% end
% for i4 = find(abs(time_vector - perturbTimes(4)) <= 1e-10) : find(abs(time_vector - perturbTimes(5)) <= 1e-10)-1 
%    Influx_vector2(i4) = original_input2(4);
% end
% for i5 = find(abs(time_vector - perturbTimes(5)) <= 1e-10) : find(time_vector == tmax)
%    Influx_vector2(i5) = original_input2(5);
% end

% % Table of time and influx rate as a function of time
% Influx_matrix = [time_vector, Influx_vector1, Influx_vector2];

% % Plot influx rate as a function of time
% figure;
% hold on;
% grid on;
% plot(time_vector, Influx_vector1, time_vector, Influx_vector2, 'LineWidth', 2);
% titleName = sprintf('Influx rates (reservoir inputs) over simulation time');
% title(titleName);
% xlabel('Time');
% ylabel('Influx rates');
% legendName1 = sprintf('Input 1: Influx of %s', S{index1});
% legendName2 = sprintf('Input 2: Influx of %s', S{index2});
% legend(legendName1, legendName2);
% hold off;

