function [conS, S, freq, R_total, iniCon_ss, iniCon_ds, orderTable, Influx, ...
    timeVec, tmax, tPerturb, numPerturb, thold, perturbTimes, index1, index2, kI_base, nS] = random_DNA_Strand_circuit_ver5(Sm_base, thold)
    % Created by: Hoang Nguyen
    % Program description:
        % Pretty much like random_DNA_Strand_circuit_ver3.m except for
        % considering the volume, initial concentration 0 --> 10^6 (molecules),
        % and update 1 (molecule) every Gillespie iteration

        % Optimal random DNA chemistry using the results from the genetic algorithm
            % (nL, nU, k) = (3, 4, 3)
        % Globally ordered type
        % Fixed partial double strand participation, fixed number of single strands
        % Full double: U1L1
        % Partial double U2L1, U1L3, U1L4, U2L1, U2L3, U2L4, U3L1, U3L2, U3L3

    %clear; clc; close all;    

    %tic

    % random seed    
    RandStream.setGlobalStream(RandStream('mt19937ar','seed',sum(100*clock)));

    % Volume of the chamber
    V = 1e-6;

    % Input parameters:
    % General parameters:
    n = 7 ; % Number of single strands [5, 10)
    p = 3/4; % ratio of upper to lower strands [0.5, 1)
    y = 1/3; % ratio of upper strands with complements [0, 1)
    a_in = 0.99; % ratio of influces to the overall number of strands [0, 1)
    a_out = 0.10; % ratio of outfluces to the overall number of strands [0, 1)
    % Random distribution parameters:
    theta.mean = 0.075; % [0.05, 0.2) % CHANGE?
    theta.variance = 0.01; % [0, 0.02) % CHANGE?
    theta_in.mean = Sm_base; % [0, 0.0006) % CHANGE?
    theta_in.variance = 0; 
    theta_out.mean = 0.0003; % [0, 0.0006) % CHANGE?
    theta_out.variance = 0;
    phi.mean = 3; % [0, 4) % k = 3 --> |phi-3| < 0.5
    phi.variance = 0.001; % [0, 0.5)

    % determine the number of lower and upper strands
    nL = round(n / (1 + p)); % number of lower strands
    nU = n - nL; % number of upper strands

    % Create upper and lower strands as nodes of network (2 partitions)
    U = {'U1', 'U2', 'U3', 'U4', 'U5', 'U6', 'U7', 'U8', 'U9', 'U10'};
    L = {'L1', 'L2', 'L3', 'L4', 'L5', 'L6', 'L7', 'L8', 'L9', 'L10'};
    U = U(1 : nU);
    L = L(1 : nL);

    % Determine the number of full double strands
    nF = y * min(nL, nU); % number of full double strands (complementary single strands)

    % Select nF complementary upper and lower DNA strands randomly
    % Connect them with solid lines
    index_full = 1; % index of the upper and lower strand to form the full double strand (keep fix at 1 for U1L1)
    F = '';
    for i = 1 : nF
        F(i, :) = [U{index_full(i)}, L{index_full(i)}];
    end

    % Draw the number of partial double strands (from distribution phi)
    phi.norm = abs(normrnd(phi.mean, phi.variance, [1, nU])); % positive normal distribution of partial double strands per upper strand
    k = round(phi.norm); %  number of partial double strands per upper strand
    k = k(1); % make k = 3 for all upper strands   

    % Impose ordering of partial double strands for strand displacemen reaction 
     ... (only DNA strands with higher order displace strands from complex)
    nP = k * nU; % total number of partial double strands
    P = ['U1L2'; 'U1L3'; 'U1L4'; 'U2L1'; 'U2L3'; 'U2L4'; 'U3L1'; 'U3L2'; 'U3L3']; P = cellstr(P); P = P';
    orderP = randperm(nP, nP); % generate an order vector randomly
    P = P(orderP); % sort P with order from low to high
    orderP = num2cell(sort(orderP)); % sort the order from low to high
    % Now, the P element i match the order element i
    orderTable = [ orderP, {nP+1} ; P, F ]; % the higher number, the higher order

    % Create a variable containing the set of all species
    S = [ U, L, P, F ];
    nS = numel(S); % Total number of species in the network

    % Choose random strands to influx and outflux 
    nI = round(a_in  * nS); % number of  species to influx
    nO = round(a_out * nS); % number of species to outflux
    I  = S(randperm(nS, nI));  % Choose numIn species to influx
    O  = S(randperm(nS, nO));  % Choose numOut species to outflux 

    % Generate vector of reaction equations
    ss = [U L]; % single strands (including upper and lower strands)
    ds = [P F]; % double strands (full and partial double strands)
    R1 = zeros(length(ds), length(ss));
    R1 = num2cell(R1);

    for i = 1 : length(ds)
        for s = 1 : length(ss)
            if (ss{s}(1) == 'U')
                for j = 1 : length(orderTable)
                    if ([ss{s}, ds{i}(3:4)] == orderTable{2, j})
                        for l = 1 : length(orderTable)
                            if (ds{i} == orderTable{2, l})
                                if (j > l)
                                    R1{i,s} = [ss{s}, ' + ', ds{i}, ' --> ', ss{s}, ds{i}(3:4), ' + ', ds{i}(1:2)];
                                end
                            end
                        end
                    end
                end
            elseif (ss{s}(1) == 'L')
                for j = 1 : length(orderTable)
                    if ([ds{i}(1:2), ss{s}] == orderTable{2, j})
                        for l = 1 : length(orderTable)
                            if (ds{i} == orderTable{2, l})
                                if (j > l)
                                    R1{i,s} = [ss{s}, ' + ', ds{i}, ' --> ', ds{i}(1:2), ss{s}, ' + ', ds{i}(3:4)];
                                end
                            end
                        end
                    end
                end
            end 
        end
    end

    R_displace = zeros(20, 1);
    R_displace = num2cell(R_displace);
    count = 1;
    for i = 1 : numel(R1)
        if (R1{i} ~= 0)
            R_displace{count} = R1{i};
            count = count + 1;
        end
    end

    R2 = zeros(length(U), length(L));
    R2 = num2cell(R2);
    for i = 1 : length(U)
        for j = 1 : length(L)
            for l = 1 : length(orderTable)
                if ([U{i}, L{j}] == orderTable{2, l})
                    R2{i,j} = [U{i}, ' + ', L{j}, ' --> ', U{i}, L{j}];
                end
            end
        end
    end

    R_bind = zeros(length(ds), 1);
    R_bind = num2cell(R_bind);
    count2 = 1;
    for i = 1 : numel(R2)
        if (R2{i} ~= 0)
            R_bind{count2} = R2{i};
            count2 = count2 + 1;
        end
    end

    R_in = zeros(length(I), 1);
    R_in = num2cell(R_in);
    for i = 1 : length(I)
        R_in{i} = ['0 --> ', I{i}];
    end

    R_out = zeros(length(O), 1);
    R_out = num2cell(R_out);
    for j = 1 : length(O)
        R_out{j} = [O{j}, ' --> 0'];
    end

    R = [R_bind; R_displace];
    R_total = [R_in ; R_out ; R_bind ; R_displace];

    % Generate random rate constants (from distribution theta)
    nR = numel(R); % number of reactions (not include influx and outflux reactions)
    k_bind = abs(normrnd(theta.mean, sqrt(theta.variance), [1, length(R_bind)])); % rate constants of reactions except for influx and outflux
    k_displace = abs(normrnd(theta.mean, sqrt(theta.variance), [1, length(R_displace)])); % rate constants of reactions except for influx and outflux

    % Generate random rate constants for influces and outfluces
    ... (from distribution theta_in, theta_out)
    kI_base = abs(normrnd(theta_in.mean , sqrt(theta_in.variance))); % rate constants of influces
    kI = kI_base * ones(nI, 1);
    kO = abs(normrnd(theta_out.mean, sqrt(theta_out.variance), [nO, 1])); % rate constants of outfluces

%     % Print out network generation results
%     fprintf('Network generation parameters:\n');
%     fprintf('\t Total number of species/strands in network: nS = %d \n', nS);
%     fprintf('\t Number of single strands: n = %d \n', n);
%     fprintf('\t Number of upper strands: nU = %d \n', nU);
%     fprintf('\t Number of lower strands: nL = %d \n', nL);
%     fprintf('\t Number of full double strands: nF = %d \n', nF);
%     fprintf('\t Number of partial double strands per upper strand: k = %d \n', k);
%     fprintf('\t Number of partial double strands: nP = %d \n', nP);
%     fprintf('\t Number of influces: nI = %d \n', nI);
%     fprintf('\t Influx rate  = %f \n', kI_base);
%     fprintf('\t Number of outfluces: nO = %d \n', nO);
%     fprintf('\t Outflux rate = %f \n', kO(1));
%     fprintf('\t Number of binding and displacement reactions: nR = %d \n', nR);
%     fprintf('\t Number of all reactions = %d \n', length(R_total));
% 
%     fprintf('\n S = all strands \n U = upper strands \n L = lower strands \n ');
%     fprintf('F = full double strands \n P = partial double strands \n ');
%     fprintf('I = influx strands \n O = outflux strands \n ');
%     fprintf('R = binding and displacement reactions \n ');
%     fprintf('k_bind = binding reaction rate constants \n ');
%     fprintf('k_displace = displacement reaction rate constants \n ');
%     fprintf('R_total = all reactions (influx, outflux, bind, displace) \n ');
%     fprintf('Initial conditions: iniCon_ss and iniCon_ds \n ');
%     fprintf('Show order of double strands: orderTable \n\n ');

    % Gillespie algorithm

    % Generate random initial concentration of the species
    iniCon_ss = randi([0, 1000], [1, n]); % M
    iniCon_ss = num2cell(iniCon_ss);
    iniCon_ds = randi([0, 1000], [1, nF + nP]); % M
    iniCon_ds = num2cell(iniCon_ds);

    con_ss = [ss; iniCon_ss]; % concentration table for single strands
    con_ds = [ds; iniCon_ds]; % concentration table for double strands

    dCon = 1;
    dInflux = 1;
    dOutflux = 1;

    t = 0;
    tmax  = 1; % maximum simulation time
    %thold = 0.2; % (s) % hold the input (influx rate) for thold seconds

    iter = 0;
    %timeVec = zeros(1000,1);
    timeVec(1) = t;

    perturbCount = 0;
    tPerturb = 0.01;
    numPerturb =  ceil((tmax - tPerturb) / thold);
    perturbTimes = zeros(numPerturb, 1);
    Influx = zeros(numPerturb, nS);
    index = randperm(nS, 2);
    index1 = index(1);
    index2 = index(2);

    r_bind = zeros(length(R_bind), 1);
    r_displace = zeros(length(R_displace), 1);
    rO = zeros(length(R_out), 1);
    prob = zeros(length(R_total), 1);

    freq = zeros(length(R_total), 1); % frequency of the reaction (how many time that the reactions happen)

%     fprintf('Perturbed species 1: %s, index %d \n ', S{index1}, index1);
%     fprintf('Perturbed species 2: %s, index %d \n ', S{index2}, index2);
%     fprintf('Simulation time interval: %d to %d \n ', t, tmax);
%     fprintf('Change in concentration for binding and displacement reactions: %f \n ', dCon);
%     fprintf('Change in concentration for influx reactions: %f \n ', dInflux);
%     fprintf('Change in concentration for outflux reactions: %f \n\n ', dOutflux);

    while (t < tmax) && (con_ss{end,1} >= 0) && (con_ss{end,2} >= 0) && (con_ss{end,3} >= 0) && ...
            (con_ss{end,4} >= 0) && (con_ss{end,5} >= 0) && (con_ss{end,6} >= 0) && (con_ss{end,7} >= 0) && ...
            (con_ds{end,1} >= 0) && (con_ds{end,2} >= 0) && (con_ds{end,3} >= 0) && (con_ds{end,4} >= 0) && ...
            (con_ds{end,5} >= 0) && (con_ds{end,6} >= 0) && (con_ds{end,7} >= 0) && (con_ds{end,8} >= 0) && ...
            (con_ds{end,9} >= 0) && (con_ds{end,10} >= 0)
        iter = iter + 1; % number of Iterations

        % perturbation (change in influx rate)
        if (t >= tPerturb)
             perturbTimes(perturbCount+1) = tPerturb;
             tPerturb = tPerturb + thold;
             perturbCount = perturbCount + 1;
             kI(index1) = kI_base * rand;
             kI(index2) = kI_base * rand;
             Influx(perturbCount, :) = kI;
        end

        % Reaction rate for each influx reactions
        rI = kI / V;

        % Reaction rate for each binding reactions
        for i = 1 : length(R_bind)
            for j1 = 1 : size(con_ss, 2)
                if (R_bind{i}(1:2) == con_ss{1,j1})
                    for j2 = 1 : size(con_ss, 2)
                        if (R_bind{i}(6:7) == con_ss{1,j2})
                            r_bind(i) = k_bind(i) * con_ss{iter + 1,j1} * con_ss{iter + 1, j2};
                            break;
                        end
                    end
                    break;
                end
            end
        end

        % Reaction rate for each displacement reactions
        for i = 1 : length(R_displace)
            for j1 = 1 : size(con_ss, 2)
                if (R_displace{i}(1:2) == con_ss{1,j1})
                    for j2 = 1 : size(con_ds, 2)
                        if (R_displace{i}(6:9) == con_ds{1,j2})
                            r_displace(i) = k_displace(i) * con_ss{iter + 1,j1} * con_ds{iter + 1, j2};
                            break;
                        end
                    end
                    break;
                end
            end
        end

        % Reaction rate for outflux reactions
        for i = 1 : length(R_out)
            space1st = find(R_out{i} == ' ');
            space1st = space1st(1);
            for j2 = 1 : size(con_ds, 2)
                for j1 = 1 : size(con_ss, 2)
                    if (space1st == 5)
                        if (R_out{i}(1:4) == con_ds{1,j2})
                            rO(i) = kO(i) * con_ds{iter + 1, j2} / V;
                            break;
                        end
                    elseif (space1st == 3)
                        if (R_out{i}(1:2) == con_ss{1,j1})
                            rO(i) = kO(i) * con_ss{iter + 1, j1} / V;
                            break;
                        end
                    end
                end
            end
        end

        % Total reaction rate
        rVec = [rI; rO; r_bind; r_displace]; % vector of all reaction rates
        rTot = sum(rVec); % sum of all reaction rates

        % Probability for each reaction to happen
        for i = 1 : length(rVec)
           prob(i) = rVec(i) / rTot; % probability
        end  

        % increase time step by dt
        mu = 1 / rTot; % mean of the exponential distribution
        dt = exprnd(mu);
        t = t + dt;
        timeVec(iter + 1) = t;

         % update species concentration
         x = sum(rand >= cumsum([0, prob']));

         freq(x) = freq(x) + 1; % frequency of the reaction (how many time that the reactions happen)

         if (x >= 1) && (x <= 17)
            space1st = find(R_total{x} == ' ');
            spacediff = length(R_total{x}) - space1st(1);
            for j2 = 1 : size(con_ds, 2)
                con_ds{(iter + 2), j2} = con_ds{(iter + 1), j2};
                for j1 = 1 : size(con_ss, 2)
                    con_ss{(iter + 2), j1} = con_ss{(iter + 1), j1};
                    if (spacediff == 8)
                        if (R_total{x}(7:10) == con_ds{1, j2})
                            con_ds{(iter + 2), j2} = con_ds{(iter + 1), j2} + dInflux;
                        end
                    elseif (spacediff == 6)
                        if (R_total{x}(7:8) == con_ss{1, j1})
                            con_ss{(iter + 2), j1} = con_ss{(iter + 1), j1} + dInflux;
                        end
                    end
                end
            end
         elseif (x >= 18) && (x <= 19)
             space1st = find(R_total{x} == ' ');
             space1st = space1st(1);
             for j2 = 1 : size(con_ds, 2)
                 con_ds{(iter + 2), j2} = con_ds{(iter + 1), j2};
                 for j1 = 1 : size(con_ss, 2)
                     con_ss{(iter + 2), j1} = con_ss{(iter + 1), j1};
                     if (space1st == 5)
                         if (R_total{x}(1:4) == con_ds{1, j2})
                             %if (con_ds{(iter + 1), j2} - dOutflux) >= 0 
                                 con_ds{(iter + 2), j2} = con_ds{(iter + 1), j2} - dOutflux;
                             %end
                         end
                     elseif (space1st == 3)
                         if (R_total{x}(1:2) == con_ss{1, j1})
                             %if (con_ss{(iter + 1), j1} - dOutflux) >= 0 
                                 con_ss{(iter + 2), j1} = con_ss{(iter + 1), j1} - dOutflux;
                             %end
                         end
                     end
                 end
             end
         elseif (x >= 20) && (x <= 29)
             for j1 = 1 : size(con_ss, 2)
                if (R_total{x}(1:2) == con_ss{1, j1})
                    for j2 = 1 : size(con_ss, 2)
                        if (R_total{x}(6:7) == con_ss{1, j2})
                            for j3 = 1 : size(con_ds, 2)
                                if (R_total{x}(13:16) == con_ds{1, j3})
                                    for i = 1 : size(con_ss, 2)
                                        con_ss{(iter + 2), i} = con_ss{(iter + 1), i};
                                    end
                                    for ii = 1 : size(con_ds, 2)
                                        con_ds{(iter + 2), ii} = con_ds{(iter + 1), ii};
                                    end
                                    %if ((con_ss{(iter + 1), j1} - dCon) >= 0)
                                        con_ss{(iter + 2), j1} = con_ss{(iter + 1), j1} - dCon;
                                    %end
                                    %if ((con_ss{(iter + 1), j2} - dCon) >= 0)
                                        con_ss{(iter + 2), j2} = con_ss{(iter + 1), j2} - dCon;
                                    %end
                                    con_ds{(iter + 2), j3} = con_ds{(iter + 1), j3} + dCon;
                                end
                            end
                        end
                    end
                end
             end
         elseif (x >= 30) && (x <= 50)
             for j1 = 1 : size(con_ss, 2)
                 if (R_total{x}(1:2) == con_ss{1, j1})
                    for j2 = 1 : size(con_ds, 2)
                        if (R_total{x}(6:9) == con_ds{1, j2})
                            for j3 = 1 : size(con_ds, 2)
                                if (R_total{x}(15:18) == con_ds{1, j3})
                                    for j4 = 1 : size(con_ss, 2)
                                        if (R_total{x}(22:23) == con_ss{1, j4})
                                            for i = 1 : size(con_ss, 2)
                                                con_ss{(iter + 2), i} = con_ss{(iter + 1), i};
                                            end
                                            for ii = 1 : size(con_ds, 2)
                                                con_ds{(iter + 2), ii} = con_ds{(iter + 1), ii};
                                            end
                                            %if (con_ss{(iter + 1), j1} - dCon) >= 0
                                                con_ss{(iter + 2), j1} = con_ss{(iter + 1), j1} - dCon;
                                            %end
                                            %if (con_ds{(iter + 1), j2} - dCon) >= 0
                                                con_ds{(iter + 2), j2} = con_ds{(iter + 1), j2} - dCon;
                                            %end
                                            con_ds{(iter + 2), j3} = con_ds{(iter + 1), j3} + dCon;
                                            con_ss{(iter + 2), j4} = con_ss{(iter + 1), j4} + dCon;
                                        end
                                    end
                                end
                            end
                        end
                    end
                 end
             end
         end

%          % Make the concentration of a species 0 if the concentration of that
%          % species end up < 0 after a Gillespie iteration
%     %      for i1 = 1 : size(con_ss, 2)
%     %          if (con_ss{end,i1} < 0)
%     %             con_ss{end,i1} = 0;
%     %          end
%     %      end
%     %      
%     %      for i2 = 1 : size(con_ds, 2)
%     %          if (con_ds{end,i2} < 0)
%     %              con_ds{end,i2} = 0;
%     %          end
%     %      end
    end

    % Create concentration vs. time plot
    timeVec = timeVec';
    conS = [con_ss, con_ds];

    figure;
    hold on;
    grid;
    %title('Random DNA Strand Circuit');
    xlabel('Time (s)');
    ylabel('Number of molecules');
    plot(timeVec, cell2mat(conS(2:end, :)));
    for i = 1 : size(conS, 2)
        legendInfo{i} = num2str(conS{1, i});
        legend(legendInfo);
    end
    hold off;

    %toc
end