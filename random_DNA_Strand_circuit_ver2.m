%function random_DNA_Strand_circuit () 
    % Created by: Hoang Nguyen
    % Program description:
        % Optimal random DNA chemistry using the results from the genetic algorithm
            % (nL, nU, k) = (3, 4, 3)
        % Globally ordered type    

    clear; clc; close all;  
    
    tic
    % random seed    
    RandStream.setGlobalStream(RandStream('mt19937ar','seed',sum(100*clock)));
    
    % Input parameters:
        % General parameters:
        n = 7 ; % Number of single strands [5, 10)
        p = 3/4; % ratio of upper to lower strands [0.5, 1)
        y = 1/3; % ratio of upper strands with complements [0, 1)
        a_in = 0.99; % ratio of influces to the overall number of strands [0, 1)
        a_out = 0.10; % ratio of outfluces to the overall number of strands [0, 1)
    
        % Random distribution parameters
        theta.mean = 0.075; % [0.05, 0.2)
        theta.variance = 0.01; % [0, 0.02)
        theta_in.mean = 0.00055; % [0, 0.0006)
        theta_in.variance = 0;
        theta_out.mean = 0.0003; % [0, 0.0006)
        theta_out.variance = 0;
        phi.mean = 3; % [0, 4) % k = 3 --> |phi-3| < 0.5
        phi.variance = 0.01; % [0, 0.5)
       
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
    index_full = randperm(min(nL, nU), nF); % index of the upper and lower strand to form the full double strand
    
    F = '';
    for i = 1 : nF
        F(i, :) = [U{index_full(i)}, L{index_full(i)}];
    end
    
    % Draw the number of partial double strands (from distribution phi)
    phi.norm = abs(normrnd(phi.mean, phi.variance, [1, nU])); % positive normal distribution of partial double strands per upper strand
    k = round(phi.norm); %  number of partial double strands per upper strand
    k = k(1); % make k = 3 for all upper strands   
    
    partial_double = zeros(nU, k + 1);
    partial_double = num2cell(partial_double);
    for j = 1 : nU % for each upper strand
        % Choose randomly number of lower strand counterparts without repetitions
        cp = randperm(nL, k); % counterpart lower strands of each upper strand to form partial double strand
        index_full_vec = ones(1, k) * index_full;
        while ((j == index_full) && sum((cp == index_full_vec) ~= 0))
            cp = randperm(nL, k); % counterpart lower strands of each upper strand to form partial double strand
        end
        
        partial_double(j,1) = U(j);
        for l = 1 : k
            partial_double(j, 1+l) = L(cp(l));
        end
    end
    
    P = zeros(nU, k);
    P = num2cell(P);
    for i1 = 1 : nU
        for i2 = 1 : k
            P{i1,i2} = [partial_double{i1,1}, partial_double{i1,1+i2}];
        end
    end
    
    % Mirror selection of partial double strands for upper strands with complements in the pool
        % No need to since there is only 1 full double strand  
    
    % Impose ordering of partial double strands for strand displacemen reaction 
     ... (only DNA strands with higher order displace strands from complex)
    nP = k * nU; % total number of partial double strands
    P = reshape(rot90(P), [1, nP]);
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
    ss = [U L]; % single strands (including upper and lower strands
    ds = [P F];
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

    RI = zeros(length(I), 1);
    RI = num2cell(RI);
    for i = 1 : length(I)
        RI{i} = ['0 --> ', I{i}];
    end
    
    RO = zeros(length(O), 1);
    RO = num2cell(RO);
    for j = 1 : length(O)
        RO{j} = [O{j}, ' --> 0'];
    end
    
    R = [R_bind; R_displace];
    R_total = [RI ; RO ; R_bind ; R_displace];
    
   	% Generate random rate constants (from distribution theta)
    nR = numel(R); % number of reactions (not include influx and outflux reactions)
    k_bind = abs(normrnd(theta.mean, sqrt(theta.variance), [1, length(R_bind)])); % rate constants of reactions except for influx and outflux
    k_displace = abs(normrnd(theta.mean, sqrt(theta.variance), [1, length(R_displace)])); % rate constants of reactions except for influx and outflux
    
    % Generate random rate constants for influces and outfluces
    ... (from distribution theta_in, theta_out)
    rI = abs(normrnd(theta_in.mean , sqrt(theta_in.variance) , [nI, 1])); % rate constants of influces
    kO = abs(normrnd(theta_out.mean, sqrt(theta_out.variance), [nO, 1])); % rate constants of outfluces
    
    % Print out network generation results
    fprintf("Network generation parameters:\n");
    fprintf("\t Total number of species/strands in network: nS = %d \n", nS);
    fprintf("\t Number of single strands: n = %d \n", n);
    fprintf("\t Number of upper strands: nU = %d \n", nU);
    fprintf("\t Number of lower strands: nL = %d \n", nL);
    fprintf("\t Number of full double strands: nF = %d \n", nF);
    fprintf("\t Number of partial double strands per upper strand: k = %d \n", k);
    fprintf("\t Number of partial double strands: nP = %d \n", nP);
    fprintf("\t Number of influces: nI = %d \n", nI);
    fprintf("\t Influx rate  = %f \n", rI(1));
    fprintf("\t Number of outfluces: nO = %d \n", nO);
    fprintf("\t Outflux rate = %f \n", kO(1));
    fprintf("\t Number of binding and displacement reactions: nR = %d \n", nR);
    
    fprintf("\n S = all strands \n U = upper strands \n L = lower strands \n ");
    fprintf("F = full double strands \n P = partial double strands \n ");
    fprintf("I = influx strands \n O = outflux strands \n ");
    fprintf("R = reactions \n r = reaction rate constants \n ");
    fprintf("Initial conditions: iniCond \n ");
    fprintf("Show order of strands: orderTable \n");
    
    % Gillespie algorithm
        
    % Generate random initial concentration of the species
    iniCon_ss = rand([1, n]); % M
    iniCon_ss = num2cell(iniCon_ss);
    iniCon_ds = rand([1, nF + nP]); % M
    iniCon_ds = num2cell(iniCon_ds);
    
    con_ss = [ss; iniCon_ss]; % concentration table for single strands
    con_ds = [ds; iniCon_ds]; % concentration table for double strands
    
    dCon = 0.01;
    dInflux = 2;
    dOutflux = dInflux;
    
    t = 0;
    tmax = 10000;
    thold = 500; % (s) % hold the input (influx rate) for thold seconds
    numPerturb = (tmax - t) / thold;
    
    iter = 0;
    %timeVec = zeros(1000,1);
    timeVec(1) = t;
    
    perturbCount = 0;
    tPerturb = t;
    
    Influx = zeros(numPerturb, nS);
    
    r_bind = zeros(length(R_bind), 1);
    r_displace = zeros(length(R_displace), 1);
    rO = zeros(length(RO), 1);
    prob = zeros(length(R_total), 1);
    
    freq = zeros(length(R_total), 1); % frequency of the reaction (how many time that the reactions happen)
    
    while (t < tmax)
        iter = iter + 1; % number of Iterations
        
        if (t >= tPerturb)
            tPerturb = tPerturb + thold;
            perturbCount = perturbCount + 1;
            rI = rI;% * rand;
            Influx(perturbCount, :) = rI;
        end
        
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
        for i = 1 : length(RO)
            space1st = find(RO{i} == ' ');
            space1st = space1st(1);
            for j2 = 1 : size(con_ds, 2)
                for j1 = 1 : size(con_ss, 2)
                    if (space1st == 5)
                        if (RO{i}(1:4) == con_ds{1,j2})
                            rO(i) = kO(i) * con_ds{iter + 1, j2};
                            break;
                        end
                    elseif (space1st == 3)
                        if (RO{i}(1:2) == con_ss{1,j1})
                            rO(i) = kO(i) * con_ss{iter + 1, j1};
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
                             if (con_ds{(iter + 1), j2} - dOutflux) >= 0 
                                 con_ds{(iter + 2), j2} = con_ds{(iter + 1), j2} - dOutflux;
                             end
                         end
                     elseif (space1st == 3)
                         if (R_total{x}(1:2) == con_ss{1, j1})
                             if (con_ss{(iter + 1), j1} - dOutflux) >= 0 
                                 con_ss{(iter + 2), j1} = con_ss{(iter + 1), j1} - dOutflux;
                             end
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
                                    if ((con_ss{(iter + 1), j1} - dCon) >= 0)
                                        con_ss{(iter + 2), j1} = con_ss{(iter + 1), j1} - dCon;
                                    end
                                    if ((con_ss{(iter + 1), j2} - dCon) >= 0)
                                        con_ss{(iter + 2), j2} = con_ss{(iter + 1), j2} - dCon;
                                    end
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
                                            if (con_ss{(iter + 1), j1} - dCon) >= 0
                                                con_ss{(iter + 2), j1} = con_ss{(iter + 1), j1} - dCon;
                                            end
                                            if (con_ds{(iter + 1), j2} - dCon) >= 0
                                                con_ds{(iter + 2), j2} = con_ds{(iter + 1), j2} - dCon;
                                            end
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
    end
    timeVec = timeVec';
    conS = [con_ss, con_ds];
    toc