clc; clear; close all;

% thold = [0.2 0.3 0.4 0.5 1];
% tPerturb = 0.01;
% tmax = 1;
% numPerturb = ceil((tmax - tPerturb) ./ thold);

% mean_velocity = [0.2574, 0.1225, 0.1787]; % mean velocity
% std_velocity = [0.3314, 0.2278, 0.2836];  % standard deviation of velocity
% figure
% hold on
% bar(1:3,mean_velocity)
% errorbar(1:3,mean_velocity,std_velocity,'.')

% a=[0,1,0,0;
% 4,3,2,1;
% 2,2,1,3;
% 1,0,0,0];
% b=[0,1,0,0;
% 1,2,1,1;
% 1,1,1,2;
% 1,0,0,0];
% ctrs = 1:4;
% data = a;
% figure(1)
% hBar = bar(ctrs, data);
% for k1 = 1:size(a,1)
%     ctr(k1,:) = bsxfun(@plus, hBar(1).XData, [hBar(k1).XOffset]');
%     ydt(k1,:) = hBar(k1).YData;
% end
% hold on
% errorbar(ctr, ydt, b, '.r')
% hold off

% x = [1998 1999 2000 2001];
% y = [75 91 105 ; 123.5 131 150 ; 179 203 226 ; 249 281.5 269];
% errhigh = [2.1 4.4 0.4 ; 3.3 2.5 0.4 ; 1.6 0.8 0.6 ; 0.8 2.2 0.9];
% errlow  = [4.4 2.4 2.3 ; 0.5 1.6 1.5 ; 4.5 1.5 0.4 ; 1.2 1.3 0.8];
% bar(x, y);
% legend('A','B','C');
% hold on
% er = errorbar(x,y,errlow,errhigh);    
% er.Color = [0 0 0];                            
% er.LineStyle = 'none';  
% hold off

% x = [1 2];
% data = [37.6 24.5 ; 14.6 18.1]';
% errhigh = [2.1 4.4 ; 0.4 3.3];
% errlow  = [4.4 2.4 ; 2.3 0.5];
% 
% bar(x,data)                
% 
% hold on
% 
% er = errorbar(x,data,errlow,errhigh);    
% er.Color = [0 0 0];                            
% er.LineStyle = 'none';  
% 
% hold off

% % Data to be plotted as a bar graph
% model_series = [10 40 50 60; 20 50 60 70; 30 60 80 90; 100 110 120 130];
% 
% %Data to be plotted as the error bars
% model_error = [1 4 8 6; 2 5 9 12; 3 6 10 13; 14 15 16 17];
% 
% % Creating axes and the bar graph
% ax = axes;
% h = bar(model_series,'BarWidth',1);
% % Set color for each bar face
% % h(1).FaceColor = 'blue';
% % h(2).FaceColor = 'yellow';
% % Properties of the bar graph as required
% ax.YGrid = 'on';
% ax.GridLineStyle = '-';
% xticks(ax,[1 2 3 4]);
% 
% % Naming each of the bar groups
% xticklabels(ax,{'0.0001', '0.0002', '0.0003', '0.0004'});
% 
% % X and Y labels
% xlabel ('Input range (Base influx rate)');
% ylabel ('Error (NRMSE)');
% 
% % Creating a legend and placing it outside the bar plot
% lg = legend('\tau = 0.2','\tau = 0.3','\tau = 0.4','\tau = 0.5','AutoUpdate','off');
% lg.Location = 'BestOutside';
% lg.Orientation = 'Horizontal';
% hold on;
% % Finding the number of groups and the number of bars in each group
% ngroups = size(model_series, 1);
% nbars = size(model_series, 2);
% % Calculating the width for each bar group
% groupwidth = min(0.8, nbars/(nbars + 1.5));
% % Set the position of each error bar in the centre of the main bar
% % Based on barweb.m by Bolu Ajiboye from MATLAB File Exchange
% for i = 1:nbars
%     % Calculate center of each bar
%     x = (1:ngroups) - groupwidth/2 + (2*i-1) * groupwidth / (2*nbars);
%     errorbar(x, model_series(:,i), model_error(:,i), 'k', 'linestyle', 'none');
% end

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
