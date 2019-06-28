clear; close all; clc;

load('/u/hoang24/MATLAB/Random_Chemistry/NRMSE_Experiments_Data/29-Mar-2019 16:10:13.mat');

% Data to be plotted as a bar graph
model_series = [0.23 , 0.11 ; 0.17 , 0.13 ; min(min(avg_NRMSE_A)) , min(min(avg_NRMSE_B))];

%Data to be plotted as the error bars
model_error = [0.05 , 0.02 ; 0.034 , 0.036 ; ...
    std_NRMSE_A(find(min(min(avg_NRMSE_A)))) , std_NRMSE_B(find(min(min(avg_NRMSE_B))))];

figure;
% Creating axes and the bar graph
ax = axes;
h = bar(model_series,'BarWidth',1);
% Set color for each bar face
% h(1).FaceColor = 'blue';
% h(2).FaceColor = 'yellow';
% Properties of the bar graph as required
ax.YGrid = 'on';
ax.GridLineStyle = '-';
xticks(ax,[1 2 3]);

% Naming each of the bar groups
xticklabels(ax,{'Goudarzi, 2013', 'Yahiro, 2018', 'Our model'});

% title, X and Y labels
xlabel ('CRN Model');
ylabel ('Error (NRMSE)');

% Creating a legend and placing it outside the bar plot
lg = legend('Short-term Memory Task', 'Long-term Memory Task', 'AutoUpdate', 'off');
lg.Location = 'BestOutside';
lg.Orientation = 'Horizontal';
hold on;
% Finding the number of groups and the number of bars in each group
ngroups = size(model_series, 1);
nbars = size(model_series, 2);
% Calculating the width for each bar group
groupwidth = min(0.8, nbars/(nbars + 1.5));
% Set the position of each error bar in the centre of the main bar
% Based on barweb.m by Bolu Ajiboye from MATLAB File Exchange
for i = 1:nbars
    % Calculate center of each bar
    x = (1:ngroups) - groupwidth/2 + (2*i-1) * groupwidth / (2*nbars);
    errorbar(x, model_series(:,i), model_error(:,i), 'k', 'linestyle', 'none');
end

better_than_Goudarzi_A = 0.23 - min(min(avg_NRMSE_A));
better_than_Goudarzi_B = 0.11 - min(min(avg_NRMSE_B));

better_than_Yahiro_A = 0.17 - min(min(avg_NRMSE_A));
better_than_Yahiro_B = 0.13 - min(min(avg_NRMSE_B));

% Plot the number of species

% figure;
% hold on;
% xlabel('CRN Model');
% ylabel('Number of Species Types');
% modelname = categorical({'Goudarzi, 2013', 'Yahiro, 2018', 'This model'});
% numspecies = [3 , 3 , 17];
% bar(modelname , numspecies);
% hold off;

numspecies = [3 ; 3 ; 17];
figure;
% Creating axes and the bar graph
ax = axes;
h = bar(numspecies,'BarWidth',0.5);
% Set color for each bar face
% h(1).FaceColor = 'blue';
% h(2).FaceColor = 'yellow';
% Properties of the bar graph as required
ax.YGrid = 'on';
ax.GridLineStyle = '-';
xticks(ax,[1 2 3]);

% Naming each of the bar groups
xticklabels(ax,{'Goudarzi, 2013', 'Yahiro, 2018', 'Our model'});

% title, X and Y labels
xlabel ('CRN Model');
ylabel ('Number of Species Types');

% grid on;

% Creating a legend and placing it outside the bar plot
% lg = legend('Short-term Memory Task', 'Long-term Memory Task', 'AutoUpdate', 'off');
% lg.Location = 'BestOutside';
% lg.Orientation = 'Horizontal';
%hold on;
% Finding the number of groups and the number of bars in each group
%ngroups = size(model_series, 1);
%nbars = size(model_series, 2);
% Calculating the width for each bar group
%groupwidth = min(0.8, nbars/(nbars + 1.5));
% Set the position of each error bar in the centre of the main bar
% Based on barweb.m by Bolu Ajiboye from MATLAB File Exchange
%for i = 1:nbars
    % Calculate center of each bar
    %x = (1:ngroups) - groupwidth/2 + (2*i-1) * groupwidth / (2*nbars);
    %errorbar(x, model_series(:,i), model_error(:,i), 'k', 'linestyle', 'none');
%end
