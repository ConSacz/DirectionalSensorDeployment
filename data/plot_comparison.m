%y.abc(1,:) is 20 node-deployment of ABC
%y.abc(2,:) is 25 node-deployment of ABC
%y.abc(3,:) is 30 node-deployment of ABC
%y.foa is node-deployment of FOA
close all;
clear;
clc;

%% Initialization 
% x  = 1:10;
% for nodes=20:5:30
%     % ABC result
%     %y.abc = zeros(1, 10);
%     for i = 1:numel(x)
%        name = ['./comparison/barrier/ABC/ABC',num2str(nodes),'nodes_', num2str(i), '.mat']; 
%        load(name, 'BestCost', 'it');
%        y.abc((nodes-20)/5+1,i) = BestCost(it);
%     end
%     clear i name BestCost it;
% 
%     % FOA result
%     %y.foa = zeros(1, 10);
%     for i = 1:numel(x)
%        name = ['./comparison/barrier/FOA/FOA',num2str(nodes),'nodes_', num2str(i), '.mat']; 
%        load(name, 'BestCost', 'it');
%        y.foa((nodes-20)/5+1,i) = BestCost(it);
%     end
%     clear i name BestCost it;
% 
%     % PSO result
%     %y.pso = zeros(1, 10);
%     for i = 1:numel(x)
%        name = ['./comparison/barrier/PSO/PSO',num2str(nodes),'nodes_', num2str(i), '.mat']; 
%        load(name, 'BestCost', 'it');
%        y.pso((nodes-20)/5+1,i) = BestCost(it);
%     end
%     clear i name BestCost it;
% 
%     % GA result
%     %y.ga = zeros(1, 10);
%     for i = 1:numel(x)
%        name = ['./comparison/barrier/GA/GA',num2str(nodes),'nodes_', num2str(i), '.mat']; 
%        load(name, 'BestCost', 'it');
%        y.ga((nodes-20)/5+1,i) = BestCost(it);
%     end
%     clear i name BestCost it;
% end
%% Results illustration
load("./zResult matrixes/barrier_comparison.mat");

%% line chart
% X-axis: Number of sensor nodes
node_counts = [20, 25, 30];

% Compute mean coverage for each algorithm across samples
mean_abc = mean(y.abc, 2)*100;  % average across rows
mean_foa = mean(y.foa, 2)*100;
mean_pso = mean(y.pso, 2)*100;
mean_ga  = mean(y.ga, 2)*100;

% Plot
figure;
plot(node_counts, mean_abc, '-o', 'DisplayName', 'ABC', 'LineWidth', 2); hold on;
plot(node_counts, mean_foa, '-s', 'DisplayName', 'FOA', 'LineWidth', 2);
plot(node_counts, mean_pso, '-^', 'DisplayName', 'PSO', 'LineWidth', 2);
plot(node_counts, mean_ga, '-d', 'DisplayName', 'GA', 'LineWidth', 2);

% Labels and formatting
xlabel('Number of Sensor Nodes');
ylabel('Coverage Rate (%)');
%title('Average Coverage vs Node Count for Different Algorithms');
legend('Location', 'best');
grid on;

%%
% Boxplot
figure;

colors = ['y', 'r', 'c', 'b'];
boxplot([(y.abc(3,:))' (y.pso(3,:))' (y.ga(3,:))' (y.foa(3,:))']*100);
h = findobj(gca,'Tag','Box');
for j=1:length(h)
    patch(get(h(j),'XData'),get(h(j),'YData'),colors(j),'FaceAlpha',.5);
end
%name = [num2str(x(i)), ' nodes deployment'];
%title(name);
xticklabels({'ABC','PSO', 'GA', 'FOA'});
ylabel('Coverage rate (%)')