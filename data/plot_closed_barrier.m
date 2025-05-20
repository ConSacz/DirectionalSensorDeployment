% y.case20 is result of case r1=20 from 20 to 45 nodes

close all;
clear;
clc;

%% Initialization 
% x  = 1:1;
% for nodes=20:5:45
%     % case20 result
%     for i = 1:numel(x)
%        name = ['./closed_barrier_coverage/barrier 20/',num2str(nodes),' nodes_', num2str(i), '.mat']; 
%        load(name, 'BestCost', 'it');
%        y.case_20((nodes-20)/5+1) = BestCost(it);
%     end
%     clear i name BestCost it;
% 
%     % case30 result
%     for i = 1:numel(x)
%        name = ['./closed_barrier_coverage/barrier 30/',num2str(nodes),' nodes_', num2str(i), '.mat'];
%        load(name, 'BestCost', 'it');
%        y.case_30((nodes-20)/5+1) = BestCost(it);
%     end
%     clear i name BestCost it;
% 
%     % case40 result
%     for i = 1:numel(x)
%        name = ['./closed_barrier_coverage/barrier 40/',num2str(nodes),' nodes_', num2str(i), '.mat']; 
%        load(name, 'BestCost', 'it');
%        y.case_40((nodes-20)/5+1) = BestCost(it);
%     end
%     clear i name BestCost it;
% 
%     % case50 result
%     for i = 1:numel(x)
%        name = ['./closed_barrier_coverage/barrier 50/',num2str(nodes),' nodes_', num2str(i), '.mat']; 
%        load(name, 'BestCost', 'it');
%        y.case_50((nodes-20)/5+1) = BestCost(it);
%     end
%     clear i name BestCost it;
% 
%     % case20-30 result
%     for i = 1:numel(x)
%        name = ['./closed_barrier_coverage/barrier 20-30/',num2str(nodes),' nodes_', num2str(i), '.mat']; 
%        load(name, 'BestCost', 'it');
%        y.case_20_30((nodes-20)/5+1) = BestCost(it);
%     end
%     clear i name BestCost it;
% 
%     % case30-40 result
%     for i = 1:numel(x)
%        name = ['./closed_barrier_coverage/barrier 30-40/',num2str(nodes),' nodes_', num2str(i), '.mat']; 
%        load(name, 'BestCost', 'it');
%        y.case_30_40((nodes-20)/5+1) = BestCost(it);
%     end
%     clear i name BestCost it;
% end
%% Results illustration
load('closed_barrier.mat');
%nodes = 20:5:45;
cases=1:numel(fieldnames(y));
bar_data = [
    y.case_20; y.case_30; y.case_40; y.case_50; y.case_20_30; y.case_30_40
]*100;

% Create a grouped bar chart
figure;
bar(cases, bar_data, 'grouped');

% Label the x-axis
%xlabel('Number of Sensor Nodes');
xticklabels({'Case 1','Case 2', 'Case 3', 'Case 4', 'Case 5', 'Case 6'});

% Label the y-axis
ylabel('Coverage rate (%)');

% Add legend for the different cases
legend({'20-node', '25-node', '30-node', '35-node', '40-node', '45-node'}, 'Location', 'best');
%%
% Boxplot
% figure;
% 
% colors = ['y', 'r', 'c', 'b'];
% boxplot([(y.case_3(1,:))' (y.case_3(2,:))' (y.case_3(3,:))' ]*100);
% h = findobj(gca,'Tag','Box');
% for j=1:length(h)
%     patch(get(h(j),'XData'),get(h(j),'YData'),colors(j),'FaceAlpha',.5);
% end
% %name = [num2str(x(i)), ' nodes deployment'];
% %title(name);
% xticklabels({'20-node deployment','25-node deployment', '30-node deployment'});
% ylabel('Coverage rate (%)')