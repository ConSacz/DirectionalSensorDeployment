% y.case1(1,:) is result of 20 nodes on case 1 
% y.case1(2,:) is result of 25 nodes on case 1 
% y.case1(3,:) is result of 30 nodes on case 1 
% y.case2 is case 2
% y.case3 is case 3
% y.case4 is case 4
close all;
clear;
clc;

%% Initialization 
% x  = 1:20;
% for nodes=20:5:30
%     % case1 result
%     for i = 1:numel(x)
%        name = ['./hetero_barrier/',num2str(nodes),' nodes/case1/hetero_', num2str(i), '.mat']; 
%        load(name, 'BestCost', 'it');
%        y.case_1((nodes-20)/5+1,i) = BestCost(it);
%     end
%     clear i name BestCost it;
% 
%     % case2 result
%     for i = 1:numel(x)
%        name = ['./hetero_barrier/',num2str(nodes),' nodes/case2/hetero_', num2str(i), '.mat']; 
%        load(name, 'BestCost', 'it');
%        y.case_2((nodes-20)/5+1,i) = BestCost(it);
%     end
%     clear i name BestCost it;
% 
%     % case3 result
%     for i = 1:numel(x)
%        name = ['./hetero_barrier/',num2str(nodes),' nodes/case3/hetero_', num2str(i), '.mat']; 
%        load(name, 'BestCost', 'it');
%        y.case_3((nodes-20)/5+1,i) = BestCost(it);
%     end
%     clear i name BestCost it;
% 
%     % case4 result
%     for i = 1:numel(x)
%        name = ['./hetero_barrier/',num2str(nodes),' nodes/case4/hetero_', num2str(i), '.mat']; 
%        load(name, 'BestCost', 'it');
%        y.case_4((nodes-20)/5+1,i) = BestCost(it);
%     end
%     clear i name BestCost it;
% end
%% Results illustration
load("./zResult matrixes/hetero_barrier.mat");

%%
% Boxplot
figure;

colors = ['y', 'r', 'c', 'b'];
boxplot([(y.case_3(1,:))' (y.case_3(2,:))' (y.case_3(3,:))' ]*100);
h = findobj(gca,'Tag','Box');
for j=1:length(h)
    patch(get(h(j),'XData'),get(h(j),'YData'),colors(j),'FaceAlpha',.5);
end
%name = [num2str(x(i)), ' nodes deployment'];
%title(name);
xticklabels({'20-node deployment','25-node deployment', '30-node deployment'});
ylabel('Coverage rate (%)')