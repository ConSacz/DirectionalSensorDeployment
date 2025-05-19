% y.case1 is case 1
% y.case2 is case 2
% y.case3 is case 3
% y.case4 is case 4
close all;
clear;
clc;

%% Initialization 
x  = 1:5;
% case1 result
for i = 1:numel(x)
   name = ['./target_coverage/case 1/case1_', num2str(i), '.mat']; 
   load(name, 'BestCost', 'it');
   y.case_1(1,i) = BestCost(it);
end
clear i name BestCost it;

% case2 result
for i = 1:numel(x)
   name = ['./target_coverage/case 2/case2_', num2str(i), '.mat']; 
   load(name, 'BestCost', 'it');
   y.case_2(1,i) = BestCost(it);
end
clear i name BestCost it;

% case3 result
for i = 1:numel(x)
   name = ['./target_coverage/case 3/case3_', num2str(i), '.mat']; 
   load(name, 'BestCost', 'it');
   y.case_3(1,i) = BestCost(it);
end
clear i name BestCost it;

% case4 result
for i = 1:numel(x)
   name = ['./target_coverage/case 4/case4_', num2str(i), '.mat']; 
   load(name, 'BestCost', 'it');
   y.case_4(1,i) = BestCost(it);
end
clear i name BestCost it;

%% Results illustration
load("./zResult matrixes/target.mat");

%%
% Boxplot
figure;

colors = ['y', 'r', 'c', 'b'];
boxplot([(y.case_1(1,:))' (y.case_2(1,:))' (y.case_3(1,:))' (y.case_4(1,:))']*100);
h = findobj(gca,'Tag','Box');
for j=1:length(h)
    patch(get(h(j),'XData'),get(h(j),'YData'),colors(j),'FaceAlpha',.5);
end
%name = [num2str(x(i)), ' nodes deployment'];
%title(name);
xticklabels({'Case 1','Case 2', 'Case 3', 'Case 4'});
ylabel('Coverage rate (%)')

