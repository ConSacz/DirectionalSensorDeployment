% y.u(1,:) is result of 20 nodes on case U barrier 
% y.u(2,:) is result of 25 nodes on case U barrier 
% y.u(3,:) is result of 30 nodes on case U barrier 
% y.v is on case V barrier

close all;
clear;
clc;

%% Initialization 
% x  = 1:10;
% for nodes=20:5:30
%     % U barrier result
%     for i = 1:numel(x)
%        name = ['./open_barrier_coverage/U_barrier/',num2str(nodes),' nodes/U_', num2str(i), '.mat']; 
%        load(name, 'BestCost', 'it');
%        y.u((nodes-20)/5+1,i) = BestCost(it);
%     end
%     clear i name BestCost it;
% 
%     % case2 result
%     for i = 1:numel(x)
%        name = ['./open_barrier_coverage/V_barrier/',num2str(nodes),' nodes/V_', num2str(i), '.mat']; 
%        load(name, 'BestCost', 'it');
%        y.v((nodes-20)/5+1,i) = BestCost(it);
%     end
%     clear i name BestCost it;
% 
% end
%% Results illustration

load("./zResult matrixes/open_barrier.mat");
% Covergence 
figure;
hold on;
xlabel('Time step');
ylabel('Coverage rate (%)');
axis([0 800 0 100]);
colors = lines(5);

%plot(x1, y.map1it, 'color', "#0000FF");
%plot(x1, y.map1_2it, 'color', "#00FFFF");
%plot(x1, y.map2it, 'color', "#D95319");
%plot(x1, y.map3it, 'color', "#EDB120");
%plot(x, mean(y.abc), ':*', 'color', "#EDB120");
for j = 5:10:40
    plot(x1, y.map3it(j,:), 'color', colors(round(j/10), :));
end
plot(x1, y.map3it(48,:), 'color', colors(round(5), :));

%legend(arrayfun(@(i) ['Line ' num2str(i)], 1:5, 'UniformOutput', false),'Location', 'best');
legend('first run','second run','third run', 'fourth run', 'fifth run' ,'Location', 'best');
grid on;
%%
% Boxplot
figure;

colors = ['y', 'r', 'c', 'b'];
boxplot([(y.u(1,:))' (y.u(2,:))' (y.u(3,:))' ]*100);
h = findobj(gca,'Tag','Box');
for j=1:length(h)
    patch(get(h(j),'XData'),get(h(j),'YData'),colors(j),'FaceAlpha',.5);
end
%name = [num2str(x(i)), ' nodes deployment'];
%title(name);
xticklabels({'20-node deployment','25-node deployment', '30-node deployment'});
ylabel('Coverage rate (%)')