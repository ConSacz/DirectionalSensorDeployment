%% Directional Sensor
%% DEPLOYMENT

%%
clc;
clear;

%% 
close all;
%% Network parameter
for trial=1:10
name=['./data/comparison/barrier/ITLBO/TLBO30nodes_',num2str(trial),'.mat' ];
% Monitor area
Covered_Area = zeros(100,100);
Obstacle_Area = gen_random_distribution_area(30,29);
%Obstacle_Area = ones(100,100);

% nodes info
MaxIt = 200;              % Maximum Number of Iterations
N = 30;
nPop=50;
rc = 20;

%homogenerous sensor
rs=ones(1,N)*10;
theta0=ones(1,N)*pi/3;
sink=[50 50];

%% Init first pop
empty_individual.Position=[]; % Empty solution
empty_individual.Cost=[];     % Empty cost function of that solution
BestSol.Position=[];	% Best Solution
BestSol.Cost=-1;	% Best solution function value
pop=repmat(empty_individual,nPop,1); % pop includes nPop solutions
BestCost=zeros(MaxIt,1); % Best solution every iteration

clear empty_individual;

for k = 1:nPop
    alpop = unifrnd(sink(1)-rc,sink(2)+rc,[N 2]);
    % gen sink node
    %alpop (1,:)= sink;
    alpop (:,3)= unifrnd(0,2*pi,[N 1]);

    pop(k).Position=alpop;
    [pop(k).Cost,~]=Cov_Func(pop(k).Position,rs,theta0,Obstacle_Area,Covered_Area);

    % check for best Solution
    if pop(k).Cost>BestSol.Cost
	    BestSol.Cost=pop(k).Cost;
	    BestSol.Position=pop(k).Position;
    end 
end
clear k alpop sink;
%%  TLBO Main Loop
for it = 1:MaxIt
    Partner = randperm(nPop);
    for i=1:nPop
        %% Teaching phase
        mean_pop = zeros(N,3);
        for k = 1:nPop
            mean_pop = mean_pop+pop(k).Position;
        end
        mean_pop = mean_pop/nPop;
        
        TF = randi([1 2],1,1);
        alpop = pop(i).Position + unifrnd(0, +1, [N 3]).*(BestSol.Position - TF*mean_pop);
        alpop(:,1:2) = min(max(alpop(:,1:2), 1),size(Obstacle_Area,1)-1);
        alpop(:,3)=mod(alpop(:,3),2*pi);

        if Connectivity_graph(Graph(alpop(:,1:2),rc),[]) == 1
            [alpop_cov,~]  = Cov_Func(alpop,rs,theta0,Obstacle_Area,Covered_Area);
            if alpop_cov >= pop(i).Cost
                pop(i).Position=alpop;
                pop(i).Cost=alpop_cov;
                if alpop_cov > BestSol.Cost
	                BestSol.Cost=alpop_cov;
	                BestSol.Position=alpop;
                end 
            end
        end

        %% Learning phase
        for j = 1:N
            alpop = pop(i).Position;
            if (pop(i).Cost > pop(Partner(i)).Cost)
                alpop(j,:) = pop(i).Position(j,:) + unifrnd(0, +1, [1 3]).*(pop(i).Position(j,:) - pop(Partner(i)).Position(j,:));
            else
                alpop(j,:) = pop(i).Position(j,:) + unifrnd(0, +1, [1 3]).*(pop(Partner(i)).Position(j,:) - pop(i).Position(j,:));
            end
            alpop(:,1:2) = min(max(alpop(:,1:2), 1),size(Obstacle_Area,1)-1);
            alpop(:,3)=mod(alpop(:,3),2*pi);
    
            if Connectivity_graph(Graph(alpop(:,1:2),rc),[]) == 1
                [alpop_cov,~]  = Cov_Func(alpop,rs,theta0,Obstacle_Area,Covered_Area);
                if alpop_cov >= pop(i).Cost
                    pop(i).Position=alpop;
                    pop(i).Cost=alpop_cov;
                    if alpop_cov > BestSol.Cost
	                    BestSol.Cost=alpop_cov;
	                    BestSol.Position=alpop;
                    end 
                end
            end  
        end
    end

    clear alpop_cov alpop TF mean_pop Partner;
    BestCost(it)=BestSol.Cost;
    disp(['Iteration ' num2str(it) ': Best Cost = ' num2str(BestCost(it))]);
    if BestCost(it)==1
        break;
    end    
end
%% plot
% %alpop=pop(4).Position;
% alpop=BestSol.Position;
% [coverage,Covered_Area] = Cov_Func(alpop,rs,theta0,Obstacle_Area,Covered_Area);
% % show map
% %imagesc(imresize(Obstacle_Area,[100 100],'bilinear'))
% figure;
% plot(BestCost);
% figure;
% hold on;
% 
% % show map
% [obs_row, obs_col] = find(Obstacle_Area == 1);
% plot(obs_row, obs_col,'.', 'MarkerSize', 10, 'Color', 'blue');
% [obs_row, obs_col] = find(Covered_Area == 1/2);
% plot(obs_row, obs_col,'.', 'MarkerSize', 2, 'Color', 'green');
% [obs_row, obs_col] = find(Covered_Area == 1);
% plot(obs_row, obs_col,'.', 'MarkerSize', 10, 'Color', 'cyan');
% %[discovered_obs_row, discovered_obs_col] = find(Covered_Area == -1);                    % show discovered map
% %plot(discovered_obs_row-1, discovered_obs_col-1,'.', 'MarkerSize', 20, 'Color', 'red');
% %colorbar;
% %{
% for i = 1:1:numel(G.Edges.EndNodes)/2
%     plot([pop(G.Edges.EndNodes(i,1)*2-1),pop(G.Edges.EndNodes(i,2)*2-1)],[pop(G.Edges.EndNodes(i,1)*2),pop(G.Edges.EndNodes(i,2)*2)],'Color','blue','linewidth',1);
% end
% %}
% for i = 1:N
%     % i-th node
%     plot (alpop(i,2) , alpop(i,1),'ro','MarkerSize', 3,'Color','red');
%     text (alpop(i,2) , alpop(i,1), num2str(i),'FontSize',10,'Color','red');
%     theta = linspace(alpop(i,3)-theta0(i)/2 , alpop(i,3) + theta0(i)/2 , 100); 
%     x = alpop(i,1) + rs(i) * cos(theta);
%     y = alpop(i,2) + rs(i) * sin(theta);
%     x_fill = [alpop(i,1), x, alpop(i,1)]; 
%     y_fill = [alpop(i,2), y, alpop(i,2)];
%     fill(y_fill, x_fill, 'g', 'FaceAlpha', 0.3, 'EdgeColor', 'g');
% 
% end
% 
% %clear i x_fill y_fill theta obs_col obs_row;
% 
% axis equal;
% xlim([0 size(Obstacle_Area,1)]);
% ylim([0 size(Obstacle_Area,2)]);
% title(['Weighted Coverage ratio:' num2str(BestSol.Cost*100) '% ' ]);
% grid on;
% drawnow;
%{
if save_fig==1
    saveas(gcf, ['Iteration' num2str(it) '.pdf' ]); % Lưu dưới dạng PDF
end
%}
%%
save(name)
end

