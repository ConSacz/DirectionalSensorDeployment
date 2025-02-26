%% random distribution ROI
%% DEPLOYMENT

%%
clc;
clear;

%% 
%for TIME=1:50
%name = ['./case study robustness/map2 v1.2 robust/map2_', num2str(TIME), '.mat']; 
%load(name);
close all;
%% Network parameter

% Monitor area
Covered_Area = zeros(100,100);
Obstacle_Area = gen_random_distribution_area();
%Obstacle_Area = ones(100,100);
% nodes info
MaxIt = 50;              % Maximum Number of Iterations
a = 1;                    % Acceleration Coefficient Upper Bound
N = 40;
nPop=50;
rc = 10;
rs = 10;
theta0=pi/3;
sink=[50 50];
Scout_bee=nPop;
Onlooker_bee=nPop;

% bee parameter
a=1;
L=zeros(1,nPop);

%% Init first pop

empty_individual.Position=[]; % Empty solution
empty_individual.Cost=[];     % Empty cost function of that solution
BestSol.Position=[];	% Best Solution
BestSol.Cost=0;	% Best solution function value
pop=repmat(empty_individual,nPop,1); % pop includes nPop solutions
BestCost=zeros(MaxIt,1); % Best solution every iteration

clear empty_individual;

for k = 1:nPop
    alpop = unifrnd(sink(1)-rc/2,sink(2)+rc/2,[N 2]);
    % gen sink node
    alpop (1,:)= sink;
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
%%     ABC Main Loop
for it = 1:MaxIt
    %% Global search
    for i=1:Scout_bee
        k=randi([1,nPop]);
        phi=a*unifrnd(-1, +1, [N 3])*(1-L(i)/MaxIt)^2;
        alpop = pop(i).Position + phi.*(pop(i).Position-pop(k).Position);
        alpop(:,1:2) = min(max(alpop(:,1:2), 1),size(Obstacle_Area,1)-1);
        alpop2 = alpop(:,3);
        alpop2(alpop2>2*pi)=alpop2(alpop2>2*pi)-2*pi;
        alpop2(alpop2<0)=alpop2(alpop2<0)+2*pi;
        alpop(:,3)=alpop2;
        if Connectivity_graph(Graph(alpop(:,1:2),rc),[]) == 1
            [pop_cov,~] = Cov_Func(pop(i).Position,rs,theta0,Obstacle_Area,Covered_Area);
            [alpop_cov,~]  = Cov_Func(alpop,rs,theta0,Obstacle_Area,Covered_Area);
            if alpop_cov >= pop_cov
                pop(i).Position=alpop;
                pop(i).Cost=alpop_cov;
            else
                L(i)=L(i)+1;
            end
            break;
        end
    end

    clear alpop_cov alpop alpop2 k phi;
    %% ranking pop
    E=zeros(1,nPop);
    for i=1:nPop
        E(i)=pop(i).Cost;
    end
    E_ave=E/sum(E);

    %% Local search
    for j=1:Onlooker_bee
        %randomly choose a pop, prioritize that have high coverage
        i=randsrc(1,1,[1:nPop;E_ave]);
        
        for k=1:N
            alpop=pop(i).Position;
            h=randi([1,N]);
            phi=a*unifrnd(-1, +1, [1 3]);
            alpop(k,1:3)  = pop(i).Position(k,1:3) + phi.*(pop(i).Position(k,1:3)-pop(i).Position(h,1:3));
            alpop(:,1:2) = min(max(alpop(:,1:2), 1),size(Obstacle_Area,1)-1);
            alpop2 = alpop(:,3);
            alpop2(alpop2>2*pi)=alpop2(alpop2>2*pi)-2*pi;
            alpop2(alpop2<0)=alpop2(alpop2<0)+2*pi;
            alpop(:,3)=alpop2;
            if Connectivity_graph(Graph(alpop(:,1:2),rc),[]) == 1
                [pop_cov,~] = Cov_Func(pop(i).Position,rs,theta0,Obstacle_Area,Covered_Area);
                [alpop_cov,~]  = Cov_Func(alpop,rs,theta0,Obstacle_Area,Covered_Area);
                if alpop_cov >= pop_cov
                    pop(i).Position=alpop;
                    pop(i).Cost=alpop_cov;
                else
                    L(i)=L(i)+1;
                end
            end
        end
    end
    clear alpop_cov alpop alpop2 k h phi;
    for i=1:nPop
        if pop(i).Cost>BestSol.Cost
            BestSol.Cost=pop(i).Cost;
            BestSol.Position=pop(i).Position;
        end 
    end
    BestCost(it)=BestSol.Cost;
    disp(['Iteration ' num2str(it) ': Best Cost = ' num2str(BestCost(it))]);

end
%% plot
alpop=BestSol.Position;
[coverage,Covered_Area] = Cov_Func(alpop(1:N,:),rs,theta0,Obstacle_Area,Covered_Area);
% show map
%imagesc(imresize(Obstacle_Area,[1000 1000],'bilinear'))
hold on;
%colorbar;
%{
for i = 1:1:numel(G.Edges.EndNodes)/2
    plot([pop(G.Edges.EndNodes(i,1)*2-1),pop(G.Edges.EndNodes(i,2)*2-1)],[pop(G.Edges.EndNodes(i,1)*2),pop(G.Edges.EndNodes(i,2)*2)],'Color','blue','linewidth',1);
end
%}
for i = 1:N
    plot (alpop(i,2) , alpop(i,1),'ro','MarkerSize', 3,'Color','red');
    text (alpop(i,2) , alpop(i,1), num2str(i),'FontSize',10,'Color','red');
end
clear i;
% show map
[obs_row, obs_col] = find(Obstacle_Area == 1);
plot(obs_row-1, obs_col-1,'.', 'MarkerSize', 1, 'Color', 'blue');
[obs_row, obs_col] = find(Covered_Area == 1/2);
plot(obs_row-1, obs_col-1,'.', 'MarkerSize', 1, 'Color', 'green');
[obs_row, obs_col] = find(Covered_Area == 1);
plot(obs_row-1, obs_col-1,'.', 'MarkerSize', 2, 'Color', 'red');
%[discovered_obs_row, discovered_obs_col] = find(Covered_Area == -1);                    % show discovered map
%plot(discovered_obs_row-1, discovered_obs_col-1,'.', 'MarkerSize', 20, 'Color', 'red');
axis equal;
xlim([0 size(Obstacle_Area,1)]);
ylim([0 size(Obstacle_Area,2)]);
title(['Weighted Coverage ratio:' num2str(BestSol.Cost*100) '% ' ]);
grid on;
drawnow;
%{
if save_fig==1
    saveas(gcf, ['Iteration' num2str(it) '.pdf' ]); % Lưu dưới dạng PDF
end
%}
%%
%save(name)


