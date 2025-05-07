%% Directional Sensor
%% DEPLOYMENT

%%
clc;
clear;

%% 
close all;
%% Network parameter
for N=30:10:50
for trial=1:20

% Monitor area
Covered_Area = zeros(100,100);
%Obstacle_Area = gen_random_distribution_area();
%Obstacle_Area = ones(100,100);
%Obstacle_Area=genarea();
Obstacle_Area = gen_target_area(1000);

% nodes info
MaxIt = 200;              % Maximum Number of Iterations
a = 1;                    % Acceleration Coefficient Upper Bound
%N = 50;
nPop=25;
rc = 20;
name=['./data/hetero_target/',num2str(N),' nodes/case4/hetero_',num2str(trial),'.mat'];
%name=['./data/comparison/target/ABC/ABC',num2str(N),'nodes_',num2str(trial),'.mat'];

% homogenerous sensor 
% case1
% rs=ones(1,N)*12;
% theta0=ones(1,N)*(pi);


%heterogenerous sensor
% case2
% rs = ones(1,N)*12;
% rs(N-9:N-5)=15;
% rs(N-4:N)=8;
% theta0=ones(1,N)*pi;
% theta0(N-9:N-5)=2*pi/3;
% theta0(N-4:N)=3*pi/2;

% case3
% v=floor(N/3);
% rs = ones(1,N)*12;
% rs(N-(2*v-1):N-v)=15;
% rs(N-(v-1):N)=8;
% theta0=ones(1,N)*pi;
% theta0(N-(2*v-1):N-v)=2*pi/3;
% theta0(N-(v-1):N)=3*pi/2;

% case4
rs = unifrnd(10,15,[1 N]);
theta0=unifrnd(pi,3*pi/2,[1 N]);

sink=[50 50];
Scout_bee=nPop;
Onlooker_bee=nPop;

% bee parameter
a=1;


%% Init first pop
empty_individual.Position=[]; % Empty solution
empty_individual.Cost=[];     % Empty cost function of that solution
BestSol.Position=[];	% Best Solution
BestSol.Cost=0;	% Best solution function value
pop=repmat(empty_individual,nPop,1); % pop includes nPop solutions
BestCost=zeros(MaxIt,1); % Best solution every iteration
L=zeros(1,nPop);

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
        phi=a*unifrnd(-1, +1, [N 3])*(1-L(i)/MaxIt)^5;
        alpop = pop(i).Position + phi.*(pop(i).Position-pop(k).Position);
        alpop(:,1:2) = min(max(alpop(:,1:2), 1),size(Obstacle_Area,1)-1);
        alpop(:,3)=mod(alpop(:,3),2*pi);
        if Connectivity_graph(Graph(alpop(:,1:2),rc),[]) == 1
            %[pop_cov,~] = Cov_Func(pop(i).Position,rs,theta0,Obstacle_Area,Covered_Area);
            [alpop_cov,~]  = Cov_Func(alpop,rs,theta0,Obstacle_Area,Covered_Area);
            if alpop_cov >= pop(i).Cost
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
        % randomly choose a pop, prioritize that have high coverage
        if all(isnan(E_ave))
            i=randi([1,nPop]);
        else
            i=randsrc(1,1,[1:nPop;E_ave]);
        end
        % change the position and angle of each node
        for k=1:N
            alpop=pop(i).Position;
            h=randi([1,N]);
            phi=a*unifrnd(-1, +1, [1 3])*(1-L(i)/MaxIt)^2;
            alpop(k,1:3)  = pop(i).Position(k,1:3) + phi.*(pop(i).Position(k,1:3)-pop(i).Position(h,1:3));
            alpop(:,1:2) = min(max(alpop(:,1:2), 1),size(Obstacle_Area,1)-1);
            alpop(:,3)=mod(alpop(:,3),2*pi);
            if Connectivity_graph(Graph(alpop(:,1:2),rc),[]) == 1
                %[pop_cov,~] = Cov_Func(pop(i).Position,rs,theta0,Obstacle_Area,Covered_Area);
                [alpop_cov,~]  = Cov_Func(alpop,rs,theta0,Obstacle_Area,Covered_Area);
                if alpop_cov >= pop(i).Cost
                    pop(i).Position=alpop;
                    pop(i).Cost=alpop_cov;
%                 else
%                     L(i)=L(i)+1;
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
    clear h j k phi alpop alpop2 alpop_cov E E_ave;
    if BestCost(it)==1
        break;
    end    
end
%% plot
% % %alpop=pop(4).Position;
alpop=BestSol.Position;
[coverage,Covered_Area] = Cov_Func(alpop,rs,theta0,Obstacle_Area,Covered_Area);
% show map
%imagesc(imresize(Obstacle_Area,[100 100],'bilinear'))
% figure;
% plot(BestCost);
figure;
hold on;

% show map

%[obs_row, obs_col] = find(Covered_Area == 0);
%plot(obs_row-1, obs_col-1,'.', 'MarkerSize', 2, 'Color', 'green');


%[discovered_obs_row, discovered_obs_col] = find(Covered_Area == -1);                    % show discovered map
%plot(discovered_obs_row-1, discovered_obs_col-1,'.', 'MarkerSize', 20, 'Color', 'red');
%colorbar;
%{
for i = 1:1:numel(G.Edges.EndNodes)/2
    plot([pop(G.Edges.EndNodes(i,1)*2-1),pop(G.Edges.EndNodes(i,2)*2-1)],[pop(G.Edges.EndNodes(i,1)*2),pop(G.Edges.EndNodes(i,2)*2)],'Color','blue','linewidth',1);
end
%}
for i = 1:N
    % i-th node
    plot (alpop(i,2) , alpop(i,1),'ro','MarkerSize', 3,'Color','red');
    text (alpop(i,2) , alpop(i,1), num2str(i),'FontSize',10,'Color','red');
    theta = linspace(alpop(i,3)-theta0(i)/2 , alpop(i,3) + theta0(i)/2 , 100); 
    x = alpop(i,1) + rs(i) * cos(theta);
    y = alpop(i,2) + rs(i) * sin(theta);
    x_fill = [alpop(i,1), x, alpop(i,1)]; 
    y_fill = [alpop(i,2), y, alpop(i,2)];
    fill(y_fill, x_fill, 'g', 'FaceAlpha', 0.3, 'EdgeColor', 'g');

end
[obs_row, obs_col] = find(Obstacle_Area == 1);
plot(obs_row, obs_col,'.', 'MarkerSize', 5, 'Color', 'red');
[obs_row, obs_col] = find(Covered_Area == 1);
plot(obs_row, obs_col,'.', 'MarkerSize', 5, 'Color', 'blue');
clear i x_fill y_fill theta obs_col obs_row;

axis equal;
xlim([0 size(Obstacle_Area,1)]);
ylim([0 size(Obstacle_Area,2)]);
title(['Coverage ratio: ' num2str(BestSol.Cost*100) '% ' ]);
grid on;
drawnow;

% save fig
exportgraphics(gcf, './data/zFigures/deployment figures/hetero_target_case4.pdf', 'ContentType', 'vector', 'BackgroundColor', 'none');

%%
save(name)
end
end
