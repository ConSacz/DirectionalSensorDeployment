function [coverage,Covered_Area] = Cov_Func_v1(pop,rs,theta0,Obstacle_Area,Covered_Area)
%This is fitness function to cal random area coverage ratio

%pop is a 1x2Ndim matrix holding position of nodes
%rs is the sensing rad of nodes                           
%Obstacle_Area= Area1;                            
%Covered_Area = zeros(size(Area1,1),size(Area1,2));

%% recover sensor covered area
[obs_row, obs_col] = find(Covered_Area ~= 0);
for i = 1:numel(obs_col)
    Covered_Area(obs_row(i), obs_col(i)) = 0;
end

%% check sensor covered area
inside_sector = false(size(Covered_Area,1),size(Covered_Area,1));
for j=1:(size(pop,1))
    %%
    % Node position
    x0 = pop(j,1);
    y0 = pop(j,2);

    % Boundary constraint
    x_ub=min(ceil(x0+rs),size(Covered_Area,1));
    x_lb=max(floor(x0-rs),1);
    y_ub=min(ceil(y0+rs),size(Covered_Area,1));
    y_lb=max(floor(y0-rs),1);

    % Local Grid
    [X, Y] = meshgrid(linspace(x_lb, x_ub, x_ub-x_lb+1), linspace(y_lb, y_ub, y_ub-y_lb+1));
    
    % node angle direction
    alpha = pop(j,3); 
    
    % Distance matrix
    D = sqrt((X - x0).^2 + (Y - y0).^2);
    
    % Angle matrix
    Theta = atan2(Y - y0, X - x0); 
    
    % Boundary constraint
    Theta(Theta < 0) = Theta(Theta < 0) + 2*pi;
    
    % In rs condition
    in_circle = D <= rs;

    % Theta in theta0 condition
    if alpha - theta0/2 < 0
        in_angle = (Theta >= alpha - theta0/2 +2*pi) | (Theta <= alpha +theta0/2); 
    elseif alpha + theta0/2 > 2*pi 
        in_angle = (Theta >= alpha - theta0/2) | (Theta <= alpha + theta0/2 - 2*pi);
    else
        in_angle = (Theta >= alpha - theta0/2) & (Theta <= alpha + theta0/2); 
    end
    
    %both conditions
    inside_sector(y_lb:y_ub,x_lb:x_ub) = inside_sector(y_lb:y_ub,x_lb:x_ub) | (in_circle & in_angle); 
    
end       
Covered_Area = inside_sector.* Obstacle_Area;
    %clear D Theta in_circle in_angle inside_sector;


%% add obstacle to covered area
[obs_row, obs_col] = find(Obstacle_Area == 0);
for i = 1:numel(obs_col)
    if Covered_Area (obs_row(i), obs_col(i)) == 1
        Covered_Area(obs_row(i), obs_col(i)) = -2;
    end
end

count1=numel(find(Covered_Area == 1));		                           % count covered points on wanted location  (wanted)
count2=numel(find(Covered_Area == -2));		                           % count covered points on unwanted location (obstacles)
count3=numel(Obstacle_Area)-numel(find(Obstacle_Area ~= 1));           % count total points on wanted location
coverage=((count1-count2)/count3);	                            % function to avoid obstacles
%coverage=(count1/count3);		                                % function to aim on wanted area

%% recover obs covered area

[obs_row, obs_col] = find(Covered_Area == -2);
for i = 1:numel(obs_col)
    Covered_Area(obs_row(i), obs_col(i)) = -1;
end	        % function to aim on wanted area
%%
%figure; hold on;
%scatter(X(:), Y(:), 5, 'k', 'filled'); % Vẽ toàn bộ grid
%scatter(X(inside_sector), Y(inside_sector), 5, 'r', 'filled'); % Tô màu đỏ điểm nằm trong cung
%plot(x_arc, y_arc, 'b', 'LineWidth', 2); % Vẽ cung tròn
%plot(x0, y0, 'ro', 'MarkerFaceColor', 'r'); % Vẽ tâm cung
