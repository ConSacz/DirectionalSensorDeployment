function [coverage,Covered_Area] = Cov_Func_v2(pop,rs,theta0,Obstacle_Area,Covered_Area)
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
    x0 = pop(j,1); % Tọa độ tâm
    y0 = pop(j,2);

    % Tạo lưới điểm
    [X, Y] = meshgrid(linspace(1, size(Covered_Area,1),100), linspace(1, size(Covered_Area,1),100));
    
    % Thông số cung tròn
    alpha = pop(j,3); % Góc hướng bắt đầu (radian)
    
    % Tính khoảng cách từ mỗi điểm đến tâm
    D = sqrt((X - x0).^2 + (Y - y0).^2);
    
    % Tính góc của mỗi điểm so với tâm
    Theta = atan2(Y - y0, X - x0); 
    
    % Chuyển góc về khoảng từ 0 đến 2*pi nếu cần
    Theta(Theta < 0) = Theta(Theta < 0) + 2*pi;
    
    % Kiểm tra điểm nào thuộc cung
    
    in_circle = D <= rs; % Điều kiện 1: Nằm trong bán kính
    if alpha - theta0/2 < 0
        in_angle = (Theta >= alpha - theta0/2 +2*pi) | (Theta <= alpha +theta0/2); 
    elseif alpha + theta0/2 > 2*pi 
        in_angle = (Theta >= alpha - theta0/2) | (Theta <= alpha + theta0/2 - 2*pi);
    else
        in_angle = (Theta >= alpha - theta0/2) & (Theta <= alpha + theta0/2); % Điều kiện 2: Nằm trong góc cung
    end

    inside_sector = inside_sector | (in_circle & in_angle); % Kết hợp hai điều kiện
    %[x_cover, y_cover] = find(inside_sector==1);
    %Covered_Area (X(x_cover,y_cover),Y(x_cover,y_cover)) = Obstacle_Area (X(x_cover,y_cover),Y(x_cover,y_cover));
    Covered_Area = inside_sector.* Obstacle_Area;
    %clear D Theta in_circle in_angle inside_sector;
end

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
