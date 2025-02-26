function [coverage,Covered_Area] = Cov_Func(pop,rs,theta0,Obstacle_Area,Covered_Area)
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
for j=1:(size(pop,1))
    % y value of sensor is pop(j,1)
    % x value of sensor is pop(j,2)
    start_point=[floor(pop(j,1)) floor(pop(j,2))];
    for i = (-rs-1):(rs+1)
        for k = (-rs-1):(rs+1)
            map_x= start_point(2)+1+i;
            map_y= start_point(1)+1+k;
            % case selected point in map
            if map_x>0 && map_x<size(Obstacle_Area,2) && map_y>0 && map_y<size(Obstacle_Area,1) 
                dist = sqrt((map_y-pop(j,1))^2+(map_x-pop(j,2))^2);
                theta=atan2(map_y-pop(j,1),map_x-pop(j,2));
                if theta < 0
                    theta = theta + 2*pi; % Đưa về khoảng [0, 2*pi]
                end
                % case dist <= rs and theta in theta0
                if all((dist <= rs) & (theta >= pop(j,3) - theta0/2) & (theta <= pop(j,3) + theta0/2))
                    Covered_Area(map_x,map_y) = Obstacle_Area(map_x,map_y);
                end
            end
        end
    end
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

