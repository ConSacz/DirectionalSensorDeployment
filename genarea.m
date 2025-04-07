function Obstacle_Area = genarea()
% Interest area value >0
% Obstacle area value =0

I_raw=imread('./maps/image6.png');          % map here
I_raw=imresize(I_raw,[100,100]);            % resize
I_raw=rgb2gray(I_raw);
I_raw=imgaussfilt(I_raw, 2);
Obstacle_Area=im2double(I_raw);
%Obstacle_Area=I_raw(:,:,1);
%Obstacle_Area=Obstacle_Area*255;
Obstacle_Area=1-Obstacle_Area;
Obstacle_Area(Obstacle_Area>=0.7)=1;
Obstacle_Area(Obstacle_Area<0.7)=0;

%[obs_row, obs_col] = find(Obstacle_Area == 0);
%figure;
%imshow(I_gray);
%figure;
%imagesc(Obstacle_Area);
%colorbar;