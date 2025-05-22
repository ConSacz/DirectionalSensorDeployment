function Obstacle_Area = gen_target_area(nT)

Obstacle_Area = zeros(100,100);
idx = randperm(100*100, nT);
Obstacle_Area(idx)=1;