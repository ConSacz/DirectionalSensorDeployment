load("firsttry.mat");
for it = 1: MaxIt
    %% plot
    pop=popIt(it,:);
    
    clf();
    hold on;
    % show map
    imagesc(Obstacle_Area)
    colorbar;
    
    G=Graph(pop,rc);
    for i = 1:1:numel(G.Edges.EndNodes)/2
        plot([pop(G.Edges.EndNodes(i,2)*2),pop(G.Edges.EndNodes(i,1)*2)],[pop(G.Edges.EndNodes(i,2)*2-1),pop(G.Edges.EndNodes(i,1)*2-1)],'Color','green','linewidth',1);
    end

    for i = 1:2:numel(pop)
        plot (pop(1,i+1) , pop(1,i),'ro','MarkerSize', 3,'Color','red');
        hold on;
        viscircles ([pop(1,i+1) pop(1,i)],rs,'LineWidth',1,'Color', 'k');
        text (pop(1,i+1) , pop(1,i), num2str(i/2+0.5),'FontSize',8,'Color','red');
    end
    clear i;

    xlim([0 size(Obstacle_Area,1)]);
    ylim([0 size(Obstacle_Area,2)]);
    title([num2str(BestCostIt(it)*100) '%  at time step:  '  num2str(it)]);
    grid on;
    drawnow;
end