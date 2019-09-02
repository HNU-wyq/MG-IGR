function PlotCosts1(pop,it)

%     Costs=[pop.Cost];
    
    plot(it,pop(:),'r*','MarkerSize',8);
    xlabel('iter');
    ylabel('2nd Objective');
    grid on;

end