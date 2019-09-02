options = gaoptimset('display', 'iter', 'Generations', 5 ,'PopulationSize', 10, 'PlotFcns', @gaplotbestf, 'TolFun',1e-4);
lb=4.5;ub=5;
fid=@interiga;
[zast,cast,exitflag] = ga(fid, 6, [], [], [], [],lb, ub, [], options);
