% 
% Copyright (c) 2016, Yarpiz (www.yarpiz.com)
% All rights reserved. Please read the "license.txt" for license terms.
% 
% Project Code: YPEA126
% Project Title: Non-dominated Sorting Genetic Algorithm III (NSGA-III)
% Publisher: Yarpiz (www.yarpiz.com)
% 
% Implemented by: S. Mostapha Kalami Heris, PhD (member of Yarpiz Team)
% 
% Contact Info: sm.kalami@gmail.com, info@yarpiz.com
% 
% Base Reference Paper:
% K. Deb and H. Jain, "An Evolutionary Many-Objective Optimization Algorithm 
% Using Reference-Point-Based Nondominated Sorting Approach, Part I: Solving
% Problems With Box Constraints,"
% in IEEE Transactions on Evolutionary Computation,
% vol. 18, no. 4, pp. 577-601, Aug. 2014.
% 
% Reference Papaer URL: http://doi.org/10.1109/TEVC.2013.2281535
% 

function y=Mutate(x,mu,sigma,VarMax,VarMin)
   while 1
    nVar=numel(x);
    
    nMu=ceil(mu*nVar);

    j=randsample(nVar,nMu);
    
    y=x;
    
    y(j)=x(j)+sigma*randn(size(j));
%        if all ((1.18<=y(1)<=1.26)&&(1.9<=y(2)<=2.3))  %add the judge sentences to select the true parameter matrix        
%           break
%        end
      if  (all(VarMin(j)<=y(j))&&all(y(j)<=VarMax(j)))  %add the judge sentences to select the true parameter matrix        
          break
      end
   end

end