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

function [y1,y2]=Crossover(x1,x2,VarMax,VarMin)
 js=0;
%    while 1
%     alpha=rand(size(x1));
%     y1=alpha.*x1+(1-alpha).*x2;
%     y2=alpha.*x2+(1-alpha).*x1;
%     
% %       if all ((1.18<=y1(1)<=1.26)&&(1.9<=y1(2)<=2.3))  %add the judge sentences to select the true parameter matrix
% %         if all (all(1.18<=y2(1)<=1.26)&&(1.9<=y2(2)<=2.3))  %add the judge sentences to select the true parameter matrix
% %           break
% %         end
% %       end
%      %% screwdriver 3组设计变量
%       if all(VarMin<=y1)&&all(y1<=VarMax)  %add the judge sentences to select the true parameter matrix
%         if all (VarMin<=y2)&&all(y2<=VarMax) %add the judge sentences to select the true parameter matrix
%           break
%         end
%       end
%    end 
    alpha=rand(size(x1));
    y1=alpha.*x1+(1-alpha).*x2;
    y2=alpha.*x2+(1-alpha).*x1;
     while 1
         js=js+1
      a1=find(VarMin<=y1);
      b1=find(y1<=VarMax);
      a2=find(VarMin<=y2);
      b2=find(y2<=VarMax);
      ab1=intersect(a1,b1);
      ab2=intersect(a2,b2);
      c=(1:10);

      buab=union((setdiff(c,ab1)),(setdiff(c,ab2)));
      if numel(buab)==0
          break
      end
      alpha1=rand(size(buab));
      y1(buab)=alpha1.*x1(buab)+(1-alpha1).*x2(buab);
      y2(buab)=alpha1.*x2(buab)+(1-alpha1).*x1(buab);
     end
    
         
end