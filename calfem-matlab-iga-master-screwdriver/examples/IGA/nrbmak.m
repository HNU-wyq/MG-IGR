function nurbs = nrbmak(coefs,knots)
%
% NRBMAK: Construct the NURBS volume given the control points
%            and the knots.
% 
nurbs.form   = 'B-NURBS';
nurbs.dim    = 4;
np = size(coefs)
dim=4;
uorder = size(knots{1},2)-np(2)
vorder = size(knots{2},2)-np(3)
 worder = size(knots{3},2)-np(4)
 
  nurbs.knots=konts;
  nurbs.order=[uorder vorder worder];