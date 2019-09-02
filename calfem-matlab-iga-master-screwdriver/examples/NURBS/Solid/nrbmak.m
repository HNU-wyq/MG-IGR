function nurbs = nrbmak(coefs,knots)
%
% NRBMAK: Construct the NURBS volume given the control points
%            and the knots.
% 
nurbs.form   = 'B-NURBS';
nurbs.dim    = 4;
np = size(coefs)
nurbs.number =np;
dim=4;
uorder = size(knots{1},2)-np(1)
vorder = size(knots{2},2)-np(2)
worder = size(knots{3},2)-np(3)
  nurbs.coefs=coefs;
  nurbs.knots=knots;
  nurbs.order=[uorder vorder worder];