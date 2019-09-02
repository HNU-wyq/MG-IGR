Xi = [0,0,1,1,2,2,3,3,4,4];
Xi = [Xi(1) Xi Xi(end)];
Eta = [0,1];
Eta = [Eta(1) Eta Eta(end)];
Zeta = [0,0,0,1,1,2,3,4,4,5,5,5];
Zeta = [Zeta(1) Zeta Zeta(end)];
k_ = length(Xi);
l_ = length(Eta);
m_ = length(Zeta);


B = cell(9,2,10);


s = 1/sqrt(2);
rr = 0.5; dy=0.4; du=0.04;
t = .01/2; %t/2
B{1,1,1} = [rr-t    0     0       1];
B{1,1,2} = [rr-t+du  0   -du      1];
B{1,1,3} = [rr-t+du  0   -dy+du   1];
B{1,1,4} = [rr-t    0   -dy       1];
B{1,1,5} = [rr-t-du  0   -dy-du   1];
B{1,1,6} = [rr-t-du  0   -2*dy+du 1];
B{1,1,7} = [rr-t    0   -2*dy     1];
B{1,1,8} = [rr-t+du  0   -2*dy-du 1];
B{1,1,9} = [rr-t+du  0   -3*dy+du 1];
B{1,1,10}= [rr-t    0   -1.2      1];

B{1,2,1} = [rr+t    0    0        1];
B{1,2,2} = [rr+t+du  0   -du      1];
B{1,2,3} = [rr+t+du  0   -dy+du   1];
B{1,2,4} = [rr+t    0    -dy      1];
B{1,2,5} = [rr+t-du  0   -dy-du   1];
B{1,2,6} = [rr+t-du  0   -2*dy+du 1];
B{1,2,7} = [rr+t    0    -2*dy    1];
B{1,2,8} = [rr+t+du  0   -2*dy-du 1];
B{1,2,9} = [rr+t+du  0   -3*dy+du 1];
B{1,2,10}= [rr+t    0     -1.2    1];


B{2,1,1} = [rr-t    rr-t   0        s];
B{2,1,2} = [rr-t+du  rr-t+du   -du      s];
B{2,1,3} = [rr-t+du  rr-t+du    -dy+du   s];
B{2,1,4} = [rr-t    rr-t    -dy      s];
B{2,1,5} = [rr-t-du  rr-t-du    -dy-du   s];
B{2,1,6} = [rr-t-du  rr-t-du    -2*dy+du s];
B{2,1,7} = [rr-t    rr-t    -2*dy    s];
B{2,1,8} = [rr-t+du  rr-t+du    -2*dy-du s];
B{2,1,9} = [rr-t+du  rr-t+du    -3*dy+du s];
B{2,1,10}= [rr-t    rr-t    -1.2    s];

B{2,2,1} = [rr+t    rr+t   0        s];
B{2,2,2} = [rr+t+du  rr+t+du    -du      s];
B{2,2,3} = [rr+t+du  rr+t+du    -dy+du   s];
B{2,2,4} = [rr+t    rr+t    -dy      s];
B{2,2,5} = [rr+t-du  rr+t-du    -dy-du   s];
B{2,2,6} = [rr+t-du  rr+t-du    -2*dy+du s];
B{2,2,7} = [rr+t    rr+t    -2*dy    s];
B{2,2,8} = [rr+t+du  rr+t+du    -2*dy-du s];
B{2,2,9} = [rr+t+du  rr+t+du    -3*dy+du s];
B{2,2,10}= [rr+t    rr+t    -1.2    s];

B{3,1,1} = [0      2*rr-t      0        1];
B{3,1,2} = [0      2*rr-t+du   -du      1];
B{3,1,3} = [0      2*rr-t+du   -dy+du   1];
B{3,1,4} = [0     2* rr-t     -dy      1];
B{3,1,5} = [0     rr-t-du   -dy-du   1];
B{3,1,6} = [0      rr-t-du   -2*dy+du 1];
B{3,1,7} = [0     2* rr-t     -2*dy    1];
B{3,1,8} = [0     2* rr-t+du   -2*dy-du 1];
B{3,1,9} = [0     2* rr-t+du   -3*dy+du 1];
B{3,1,10}= [0     2* rr-t     -1.2    1];

B{3,2,1} = [0    2*rr+t   0        1];
B{3,2,2} = [0  2*rr+t+du    -du      1];
B{3,2,3} = [0  2*rr+t+du    -dy+du   1];
B{3,2,4} = [0    2*rr+t    -dy      1];
B{3,2,5} = [0  rr+t-du    -dy-du   1];
B{3,2,6} = [0  rr+t-du    -2*dy+du 1];
B{3,2,7} = [0    2*rr+t    -2*dy    1];
B{3,2,8} = [0  2*rr+t+du    -2*dy-du 1];
B{3,2,9} = [0  2*rr+t+du    -3*dy+du 1];
B{3,2,10}= [0   2* rr+t    -1.2   1];

B{4,1,1} = [-(rr-t)    rr-t   0        s];
B{4,1,2} = [-(rr-t+du)  rr-t+du   -du      s];
B{4,1,3} = [-(rr-t+du)  rr-t+du    -dy+du   s];
B{4,1,4} = [-(rr-t)    rr-t    -dy      s];
B{4,1,5} = [-(rr-t-du)  rr-t-du    -dy-du   s];
B{4,1,6} = [-(rr-t-du)  rr-t-du    -2*dy+du s];
B{4,1,7} = [-(rr-t)    rr-t    -2*dy    s];
B{4,1,8} = [-(rr-t+du)  rr-t+du    -2*dy-du s];
B{4,1,9} = [-(rr-t+du)  rr-t+du    -3*dy+du s];
B{4,1,10}= [-(rr-t)    rr-t    -1.2    s];

B{4,2,1} = [-(rr+t)    rr+t   0        s];
B{4,2,2} = [-(rr+t+du)  rr+t+du   -du      s];
B{4,2,3} = [-(rr+t+du)  rr+t+du    -dy+du   s];
B{4,2,4} = [-(rr+t)    rr+t    -dy      s];
B{4,2,5} = [-(rr+t-du)  rr+t-du    -dy-du   s];
B{4,2,6} = [-(rr+t-du)  rr+t-du    -2*dy+du s];
B{4,2,7} = [-(rr+t)    rr+t    -2*dy    s];
B{4,2,8} = [-(rr+t+du)  rr+t+du    -2*dy-du s];
B{4,2,9} = [-(rr+t+du)  rr+t+du    -3*dy+du s];
B{4,2,10}= [-(rr+t)    rr+t    -1.2    s];

B{5,1,1} = [-(rr-t)    0   0        1];
B{5,1,2} = [-(rr-t+du)  0   -du      1];
B{5,1,3} = [-(rr-t+du)  0   -dy+du   1];
B{5,1,4} = [-(rr-t)    0   -dy      1];
B{5,1,5} = [-(rr-t-du)  0   -dy-du   1];
B{5,1,6} = [-(rr-t-du)  0   -2*dy+du 1];
B{5,1,7} = [-(rr-t)    0   -2*dy    1];
B{5,1,8} = [-(rr-t+du)  0   -2*dy-du 1];
B{5,1,9} = [-(rr-t+du)  0   -3*dy+du 1];
B{5,1,10}= [-(rr-t)    0   -1.2    1];

B{5,2,1} = [-(rr+t)    0   0        1];
B{5,2,2} = [-(rr+t+du)  0   -du      1];
B{5,2,3} = [-(rr+t+du)  0   -dy+du   1];
B{5,2,4} = [-(rr+t)    0   -dy      1];
B{5,2,5} = [-(rr+t-du)  0   -dy-du  1];
B{5,2,6} = [-(rr+t-du)  0   -2*dy+du 1];
B{5,2,7} = [-(rr+t)    0   -2*dy    1];
B{5,2,8} = [-(rr+t+du)  0   -2*dy-du 1];
B{5,2,9} = [-(rr+t+du)  0   -3*dy+du 1];
B{5,2,10}= [-(rr+t)    0   -1.2    1];


B{6,1,1} =[-(rr-t)    -(rr-t)   0        s];
B{6,1,2} = [-(rr-t+du)  -(rr-t+du)   -du      s];
B{6,1,3} = [-(rr-t+du)  -(rr-t+du)   -dy+du   s];
B{6,1,4} = [-(rr-t)    -(rr-t)   -dy      s];
B{6,1,5} = [-(rr-t-du)  -(rr-t-du)   -dy-du   s];
B{6,1,6} = [-(rr-t-du)  -(rr-t-du)   -2*dy+du s];
B{6,1,7} = [-(rr-t)    -(rr-t)   -2*dy    s];
B{6,1,8} = [-(rr-t+du)  -(rr-t+du)   -2*dy-du s];
B{6,1,9} = [-(rr-t+du)  -(rr-t+du)   -3*dy+du s];
B{6,1,10}= [-(rr-t)    -(rr-t)   -1.2   s];

B{6,2,1} = [-(rr+t)    -(rr+t)   0       s];
B{6,2,2} = [-(rr+t+du)  -(rr+t+du)   -du      s];
B{6,2,3} = [-(rr+t+du)  -(rr+t+du)   -dy+du   s];
B{6,2,4} = [-(rr+t)    -(rr+t)   -dy      s];
B{6,2,5} = [-(rr+t-du)  -(rr+t-du)   -dy-du   s];
B{6,2,6} = [-(rr+t-du)  -(rr+t-du)   -2*dy+du s];
B{6,2,7} = [-(rr+t)    -(rr+t)   -2*dy    s];
B{6,2,8} = [-(rr+t+du)  -(rr+t+du)   -2*dy-du s];
B{6,2,9} = [-(rr+t+du)  -(rr+t+du)   -3*dy+du s];
B{6,2,10}= [-(rr+t)    -(rr+t)   -1.2   s];

B{7,1,1} = [0    -(rr-t)   0        1];
B{7,1,2} = [0  -(rr-t+du)   -du      1];
B{7,1,3} = [0  -(rr-t+du)   -dy+du   1];
B{7,1,4} = [0    -(rr-t)   -dy      1];
B{7,1,5} = [0  -(rr-t-du)   -dy-du   1];
B{7,1,6} = [0  -(rr-t-du)   -2*dy+du 1];
B{7,1,7} = [0    -(rr-t)   -2*dy    1];
B{7,1,8} = [0  -(rr-t+du)   -2*dy-du 1];
B{7,1,9} = [0  -(rr-t+du)   -3*dy+du 1];
B{7,1,10}= [0    -(rr-t)   -1.2    1];

B{7,2,1} = [0    -(rr+t)   0        1];
B{7,2,2} = [0  -(rr+t+du)   -du      1];
B{7,2,3} = [0  -(rr+t+du)   -dy+du   1];
B{7,2,4} = [0    -(rr+t)   -dy      1];
B{7,2,5} = [0  -(rr+t-du)   -dy-du   1];
B{7,2,6} = [0  -(rr+t-du)   -2*dy+du 1];
B{7,2,7} = [0    -(rr+t)   -2*dy    1];
B{7,2,8} = [0  -(rr+t+du)   -2*dy-du 1];
B{7,2,9} = [0  -(rr+t+du)   -3*dy+du 1];
B{7,2,10}= [0    -(rr+t)   -1.2    1];

B{8,1,1} = [(rr-t)    -(rr-t)   0        s];
B{8,1,2} = [rr-t+du  -(rr-t+du)   -du      s];
B{8,1,3} = [rr-t+du  -(rr-t+du)    -dy+du   s];
B{8,1,4} = [rr-t    -(rr-t)    -dy      s];
B{8,1,5} = [rr-t-du  -(rr-t-du)    -dy-du   s];
B{8,1,6} = [rr-t-du  -(rr-t-du)    -2*dy+du s];
B{8,1,7} = [rr-t    -(rr-t)    -2*dy    s];
B{8,1,8} = [rr-t+du  -(rr-t+du)    -2*dy-du s];
B{8,1,9} = [rr-t+du  -(rr-t+du)    -3*dy+du s];
B{8,1,10}= [rr-t    -(rr-t)    -1.2    s];

B{8,2,1} = [rr+t    -(rr+t)   0        s];
B{8,2,2} = [rr+t+du  -(rr+t+du)    -du      s];
B{8,2,3} = [rr+t+du  -(rr+t+du)    -dy+du   s];
B{8,2,4} = [rr+t    -(rr+t)    -dy      s];
B{8,2,5} = [rr+t-du  -(rr+t-du)    -dy-du   s];
B{8,2,6} = [rr+t-du  -(rr+t-du)    -2*dy+du s];
B{8,2,7} = [rr+t    -(rr+t)    -2*dy    s];
B{8,2,8} = [rr+t+du  -(rr+t+du)    -2*dy-du s];
B{8,2,9} = [rr+t+du  -(rr+t+du)    -3*dy+du s];
B{8,2,10}= [rr+t    -(rr+t)    -1.2   s];

B{9,1,1} =[rr-t    0   0        1];
B{9,1,2} = [rr-t+du  0   -du      1];
B{9,1,3} = [rr-t+du  0   -dy+du   1];
B{9,1,4} = [rr-t    0   -dy      1];
B{9,1,5} = [rr-t-du  0   -dy-du   1];
B{9,1,6} = [rr-t-du  0   -2*dy+du 1];
B{9,1,7} = [rr-t    0   -2*dy    1];
B{9,1,8} = [rr-t+du  0   -2*dy-du 1];
B{9,1,9} = [rr-t+du  0   -3*dy+du 1];
B{9,1,10}= [rr-t    0   -1.2   1];

B{9,2,1} = [rr+t    0   0        1];
B{9,2,2} = [rr+t+du  0   -du      1];
B{9,2,3} = [rr+t+du  0   -dy+du   1];
B{9,2,4} = [rr+t    0   -dy      1];
B{9,2,5} = [rr+t-du  0   -dy-du   1];
B{9,2,6} = [rr+t-du  0   -2*dy+du 1];
B{9,2,7} = [rr+t    0   -2*dy    1];
B{9,2,8} = [rr+t+du  0   -2*dy-du 1];
B{9,2,9} = [rr+t+du  0   -3*dy+du 1];
B{9,2,10}= [rr+t    0   -1.2    1];





%Ğı×ª
% th=pi/3;
% R = [ 1 0 0 0; 0 cos(th) -sin(th) 0; 0 sin(th) cos(th) 0; 0 0 0 1];
% 
% for j = 1 : size(B,2)
%     for k = 1 : size(B,3)
%         B{2,j,k} = (R*B{2,j,k}')';
%     end
% end
% 
% for j = 1 : size(B,2)
%     for k = 1 : size(B,3)
%         B{3,j,k} = (R*R*B{3,j,k}')';
%     end
% end
% 
% for j = 1 : size(B,2)
%     for k = 1 : size(B,3)
%         B{4,j,k} = (R*R*R*B{4,j,k}')';
%     end
% end

% Number of weights
n = size(B,1);
m = size(B,2);
l = size(B,3);

% Order of basis
p = k_-n-1;
q = l_-m-1;
r = m_ -l -1;
