%%薄壁结构，长宽高（3,1,0.1）

p = 1;
q = 1;
r = 1;

%---参数坐标
Xi = [0,1];
Xi = [Xi(1) Xi Xi(end)];
Eta = [0,1];
Eta = [Eta(1) Eta Eta(end)];
Zeta = [0,1];
Zeta = [Zeta(1) Zeta Zeta(end)];
k_ = length(Xi);
l_ = length(Eta);
m_ = length(Zeta);


B = cell(2,2,2);

s = 1/sqrt(2);
%--控制点
B{1,1,1} = [0 0 0 1];
B{1,1,2} = [0 0 .1 1];
B{1,2,1} = [0 1 0 1];
B{1,2,2} = [0 1 .1 1];

B{2,1,1} = [3 0 0 1];
B{2,1,2} = [3 0 .1 1];
B{2,2,1} = [3 1 0 1];
B{2,2,2} = [3 1 .1 1];


% Number of weights
n = size(B,1);
m = size(B,2);
l = size(B,3);

% Order of basis
p = k_-n-1;
q = l_-m-1;
r = m_ -l -1;
