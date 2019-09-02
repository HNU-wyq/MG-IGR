Xi = [0,0,1,1];
Xi = [Xi(1) Xi Xi(end)];
Eta = [0,0,1,1];
Eta = [Eta(1) Eta Eta(end)];
Zeta = [0,0,0,1,2,2,2];
Zeta = [Zeta(1) Zeta Zeta(end)];
k_ = length(Xi);
l_ = length(Eta);
m_ = length(Zeta);

B = cell(2,3,5);s=0.7;
B{1,1,1}=[0,       0,  0,  1];
B{1,1,2}=[-0.4,  0.3,  0,  s];
B{1,1,3}=[-1,    0.3,  0,  s];
B{1,1,4}=[-1.3,  -0.3, 0,  s+0.1];
B{1,1,5}=[0,      -1.3,  0,  1];

B{1,2,1}=[0.1,       0,  0,  1];
B{1,2,2}=[0.1,     -0.3,  0,  s];
B{1,2,3}=[0.1,     -0.6,  0,  s];
B{1,2,4}=[0.1,      -0.8, 0,  s+0.1];
B{1,2,5}=[0.1,      -1.3,  0,  s+0.1];

B{1,3,1}=[0.2,      0,  0,  1];
B{1,3,2}=[0.7,  0.3,  0,  s];
B{1,3,3}=[1.2,    0.3,  0,  s];
B{1,3,4}=[1.5,  -0.3, 0,  s+0.1];
B{1,3,5}=[0.2,      -1.3,  0,  1];

B{2,1,1}=[0,       0,  0.5,  1];
B{2,1,2}=[-0.5,  0.3,  0.5,  s];
B{2,1,3}=[-1,    0.3,  0.5,  s];
B{2,1,4}=[-1.3,  -0.3, 0.5,  s+0.1];
B{2,1,5}=[0,      -1.3,  0.5,  1];

B{2,2,1}=[0.1,       0,  0.5,  1];
B{2,2,2}=[0.1,     -0.3,  0.5,  s];
B{2,2,3}=[0.1,     -0.6,  0.5,  s];
B{2,2,4}=[0.1,      -0.8, 0.5,  s+0.1];
B{2,2,5}=[0.1,      -1.3,  0.5,  s+0.1];

B{2,3,1}=[0.2,      0,  0.5,  1];
B{2,3,2}=[0.7,  0.3,  0.5,  s];
B{2,3,3}=[1.2,    0.3,  0.5,  s];
B{2,3,4}=[1.5,  -0.3, 0.5,  s+0.1];
B{2,3,5}=[0.2,      -1.3,  0.5,  1];

B{3,1,1}=[0,       0,  1,  1];
B{3,1,2}=[-0.5,  0.3,  1,  s];
B{3,1,3}=[-1,    0.3,  1,  s];
B{3,1,4}=[-1.3,  -0.3, 1,  s+0.1];
B{3,1,5}=[0,      -1.3,  1,  1];

B{3,2,1}=[0.1,       0,  1,  1];
B{3,2,2}=[0.1,     -0.3,  1,  s];
B{3,2,3}=[0.1,     -0.6,  1,  s];
B{3,2,4}=[0.1,      -0.8, 1,  s+0.1];
B{3,2,5}=[0.1,      -1.3,  1,  s+0.1];

B{3,3,1}=[0.2,      0,  1,  1];
B{3,3,2}=[0.7,  0.3,  1,  s];
B{3,3,3}=[1.2,    0.3,  1,  s];
B{3,3,4}=[1.5,  -0.3, 1,  s+0.1];
B{3,3,5}=[0.2,      -1.3,  1,  1];


n = size(B,1);
m = size(B,2);
l = size(B,3);

% Order of basis
p = k_-n-1;
q = l_-m-1;
r = m_ -l -1;