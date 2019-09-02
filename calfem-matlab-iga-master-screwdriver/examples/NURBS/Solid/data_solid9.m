p = 1;
q = 1;
r = 1;

% Xi = [0,0,1,1,2,2,3,3,4,4];
% Xi = [Xi(1) Xi Xi(end)];
% Eta = [0,1];
% Eta = [Eta(1) Eta Eta(end)];
% Zeta = [0,1];
% Zeta = [Zeta(1) Zeta Zeta(end)];
% k_ = length(Xi);
% l_ = length(Eta);
% m_ = length(Zeta);


B = cell(9,2,3);

s = 1/sqrt(2);
r = 2;
t = 1/2; %t/2

B{1,1,1} = [r-t 0 0 1];
B{1,1,2} = [r-t 0 2.5 1];
B{1,1,3} = [r-t 0 5 1];
B{1,2,1} = [r 0 0 1];
B{1,2,2} = [r 0 2.5 1];
B{1,2,3} = [r 0 5 1];

B{2,1,1} = [r-t r-t 0 s];
B{2,1,2} = [r-t r-t 2.5 s];
B{2,1,3} = [r-t r-t 5 s];
B{2,2,1} = [r r 0 s];
B{2,2,2} = [r r 2.5 s];
B{2,2,3} = [r r 5 s];

B{3,1,1} = [0 r-t 0 1];
B{3,1,2} = [0 r-t 2.5 1];
B{3,1,3} = [0 r-t 5 1];
B{3,2,1} = [0 r 0 1];
B{3,2,2} = [0 r 2.5 1];
B{3,2,3} = [0 r 5 1];

B{4,1,1} = [-(r-t) r-t 0 s];
B{4,1,2} = [-(r-t) r-t 2.5 s];
B{4,1,3} = [-(r-t) r-t 5 s];
B{4,2,1} = [-(r) r 0 s];
B{4,2,2} = [-(r) r 2.5 s];
B{4,2,3} = [-(r) r 5 s];

B{5,1,1} = [-(r-t) 0 0 1];
B{5,1,2} = [-(r-t) 0 2.5 1];
B{5,1,3} = [-(r-t) 0 5 1];
B{5,2,1} = [-(r) 0 0 1];
B{5,2,2} = [-(r) 0 2.5 1];
B{5,2,3} = [-(r) 0 5 1];

B{6,1,1} = [-(r-t) -(r-t) 0 s];
B{6,1,2} = [-(r-t) -(r-t) 2.5 s];
B{6,1,3} = [-(r-t) -(r-t) 5 s];
B{6,2,1} = [-(r) -(r) 0 s];
B{6,2,2} = [-(r) -(r) 2.5 s];
B{6,2,3} = [-(r) -(r) 5 s];

B{7,1,1} = [0 -(r-t) 0 1];
B{7,1,2} = [0 -(r-t) 2.5 1];
B{7,1,3} = [0 -(r-t) 5 1];
B{7,2,1} = [0 -(r) 0 1];
B{7,2,2} = [0 -(r) 2.5 1];
B{7,2,3} = [0 -(r) 5 1];

B{8,1,1} = [r-t -(r-t) 0 s];
B{8,1,2} = [r-t -(r-t) 2.5 s];
B{8,1,3} = [r-t -(r-t) 5 s];
B{8,2,1} = [r -(r) 0 s];
B{8,2,2} = [r -(r) 2.5 s];
B{8,2,3} = [r -(r) 5 s];

B{9,1,1} = [r-t 0 0 1];
B{9,1,2} = [r-t 0 2.5 1];
B{9,1,3} = [r-t 0 5 1];
B{9,2,1} = [r 0 0 1];
B{9,2,2} = [r 0 2.5 1];
B{9,2,3} = [r 0 5 1];
B=permute(B,[2 3 1]);
% Number of weights
n = size(B,1);
m = size(B,2);
l = size(B,3);

Zeta = [0,0,1,1,2,2,3,3,4,4];
Zeta= [Zeta(1) Zeta Zeta(end)];
Eta = [0,0,1,1];
Eta = [Eta(1) Eta Eta(end)];
Xi = [0,1];
Xi = [Xi(1) Xi Xi(end)];
k_ = length(Xi);
l_ = length(Eta);
m_ = length(Zeta);

% Order of basis
p = k_-n-1;
q = l_-m-1;
r = m_ -l -1;
