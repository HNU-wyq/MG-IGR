clear all; close all; clc;

%% Load data   (3个自由度的solid单元)
addpath(genpath('C:\Users\Administrator\Desktop\IGA\calfem-matlab-iga-master'));
count=3;
run ./../NURBS/Solid/data_solid8
 %run ./../NURBS/solid/data_solid7

%判断二维或三维
if count==2
   deg.p = p; deg.q=q; clear p q
elseif count==3
   deg.p = p; deg.q=q; deg.r = r; clear p q r
end

% Perform 3 subdivisons in Xi, and 1 in  Eta & Zeta
[ X_ ] = getSubDivKVValues( Xi, 3 );
[B,Xi]=RefineKnotSolidXi( B,deg.p,Xi,X_ );
[ E_ ] = getSubDivKVValues( Eta, 1);
[B,Eta]=RefineKnotSolidEta( B,deg.q,Eta,E_ );
[ Z_ ] =getSubDivKVValues( Zeta, 1);
[B,Zeta]=RefineKnotSolidZeta( B,deg.r,Zeta,Z_ );

% Number of control points after refinement
n = size(B,1);
m = size(B,2);
l = size(B,3);

%% Build connectivity arrays
% Knot vectors in analysis
KV.Xi =Xi; KV.Eta = Eta; KV.Zeta = Zeta; clear Xi Eta Zeta

% Build connectivity arrays
nel = (n-deg.p) * (m-deg.q) * (l-deg.r);       % number of elements-三维
%nel = (n-deg.p) * (m-deg.q);

nnp = n*m*l; % number of global basis functions
nen = (deg.p+1)*(deg.q+1)*(deg.r+1);           % number of local basis functions
%nen = (deg.p+1)*(deg.q+1);

ndof = nnp*3; % number of global degrees of freedom
ldof = nen*3; % number of local degrees of freedom
% Build connectivity arrays
[INN,IEN] = BldINCIEN( deg.p,deg.q,deg.r,n,m,l ); % = INn (is unique for IGA), IEN(=ENOD', if there is a CALFEM counterpart it might be called ENOD, but transposed)
ID = reshape(1:max(max(IEN))*3,3,max(max(IEN))); % ~= DOF' (Similar to DOF in CALFEM, but transposed)
LM = zeros(3*nen,nel); % ~= EDOF' (Similar to EDOF in CALFEM, but transposed)
for i = 1 : nel
    LM(:,i)=reshape(ID(:,IEN(:,i)),3*nen,1);
end

%% Material parameters:
E_Y = 210e9;
nu_P = 0.3;
denisty=785e3;
damp=0.1;
D=hooke_strain(5,E_Y,nu_P);

%% Gauss-Legendre quadrature points:
[ gp_x,w_x ] = getGP( deg.p );
[ gp_y,w_y ] = getGP( deg.q );
[ gp_z,w_z ] = getGP( deg.r );
NQUADx = size(gp_x,2);
NQUADy = size(gp_y,2);
NQUADz = size(gp_z,2);

%% Stiffness matrix and load vector computation

% Element loop
K = zeros(ndof); % Needs to be changed to sparse for large problems!!
M = zeros(ndof);
C = zeros(ndof);
F = zeros(ndof,1);
for e = 1 : nel
    % NURBS coordinates; convention consistent with Algorithm 7 in Hughes
    ni = INN(IEN(1,e),1);
    nj = INN(IEN(1,e),2);
    nk = INN(IEN(1,e),3);
    
    % Check if element has zero measure
    if (KV.Xi(ni+1) == KV.Xi(ni)) || (KV.Eta(nj+1) == KV.Eta(nj)) || (KV.Zeta(nk+1) == KV.Zeta(nk))
        continue
    end
    
    Ke = zeros(nen*3);
    Fe = zeros(nen*3,1);
    Me = zeros(nen*3);
    Ce = zeros(nen*3);jishu=0;
    for i = 1 : NQUADx % Loop trough Gauss points高斯积分
        for j = 1 : NQUADy
            for k = 1 : NQUADz
                % Gauss point
                GP.xi_tilde = gp_x(i);
                GP.eta_tilde = gp_y(j);
                GP.zeta_tilde = gp_z(k);
                jishu=jishu+1
                % Get Basis, derivatives, and det(J) for current gauss pt
                [ R,dR_dx,Jdet ] = Shape_function( GP,e,deg,B,KV,INN,IEN);
                
                % Combine quadrature weights with det(J)
                Jmod = abs(Jdet)*w_x(i)*w_y(j)*w_z(k);
                
                % Build Ke
                [ Ke_ ] = Build_K_Local( dR_dx,Jmod,D,nen );   %solid 
               % [ Ke_ ] = Build_K_Local_shell( dR_dx,Jmod,D,nen );   %超参shell 
                Ke = Ke + Ke_;
                
                % Build Me
                 Ni=eye(3);
                 for s=1:length(R)
                     shamtx(:,(s-1)*3+1:s*3)=R(s)*Ni;
                 end
                 Me_=denisty*shamtx'*shamtx*Jmod;
                 Ce_=damp*shamtx'*shamtx*Jmod;
                 Me = Me + Me_;
                 Ce = Me + Me_;
            end
        end
    end

    % Global Assembly
    idx = LM(:,e)';
    K(idx,idx) = K(idx,idx) + Ke;
    M(idx,idx) = M(idx,idx) + Me;
    C(idx,idx) = C(idx,idx) + Ce;
    F(idx) = F(idx)+Fe;
end

%{
% Apply load at two corner nodes
loc_z=0.1;
loc_x=3.0;
constNod = [];
for i = 1 : numel(B)
    if B{i}(3) == loc_z
        if (B{i}(2) == 0.0) || (B{i}(2) == 0.1)
            if B{i}(1) == loc_x
                constNod=[constNod i];
            end
        end
        
    end
end
% Apply load in vertical direction on identified nodes that are listed in
% constNod:
F(ID(3,constNod)) = -3e4;
%}
%% Boundary condiitons
% Find controlpoints with x = 0 to constrain
loc_x=0;
constNod = [];
for i = 1 : numel(B)
    if B{i}(1) == loc_x
        constNod=[constNod i];
    end
end
bc1=reshape(ID(:,constNod),numel(ID(:,constNod)),1);
bc=[bc1,zeros(length(bc1),1)];

loc_f=3;
constf = [];
for i = 1 : numel(B)
    if B{i}(1) == loc_f
        constf=[constf i];
    end
end
bf1=reshape(ID(:,constf),numel(ID(:,constf)),1);


%% 时域初始计算 Newmark法

%给定初始值

deltt=0.01;alf=0.25;delt=0.5;
c0=1/(alf*deltt^2);c1=delt/(alf*deltt);c2=1/(alf*deltt); c3=1/(2*alf)-1;
c4=delt/alf-1; c5=deltt/2*(delt/alf-2); c6=deltt*(1-delt);  c7=delt*deltt;

%形成有效的质量矩阵
Md=K+c0*M+c1*C;
U=chol(Md);

dt= zeros(ndof,1); vt = zeros(ndof,1); at = zeros(ndof,1); att = zeros(ndof,1); t=0; Q=F;
for n=1:1/deltt
   t=t+deltt;
   %claculate the payload 计算有效载荷
   Qe=Q+M*(c0*dt+c2*vt++c3*at)+C*(c1*dt+c4*vt+c5*at);
    %claculate the displacement  U*U'*dtt= Qe
    mid=U\Qe;
    dtt=U'\mid;
     %claculate the velocity and acceleration
    att=c0*(dtt-dt)-c2*vt-c3*att;
    for j=1:length(bf1)
    vt(j,1)=5000;   %常速度载荷
    end
    %赋给下次迭代 
    dt=dtt;
    at=att;
    
end

% Rewrite to cell structure. (For plotting)
u = cell(size(B));
comb = cell(size(B));
for i = 1 : size(ID,2)
    u{i} = [dt(ID(:,i)); 0];
    comb{i} = B{i} + u{i};
end




%% plot
plotNurbsSolidElementSimple( KV, B )
title('Initial geometry')
plotNurbsSolidElementSimple( KV, comb )
title('Displacements')