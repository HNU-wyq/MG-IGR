function    obj= REiga(x)



tic
%% Load data
addpath(genpath('C:\Users\Administrator\Desktop\IGA\calfem-matlab-iga-master'));

% count=1;
%  run ./../NURBS/Solid/data_solid_chilun
%  run ./../NURBS/solid/data_solid10
 
% for i=1:size(B,1)
%   for j=1:size(B,2)
%         B{i,j,3}(2)=B{i,j,3}(2)+0.5;
%   end
% end
   
% for i=1:size(B,1)
%  
%         B{i,2,5}(1:2)=B{i,2,5}(1:2)*1.5;
%    
% end
load initial

%{
 pro;

% Xi =[0,0,1,1];
% 
% Eta = [0,0,0,1,2,3,4,4,4,];
% Eta = [Eta(1) Eta Eta(end)];
% 
% Zeta = [0,0,0,1,1,2,2,3,3,4,4,4];
deg.p = p; deg.q=q; deg.r = r; clear p q r

%% h-refinement
% Perform 3 subdivisons in Xi, and 1 in  Eta & Zeta
[ X_ ] = getSubDivKVValues( Xi, 1); 
[B,Xi]=RefineKnotSolidXi( B,deg.p,Xi,X_ );
[ E_ ] = getSubDivKVValues( Eta,2);
[B,Eta]=RefineKnotSolidEta( B,deg.q,Eta,E_ );
[ Z_ ] =getSubDivKVValues( Zeta, 2);
[B,Zeta]=RefineKnotSolidZeta( B,deg.r,Zeta,Z_ );
%}

  loc_y1=-19; 
%   loc_x=30;
modNod=[];
for i = 1 : numel(B)
 
%      if ( B{i}(2) == loc_y && B{i}(1) == loc_x)
           if (  B{i}(2)-loc_y1<=1e-10)
            modNod=[modNod i];
    
           end
   
end

% x=BestSol.Position;
  center1=[30,loc_y1,0]; center2=[30,B{1,(size(B,3)-2)}(2),1];
% for i=1:size(modNod,2)
% 
%  B{modNod(i)}(1:3)=(x(1).*(B{modNod(i)}(1:3)-center1')+center1');
% end

for i=1:size(B,1)
    for j=1:size(B,3)
      B{i,size(B,2),j}(1:3)=(x(1).*(B{i,size(B,2),j}(1:3)-center1')+center1');
      B{i,size(B,2)-2,j}(1:3)=(x(2).*(B{i,size(B,2)-2,j}(1:3)-center2')+center2');
    end
end
 
  %%  
BJ=cell2mat(B);
 n1=size(B,1);m1=size(B,2);l1=size(B,3);
 [INN1,IEN1] = BldINCIEN( deg,n1,m1,l1 );

    B1=B(:,:,1:(l1-1));
% Number of control points after refinement
n = size(B1,1);
m = size(B1,2);
l = size(B1,3);

%% Build connectivity arrays
% Knot vectors in analysis

% KV.Xi =Xi; KV.Eta = Eta; KV.Zeta = Zeta; clear Xi Eta Zeta

% Build connectivity arrays
nel = (n-deg.p) * (m-deg.q) * (l-deg.r); % number of elements
nnp = n*m*l; % number of global basis functions

nen = (deg.p+1)*(deg.q+1)*(deg.r+1); % number of local basis functions
ndof = nnp*3; % number of global degrees of freedom
ldof = nen*3; % number of local degrees of freedom
% Build connectivity arrays
[INN,IEN] = BldINCIEN( deg,n,m,l ); % = INC (is unique for IGA), IEN(=ENOD', if there is a CALFEM counterpart it might be called ENOD, but transposed)
ID = reshape(1:max(max(IEN))*3,3,max(max(IEN))); % ~= DOF' (Similar to DOF in CALFEM, but transposed)
LM = zeros(3*nen,nel); % ~= EDOF' (Similar to EDOF in CALFEM, but transposed) 在系统中单元节点或者说是自由度编号。
for i = 1 : nel
    LM(:,i)=reshape(ID(:,IEN(:,i)),3*nen,1);
end

%% Material parameters:
E_Y = 210e4;
nu_P = 0.3;
denisty=785;
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
% K = zeros(ndof);
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
    Ce = zeros(nen*3);
    
    for i = 1 : NQUADx % Loop trough Gauss points
        for j = 1 : NQUADy
            for k = 1 : NQUADz
                % Gauss point
                GP.xi_tilde = gp_x(i);
                GP.eta_tilde = gp_y(j);
                GP.zeta_tilde = gp_z(k);
                
                % Get Basis, derivatives, and det(J) for current gauss pt
                [ R,dR_dx,Jdet ] = Shape_function( GP,e,deg,B1,KV,INN,IEN);
                
                % Combine quadrature weights with det(J)
                Jmod = abs(Jdet)*w_x(i)*w_y(j)*w_z(k);
                
                % Build Ke
                [ Ke_,kinmtx ] = Build_K_Local( dR_dx,Jmod,D,nen );
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
%     M(idx,idx) = M(idx,idx) + Me;
%     C(idx,idx) = C(idx,idx) + Ce;
%     F(idx) = F(idx)+Fe;
end

% Apply load at two corner nodes
%     loc_z=0;loc_x=0;   %data_10的
%     loc_y=0;loc_x=0;
%   loc_y=8;loc_x=3;  %chilun
  loc_y=2.5;loc_x=30;  %guan
constNod = [];
for i = 1 : numel(B1)
%        if B{i}(3) == loc_z
%            if (B{i}(2) == 0.495||B{i}(2) == 0.505)
          if B1{i}(2) == loc_y
             if B1{i}(1) == loc_x
                constNod=[constNod i];
             end
          end
        
%         end
end
% Apply load in vertical direction on identified nodes that are listed in
% constNod:
F(ID(3,constNod)) = 3e5;

%% Boundary condiitons
% Find controlpoints with x = 0 to constrain
%   loc_z=-1.2;
%   loc_y=-1.3;
   loc_y=20;
constNod = [];
for i = 1 : numel(B1)
    if B1{i}(2) == loc_y
        constNod=[constNod i];
    end
end
bc1=reshape(ID(:,constNod),numel(ID(:,constNod)),1);
bc=[bc1,zeros(length(bc1),1)];

nn=length(bc1);
for i=1:nn
    c=bc1(i);
    K(c,:)=0;
    K(:,c)=0;
    K(c,c)=1;
  
end
% [a1,dx,d1,d2,d,Fr,sd] = R_IFU(K,F,L,K0,disp0);
a1 = R_IC(K,F,disp0);

%% Solve system
% [a,K]=solveq(K,F,bc);

%%
ai=zeros(3*nen,nel);num_node=zeros(nnp,1);
for q1 = 1 : nel
     ai(:,q1)=a1(LM(:,q1)) ;
     
      ni = INN(IEN(1,q1),1); nj = INN(IEN(1,q1),2);nk = INN(IEN(1,q1),3);
    
    % Check if element has zero measure
    if (KV.Xi(ni+1) == KV.Xi(ni)) || (KV.Eta(nj+1) == KV.Eta(nj)) || (KV.Zeta(nk+1) == KV.Zeta(nk))
        continue
    end
    gp=1;np=1;
    for i = 1 : NQUADx % Loop trough Gauss points
        for j = 1 : NQUADy
            for k = 1 : NQUADz
                % Gauss point
                GP.xi_tilde = gp_x(i);
                GP.eta_tilde = gp_y(j);
                GP.zeta_tilde = gp_z(k);
                
                % Get Basis, derivatives, and det(J) for current gauss pt
                [ R,dR_dx,Jdet ] = Shape_function( GP,e,deg,B1,KV,INN,IEN);
                
                % Combine quadrature weights with det(J)
                Jmod = abs(Jdet)*w_x(i)*w_y(j)*w_z(k);
                [ Ke_, kinmtx ] = Build_K_Local( dR_dx,Jmod,D,nen );
                strain(:,gp,q1)=kinmtx*ai(:,q1);
                siga=D*strain(:,gp,q1);
                stress(gp,q1,1:3)=siga(1:3);
          % von Mises stress
                stress(gp,q1,4)  = sqrt(siga(1)^2+siga(2)^2-...
                                   siga(1)*siga(2)+3*siga(3)^2); 
                stress_n(IEN(gp,q1),:)=stress(gp,q1,:);
                num_node(IEN(gp,q1))= num_node(IEN(gp,q1))+1;
                gp=gp+1;
            end
        end
    end
    
end



for i=1:n1*m1*l1
  h=INN1(i,:)  ;
gcoord(i,:)=B{h(1),h(2),h(3)}(1:3);
end
for i=1:nnp
    num=num_node(i);
    stress_n(i,:)=stress_n(i,:)./num;
end


%%------------环形模型需增加的的应力操作---------
%{
 mise=zeros(n1*m1*l1,4);ip=size(B1,2)*size(B1,1);
 for i=1:ip
     mise((l1*(i-1)+1):(l1*i-1),:)=stress_n(((l1-1)*(i-1)+1):(l1-1)*i,:);
     mise(l1*i,:)=stress_n((l1-1)*(i-1)+1,:);
 end
%}
 mise=zeros(n1*m1*l1,4);ip=size(B1,2)*size(B1,1);
 
 mise(1:n1*m1*(l1-1),:)=stress_n;
 mise((n1*m1*(l1-1)+1):end,:)=stress_n(1:ip,:);

  
 objiga1 = max( mise(:,4))
%  [objiga2,id] = max (abs(a1))
%  obj=[objiga1 objiga2]';
obj=objiga1;

%{
A = 0; e=0;
    
    for k = 1 : l
        for j = 1 : m
            for i = 1 : n
                
                A = A + 1;
                
                if i >= (1+1) & j >= (1+1) & k >= (1+1)
                    e = e + 1;
                    for kloc = 0 : 1
                        for jloc = 0 : 1
                            for iloc = 0 : 1
                                B_ = A - kloc*n*m - jloc*n - iloc; % global function number
                                b = kloc*(1+1)*(1+1) + jloc*(1+1)+ iloc + 1; % local function number
                                BR(b,e) = B_; % assign connectivity
                            end
                        end
                    end
                end
            end
        end
    end
    BR1=BR;
    BR1(1,:)=BR(2,:);
    BR1(2,:)=BR(1,:);
    BR1(5,:)=BR(6,:);
    BR1(6,:)=BR(5,:);
  tecplot(BR1,gcoord,stress_n);
  
 
% Rewrite to cell structure. (For plotting)




u = cell(size(B));
comb = cell(size(B));
for i = 1 : size(ID,2)
    u{i} = [a1(ID(:,i)); 0];
    comb{i} = B{i} + u{i};
end
%    comb(21,:,:)=comb(1,:,:);
%    comb1=cell2mat(comb);
%% plot
plotNurbsSolidElementSimple( KV, B )
title('Initial geometry_the circular plate coupled to a doubly curved barrel shape')   %data_solid10
% plotNurbsSolidElementSimple( KV, comb )
% title('Displacements')
%}
toc

end
