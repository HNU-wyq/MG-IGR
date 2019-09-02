clc

clear all
clf
tic
%%% 层合壳元
addpath(genpath('C:\Users\Administrator\Desktop\calfem-matlab-iga-master - (2)'));
count=1;

%  data_l;
%      run ./../NURBS/Solid/data_solid_chilun
%           run ./../NURBS/solid/data_solid9
%           run ./../NURBS/Solid/data_part2
%             run ./../NURBS/Solid/data_multi_part
      
%  for i=1:size(B,1)
%   for j=1:size(B,2)
%         B{i,j,3}(2)=B{i,j,3}(2)+0.5;
%   end
% end   
%}

%  pro5;
 B=cell(3,3);
 B{1,1}=[-1,-1,0,1];
 B{1,2}=[0,-1,0,1];
 B{1,3}=[1,-1,0,1];
 B{2,1}=[-1,0,0,1];
 B{2,2}=[0,0,0,1];
  B{2,3}=[1,0,0,1];
   B{3,1}=[-1,1,0,1];
 B{3,2}=[0,1,0,1];
  B{3,3}=[1,1,0,1];
 Xi=[0,0,0,1,1,1];
 Eta=[0,0,0,1,1,1];
%   Eta=[0,0,1,1];
 Zeta=[0,0,1,1];
 p=numel(Xi)-size(B,1)-1;
  q=numel(Eta)-size(B,2)-1;
% load Pso BestSol
% for i=1:numel(B)
%     B{i}(1:2:3)=B{i}(1:2:3)./10;
% end
% B=permute(B,[3 2 1]);
deg.p = p; deg.q=q; 
% deg.r = r; 
clear p q 

%% h-refinement (and solve the transform matrix for mmultigrid) 

[ X_ ] = getSubDivKVValues( Xi, 3); 
[B,Xi,Rx]=RefineKnotSolidXi( B,deg.p,Xi,X_ );
[ E_ ] = getSubDivKVValues( Eta,3);
[B,Eta,Ry]=RefineKnotSolidEta( B,deg.q,Eta,E_ );
% [ Z_ ] =getSubDivKVValues( Zeta, 1);
% [B,Zeta,Rz]=RefineKnotSolidZeta( B,deg.r,Zeta,Z_ );
KV.Xi=Xi; KV.Eta=Eta;
% 
% i=1:size(Rx,1);j=1:size(Ry,1); k=1:size(Rz,1);
% s=1:size(Rx,2);t=1:size(Ry,2);l=1:size(Rz,2);
% [ii,jj,kk]=meshgrid(i,j,k); [ss,tt,ll]=meshgrid(s,t,l);
% 
% cn=ii+(jj-ones(size(jj)))*size(Rx,1)+(kk-ones(size(jj)))*size(Ry,1)*size(Rx,1); 
% xn=ss+(tt-ones(size(tt)))*size(Rx,2)+(ll-ones(size(ll)))*size(Ry,2)*size(Rx,2);
% zz=kron(Rx(i,s,size(Rx,3)),Ry(j,t,size(Ry,3)));
% RR=kron(zz,Rz(:,:,size(Rz,3)));
%  RR=kron(RR,eye(3));  %疏密网格间的限制矩阵
%  Rs=sparse(RR);
%  P1=RR';          %疏密网格间的插值矩阵
% clear RR cn xn zz ii jj kk ss tt ll
% numel(Rs(Rs~=0))

% % [ X_ ] = getSubDivKVValues( Xi, 1); 
% % [B,Xi,Rx]=RefineKnotSolidXi( B,deg.p,Xi,X_ );
% %  [ E_ ] = getSubDivKVValues( Eta,1);
% %  [B,Eta,Ry]=RefineKnotSolidEta( B,deg.q,Eta,E_ );
% [ Z_ ] =getSubDivKVValues( Zeta, 1);
% [B,Zeta,Rz]=RefineKnotSolidZeta( B,deg.r,Zeta,Z_ );
%  Rx=eye(size(Rx,2)); Ry=eye(size(Ry,2));
% KV.Xi=Xi; KV.Eta=Eta; KV.Zeta=Zeta;
% %对疏密网格自由度进行全局编号《得到总体自由度的转换矩阵的元素为
% i=1:size(Rx,1);j=1:size(Ry,1);k=1:size(Rz,1);
% s=1:size(Rx,2);t=1:size(Ry,2);l=1:size(Rz,2);
% [ii,jj,kk]=meshgrid(i,j,k);
% [ss,tt,ll]=meshgrid(s,t,l);
% % [rx,ry,rz]=meshgrid(Rx(:,:,size(Rx,3)),Ry(:,:,size(Ry,3)),Rz(:,:,size(Rz,3)));
% 
% cn=ii+(jj-ones(size(jj)))*size(Rx,1)+(kk-ones(size(jj)))*size(Ry,1)*size(Rx,1); 
% xn=ss+(tt-ones(size(tt)))*size(Rx,2)+(ll-ones(size(ll)))*size(Ry,2)*size(Rx,2);
% zz=kron(Rx(i,s,size(Rx,3)),Ry(j,t,size(Ry,3)));
% RR=kron(zz,Rz(:,2:size(Rz,2),size(Rz,3)));
%  RR=kron(RR,eye(3));  %疏密网格间的限制矩阵
%  Rs2=sparse(RR);
% P2=RR';          %疏密网格间的插值矩阵
% clear RR cn xn zz ii jj kk ss tt ll
% numel(Rs(Rs~=0))
% 
% %{
% % Perform 3 subdivisons in Xi, and 1 in  Eta & Zeta
% [ X_ ] = getSubDivKVValues( Xi, 1); 
% [B,Xi]=RefineKnotSolidXi( B,deg.p,Xi,X_ );
% [ E_ ] = getSubDivKVValues( Eta,2);
% [B,Eta]=RefineKnotSolidEta( B,deg.q,Eta,E_ );
% [ Z_ ] =getSubDivKVValues( Zeta, 2);
% [B,Zeta]=RefineKnotSolidZeta( B,deg.r,Zeta,Z_ );
% %}
% 
% KV.Xi=Xi; KV.Eta=Eta; KV.Zeta=Zeta;
% BJ=cell2mat(B);
% for i=1:size(B,1)
%     for j=1:size(B,2)
%      Norm(i,j)=   norm(B{i,j,1}(2:3));
%     end
% end
% Max=round(max(Norm(1,:)));
% if Max<max(Norm(1,:))
%  Max=Max+1;
% end
% Min=round(min(Norm(3,:)));
% if Min>min(Norm(3,:))
%     Min=Min-1;
% end
%  Xmax=max(max(max(BJ(1:4:end,:,:))))+1;
%  Xmin=min(min(min(BJ(1:4:end,:,:))))-1;
%  Ymax=max(max(max(BJ(2:4:end,:,:))))+1;
%  Ymin=min(min(min(BJ(2:4:end,:,:))));
%  Zmax=max(max(max(BJ(3:4:end,:,:))))+1;
%  Zmin=min(min(min(BJ(3:4:end,:,:))));
   n1=size(B,1);m1=size(B,2);l1=size(B,3);
  [INN1,IEN1] = BldINCIEN( deg,n1,m1,l1 );
 B1=B(:,:,1:size(B,3));  %(9*3*12)  %无重合点时


% nVar=180;
% parent_Position1=1.2*ones(1,180);
% for i=1:nVar
%     B1{1+3*(i-1)}(2:3)=parent_Position1(i).*B1{1+3*(i-1)}(2:3);
%     B1{2+3*(i-1)}(2:3)=parent_Position1(i).*B1{2+3*(i-1)}(2:3);
%     B1{3*(i)}(2:3)=parent_Position1(i).*B1{3*(i)}(2:3);
% end

% Number of control points after refinement
n = size(B1,1);
m = size(B1,2);
l = size(B1,3);
%   B1=B(:,:,1:l);
%% Build connectivity arrays
% Knot vectors in analysis
KV.Xi =Xi; KV.Eta = Eta; 
% KV.Zeta = Zeta; 
% clear Xi Eta

% Build connectivity arrays
% nel = (n-deg.p) * (m-deg.q) * (l1-deg.r); % number of elements
nel = (n-deg.p) * (m-deg.q) ; % number of elements

nnp = n*m; % number of global basis functions

nen = (deg.p+1)*(deg.q+1); % number of local basis functions
sdof = nnp*6; % number of global degrees of freedom
ldof = nen*6; % number of local degrees of freedom
% Build connectivity arrays
% [INN,IEN,AA] = BldINCIEN( deg,n,m,l ); % = INC (is unique for IGA), IEN(=ENOD', if there is a CALFEM counterpart it might be called ENOD, but transposed)
%   [INN,IEN,INB] = BLDINCIEN( deg,n,m,l );
   [INN,IEN,INB] =BldINCIEN( deg,n,m,l );

ID = reshape(1:max(max(IEN))*6,6,max(max(IEN))); % ~= DOF' (Similar to DOF in CALFEM, but transposed)
LM = zeros(6*nen,nel); % ~= EDOF' (Similar to EDOF in CALFEM, but transposed) 在系统中单元节点或者说是自由度编号。
for i = 1 : nel
    LM(:,i)=reshape(ID(:,IEN(:,i)),6*nen,1);
end
ndof=size(IEN,1);  %单元节点数
%% Material parameters:
thickness=0.1;   %plate thickness
density=1000;
nc=4 ;        %层数number of layers
El = 181e9;       % 增大
Et = 10.270e9;        % 增大
vlt = 0.28;          % 减小
vtl = Et*vlt/El;
% vtl = 0.14;
Glt = 7.17e9;       % 增大
Glw = 7.17e9;
Gtw = 3.780e9;        % 增大
E=[El Et] ;   %沿纤维方向弹模1和横向弹模2
v=[vlt vtl];   %主泊松比vlt，
Q = [El/(1-vlt*vtl) vtl*El/(1-vlt*vtl) 0;
    vtl*El/(1-vlt*vtl) Et/(1-vlt*vtl)  0;
    0     0         Glt];
% %% Material parameters:
% E_Y = 210e4;
% nu_P = 0.3;
% denisty=785;
% damp=0.1;
% %D=hooke_strain(5,E_Y,nu_P);

%% Gauss-Legendre quadrature points: 
[ gp_x,w_x ] = getGP( deg.p );     
[ gp_y,w_y ] = getGP( deg.q );
% [ gp_z,w_z ] = getGP( deg.r );
NQUADx = size(gp_x,2);
NQUADy = size(gp_y,2);
% NQUADz = size(gp_z,2);

%% Stiffness matrix and load vector computation
Bmatrix=cell(size(IEN));
% Element loop
K = zeros(sdof); % Needs to be changed to sparse for large problems!!
 M = zeros(sdof);
F = zeros(sdof,1);
h=0.2;
% z1=[-h/4 ;   0; h/4; h/2];  %upper
% z2=[-h/2 ;-h/4;   0; h/4];   %lower
% Zt=[z1 z2];

z= 0:h/nc:h;
zm = h/2;
zi = z-zm;
Zt = zeros(nc,2);
Zt(:,1)=zi(1:nc);
Zt(:,2)=zi(2:nc+1);

thetai=zeros(1,4);
for e =1:nel
    je=1;
    % NURBS coordinates; convention consistent with Algorithm 7 in Hughes
    ni = INN(IEN(1,e),1);
    nj = INN(IEN(1,e),2);
%     nk = INN(IEN(1,e),3);

    % Check if element has zero measure
    if ((KV.Xi(ni+1) == KV.Xi(ni)) || (KV.Eta(nj+1) == KV.Eta(nj)))
        continue
    end
   
    Ke = zeros(nen*6);
    Fe = zeros(nen*6,1);
     Me = zeros(nen*6);
%     Ce = zeros(nen*3);

%      for k = 1 : NQUADz
%      for i = 1 : NQUADx % Loop trough Gauss points
         for j = 1 : NQUADy
             for i = 1 : NQUADx 
%               for k = 1 : NQUADz
                % Gauss point
                GP.xi_tilde = gp_x(i);
                GP.eta_tilde = gp_y(j);
%                 GP.zeta_tilde = gp_z(k);
                
                % Get Basis, derivatives, and det(J) for current gauss pt
                [ R,dR_dx,Jdet ] =D2_Shape_function( GP,e,deg,B,KV,INN,IEN);
                % Combine quadrature weights with det(J)
                Jmod = abs(Jdet)*w_x(i)*w_y(j);
                
                % Build  Ke from membrane , dending ,shear components  
                [ Ke_ ,B_m,dme] = shell_Build_K_Local( dR_dx,Jmod,R ,ndof,Zt,Q,nc,Gtw,Glw,h,density);
%               [ Kg_ ,kinmtx ] = shell_Build_Kg_Local( dR_dx,Jmod,nen,thetai,R ,ndof);

                ns=[] ;
                for ii=1:ndof         % 一个节点所在多个单元处于同一平面内时，刚度矩阵会奇异，这里给第六自由度(对角)赋值消除奇异
                 ns=[ns, 1+6*(ii-1):5+6*(ii-1)];
                  Ke(6*ii,6*ii) = 0.05;  % 圆柱壳在0.05附近,越小频率越小,小于0.01不稳定；板对大小不敏感
                end
                 Ke(ns,ns)= Ke(ns,ns)+Ke_;
                 Me(ns,ns)= Me(ns,ns)+dme;
                 Bm{je,e}=B_m; %membrane B matrix
                 DR_dx{je,e}=dR_dx;
%                 
%                 a = 1:6:6*(ndof);
% aax = bsxfun(@plus,a,(0:4)')';
% bbx = bsxfun(@plus,a,(0:4)');
% aax = aax(:)';
% bbx = bbx(:)';
% Ke(aax,aax) = Ke(aax,aax)+Ke_;
% for ii=6:6:6*(ndof)         % 一个节点所在多个单元处于同一平面内时，刚度矩阵会奇异，这里给第六自由度(对角)赋值消除奇异
%   Ke(ii,ii) = 0.05;  % 圆柱壳在0.05附近,越小频率越小,小于0.01不稳定；板对大小不敏感
% end

% te=zeros(6);
% te(1:3,4:6)=eye(3);
% te(4:6,1:3)=eye(3);
% % te(6,6)=1;
% Te=blkdiag(te,te,te,te);
% ke=Te*ke*Te';
% me=Te*me*Te';
% te=zeros(6);
% te(1,4)=1;
% te(2,5)=1;
% te(3,1)=1;
% te(4,3)=1;
% te(5,2)=1;
% te(6,6)=1;
% Te=blkdiag(te,te,te,te,te,te);
% Ke=Te*Ke*Te';
% me(bbx,bbx) = me(bbx,bbx)+me1;
        je=je+1;
            end
        end
%     end

    % Global Assembly
    idx = LM(:,e)';
    K(idx,idx) = K(idx,idx) + Ke;
%     KG(idx,idx) = KG(idx,idx) + Kg;
     M(idx,idx) = M(idx,idx) + Me;
%     C(idx,idx) = C(idx,idx) + Ce;
%     F(idx) = F(idx)+Fe;
end

%% Apply load in vertical direction on identified nodes that are listed in
% constNod:
constNod=[];
constNod=[constNod reshape(INB(2:size(INB,1),1,:),1,numel(INB(2:size(INB,1),1,:)))];
F(ID(1,constNod)) = 8e8;
constNod=[];
constNod=[constNod reshape(INB((2:size(INB,1)),size(INB,2),:),1,numel(INB((2:size(INB,1)),size(INB,2),:)))];
F(ID(1,constNod)) = -8e8;

%% Boundary condiitons
% Find controlpoints with x = 0 to constrain
constNod = [];
constNod=[constNod reshape(INB(:,size(INB,2),:),1,numel(INB(:,size(INB,2),:)))];
constNod=[constNod reshape(INB(:,1,:),1,numel(INB(:,1,:)))];
constNod=[constNod reshape(INB(1,:,:),1,numel(INB(1,:,:)))];
constNod=[constNod reshape(INB(size(INB,1),:,:),1,numel(INB(size(INB,1),:,:)))];
constNod=unique(constNod);
bc1=reshape(ID(3,constNod),numel(ID(3,constNod)),1);
constNod = [];
constNod=[constNod INB(1,1,:)];
bc2=reshape(ID(2,constNod),numel(ID(2,constNod)),1);
constNod = [];
constNod=[constNod INB(1,size(INB,2))];
bc3=reshape(ID(1:2,constNod),numel(ID(1:2,constNod)),1);
bc1=unique([bc1;bc2;bc3]);
bc=[bc1,zeros(length(bc1),1)];
% n=length(bc1);
%  K=sparse(K);


%% Solve system static
[a,K]=solveq(K,F,bc);
% ex4(Xi,Eta,B,deg.p,deg.q )   

%% calculate the contribuation of in-plane stress and construct the geometric matrix KG

ai=zeros(5*nen,nel);num_node=zeros(nnp,1);
stress_n=zeros(nnp,3);
 
 KG=zeros(sdof);
for sq = 1:nel
     ai(:,sq)=a(LM(ns,sq)) ;
    Gb=zeros(2,5*ndof);
 Gs1=zeros(2,5*ndof);
 Gs2=zeros(2,5*ndof);
 KGb=zeros(5*ndof);
 KGs=zeros(5*ndof);
      ni = INN(IEN(1,sq),1); nj = INN(IEN(1,sq),2);
%       nk = INN(IEN(1,sq),3);
    
    % Check if element has zero measure
    if (KV.Xi(ni+1) == KV.Xi(ni)) || (KV.Eta(nj+1) == KV.Eta(nj)) 
        continue
    end
     stress_g = [];
     for j = 1 : NQUADy
         for i = 1 : NQUADx % Loop trough Gauss points
        
          
                % Gauss point
                GP.xi_tilde = gp_x(i);
                GP.eta_tilde = gp_y(j);
                GP.zeta_tilde = gp_z(k);
                
                % Get Basis, derivatives, and det(J) for current gauss pt
                [ R,dR_dx,Jdet ] = Shape_function( GP,q,deg,B,KV,INN,IEN);
                
                % Combine quadrature weights with det(J)
                Jmod = abs(Jdet)*w_x(i)*w_y(j)*w_z(k);
                
                % Build Ke
                [ Ke_,kinmtx ] = Build_K_Local( dR_dx,Jmod,D,nen );
                
                % Stress
                strain(:,gp,q)=kinmtx*ai(:,q);
                siga=D*strain(:,gp,q);
                gstress=zeros(1,4);
                gstress(1:3)=siga(1:3)';
                % von Mises stress
                gstress(4) = (((siga(1)-siga(2))^2+(siga(2)-siga(3))^2+(siga(3)-siga(1))^2)/2)^0.5;
                
                stress_g = [stress_g;gstress];
                num_node(IEN(gp,q))= num_node(IEN(gp,q))+1;
                gp=gp+1;
            end
        end
    end
    
    % stress recovery
    gpx = [1/gp_x(1) gp_x(2) 1/gp_x(3)];
    gpy = [1/gp_y(1) 1/gp_y(2)];
    gpz = [1/gp_z(1) gp_z(2) 1/gp_z(3)];
    H = [];
    for gpi = 1:NQUADx
        for gpj = 1:NQUADy
            for gpk = 1:NQUADz
                GP1.xi_tilde = gpx(gpi);
                GP1.eta_tilde = gpy(gpj);
                GP1.zeta_tilde = gpz(gpk);
                H1 = Shape_function( GP1,q,deg,B,KV,INN,IEN)';
                H = [H;H1];
            end
        end
    end
    stress = H*stress_g;
    
    
    
    gp=1;np=1;
    Kg=zeros(6*ndof);
%     for k = 1 : NQUADz
%         for i = 1 : NQUADx % Loop trough Gauss points
         for j = 1 : NQUADy
             for i = 1 : NQUADx 
                
                kinmtx=Bm{gp,sq};
                strain(:,gp,sq)=kinmtx*ai(:,sq);
                siga=zeros(3,1);
                dN=DR_dx{gp,sq};
                
                
                Gb(1,1:5:5*ndof)=dN(:,1)';
                Gb(2,1:5:5*ndof)=dN(:,2)';
                Gs1(1,2:5:5*ndof)=dN(:,1)';
                Gs1(2,2:5:5*ndof)=dN(:,2)';
                Gs2(1,3:5:5*ndof)=dN(:,1)';
                Gs2(2,3:5:5*ndof)=dN(:,2)';

                for il = 1:nc    %layers
                  thi = thetai(il);   %angles
                   Ti=[cos(thi)^2 sin(thi)^2 -2*cos(thi)*sin(thi);
                       sin(thi)^2 cos(thi)^2 2*cos(thi)*sin(thi);
                       cos(thi)*sin(thi) -cos(thi)*sin(thi) cos(thi)^2-sin(thi)^2];
                         z1 = Zt(il,1);   %Zt
                         z2 = Zt(il,2);
                          Qb=Ti*Q*Ti';
                    siga=siga+Qb*strain(:,gp,sq);
                    Sigma=[siga(1),siga(3); siga(3),siga(2)];
                    KGb=KGb+Gb'*Sigma*Gb*det(Jmod)*(z2-z1);
                    KGs=KGs+(Gs1'*Sigma*Gs1+Gs2'*Sigma*Gs2)*det(Jmod)*(z2-z1)^3/12;
                    
                end
               ns=[] ;
                for ii=1:ndof         % 一个节点所在多个单元处于同一平面内时，刚度矩阵会奇异，这里给第六自由度(对角)赋值消除奇异
                 ns=[ns, 1+6*(ii-1):5+6*(ii-1)];
%                  Kg(6*ii,6*ii) = 0.05;  % 圆柱壳在0.05附近,越小频率越小,小于0.01不稳定；板对大小不敏感
                end
                Kg(ns,ns)= Kg(ns,ns)+KGb+KGs;

                stress_n(IEN(gp,sq),:)=stress_n(IEN(gp,sq),:)+siga';
%                  num_node(IEN(gp,sq))= num_node(IEN(gp,sq))+1;
                gp=gp+1;
            end
        end
%     end   
                  idx = LM(:,sq)';
                  KG(idx,idx) = KG(idx,idx) + Kg;
end

for ii=1:nel
    for jj=1:size(IEN,1)
        num_node(IEN(jj,ii))= num_node(IEN(jj,ii))+1;
    end
 end
for i=1:n*m
  h=INN(i,:)  ;
gcoord(i,:)=B1{h(1),h(2)}(1:3);
end
for i=1:nnp
    num=num_node(i);
    stress_n(i,:)=stress_n(i,:)./num;
end

[V, D]=eig(K,KG);
 [D,Nu]=sort(diag(D));
% D=sort(D);
Vz=V(3:6:sdof,Nu);


%% TRANSFORM the matrix of element(p+1)*(q+1)(r+1) into 2*2*2
 e=0;bi=0;
    p= deg.p ; q=deg.q;
%     r=deg.r ;
% for k = 1 : (l)
        for j = 1 : (m-1)
            for i = 1 : (n-1)
                e=e+1;
%                for kloc = 0 : 1   %1 ,  3
                   for jloc = 0 : 1   %2 ,  12
                           for iloc = 0 : 1    %3 ,  55                     
%                               if (kloc+k)==(l+1)
%                                
%                                 B_j=INB((iloc+i), jloc+j,1); 
%                                else
                                B_j=INB((iloc+i), jloc+j);   % global function number
%                               end
                              
                              b = (1+1)*(1+1)- jloc*(1+1)- iloc ; % local function number
                              BR(b,e) = B_j; % assign connectivity
                           end
                   end
%                end
             
            end
        end
% end
    BR1=BR;
    BR1(1,:)=BR(2,:);
    BR1(2,:)=BR(1,:);
%     BR1(5,:)=BR(6,:);
%     BR1(6,:)=BR(5,:);
    
%   tecplot(BR1,gcoord,mise);
     tecplot2(BR1,gcoord,stress_n,a,Vz);%     tecplot(BR1,gcoord,stress_n);
%}
%}

%%
u = cell(size(B1)); 
comb = cell(size(B1));
for i = 1 : size(ID,2) %
u{i} = [a(ID(1:3,i)); 0]; 
comb{i} = B1{i} + u{i}; 
end
 ex4(Xi,Eta,comb,deg.p,deg.q )   
%        comb(:,:,l1)=comb(:,:,1);
%        comb1=cell2mat(comb);
%      tecplot(BR1,gcoord,stress_n,a);%     tecplot(BR1,gcoord,stress_n);
%% plot

plotNurbsSolidElementSimple( KV, B )
title('Initial geometry_the circular plate coupled to a doubly curved barrel shape')   %data_solid10

 plotNurbsSolidElementSimple( KV, comb )
 title('Displacements')
% 
 disp0=a;K0=K; ndof0=sdof;B0=B;

  save initial  disp0 ndof0  K0  B0 deg KV D B1 P1 P2 Max Min D LM  ID  INN IEN INB F
toc

