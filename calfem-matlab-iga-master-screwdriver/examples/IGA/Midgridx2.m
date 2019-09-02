function [INBz1,INNz1,Br1,Bz11,Bz] =Midgridx2(B,lb,ub,KV,deg) 
Bz=cell(size(B));

Bz=B;

 Aj=0;
 for i=(size(Bz,1)-6):size(Bz,1) 
     Aj=Aj+1;
%     for j=1:size(B,2)
%      if  B{i,j,1}(2)>=10e-10
%          B{i,j,1}(2)=2;
%      elseif B{i,j,1}(2)<=-10e-10
%          B{i,j,1}(2)=-2;
%      end
%     end 
    mid=cell2mat(B(i,:,1));
    nub=find(abs(abs(mid(3,:))-max(mid(3,:)))<1e-8);
    for k=1:numel(nub)
     if  Bz{i,nub(k),1}(2)>=10e-10    
         Bz{i,nub(k),1}(2)=ub(1);             %刀头处优化为凸边形（8边形） 一个设计变量改变6组控制点 （0.8-1.8）
     elseif Bz{i,nub(k),1}(2)<=-10e-10
         Bz{i,nub(k),1}(2)=-ub(1);
     end
    end
end

%% optimization the handle of screwdriver
% x=[0.9,1.10,1.15,0.95,1.15,1.1,1.25,0.85,1.35,1.1,0.7,0.75];
Ajs=2;
for i=(size(Bz,1)-12):(size(Bz,1)-7)  %刀杆（不连刀头）处优化  6组控制点，6个设计变量 （0.8-1.2）
     
    for j=1:size(Bz,2)
        
         Bz{i,j,1}(2:3)=ub(Ajs).* B{i,j,1}(2:3);
    end    
    Ajs=Ajs+1;
end
for i=12:14   %刀柄内凹处   3个设计变量，3组控制点  （0.5,1.5）
   
    for j=1:size(Bz,2) 
         Bz{i,j,1}(2:3)=ub(Ajs).* B{i,j,1}(2:3);    
    end    
      Ajs=Ajs+1;
end
[ X_ ] = getSubDivKVValues( KV.Xi, 1); 
[Bz,KV.Xi]=RefineKnotSolidXi( Bz,deg.p,KV.Xi,X_ );
% [ Z_ ] =getSubDivKVValues( KV.Zeta, 1);
% [Bz,KV.Zeta]=RefineKnotSolidZeta( Bz,deg.r,KV.Zeta,Z_ );
%{
for i=1:size(B,2)
    for j=1:size(B,3)
        Bz{2,i,j}(2:3)= ub(1)*B{2,i,j}(2:3);
        Bz{1,i,j}(2:3)= ub(1)*B{1,i,j}(2:3);
%         Bz{2,i,j}(2:3)=( Bz{1,i,j}(2:3)+ Bz{3,i,j}(2:3))/2 ;
         Bz{3,i,j}= B{3,i,j};  
         Bz{2,i,j}(4)= B{2,i,j}(4);
         Bz{1,i,j}(4)= B{1,i,j}(4);

         Bz{2,i,j}(1)=B{2,i,j}(1);
         Bz{1,i,j}(1)=B{1,i,j}(1);
    end
end
%}
Xi=KV.Xi; Eta=KV.Eta; Zeta=KV.Zeta;
p=deg.p ; p=deg.q; r=deg.r ;
%{
p = 1;
q = 1;
r = 2;

% Xi = [0,0.5,1];
Xi = [0,1];
Xi = [Xi(1) Xi Xi(end)];
Eta = [0,1];
Eta = [Eta(1) Eta Eta(end)];
   Zeta = [0,0,1,1,2,2,3,3,4,4];
%    Zeta = [0,1,2,3,4,5,6,7,8];
%  Zeta = [0,0.5,1];
Zeta = [Zeta(1) Zeta Zeta(end)];
k_ = length(Xi);
l_ = length(Eta);
m_ = length(Zeta);
in1=cell2mat(B(size(B,1),:,:)); in0=cell2mat(B01(size(B01,1),:,:));
out1=cell2mat(B(1,:,:)); out0=cell2mat(B01(1,:,:));
fun1=in1(1,:,:).^2+in1(3,:,:).^2;
fun2=in0(1,:,:).^2+in0(3,:,:).^2;
c1=find(fun1<=min(min(fun1)));c0=find(fun2<=min(min(fun2)));
m1=mod(c1(1),size(B,2));z1=floor(c1(1)/size(B,2));
m0=mod(c0(1),size(B01,2));z0=floor(c0(1)/size(B01,2));
if m1==0
    m1=size(B,2);
else 
    z1=z1+1;
end
if m0==0
    m0=size(B,2);
else 
    z0=z0+1;
end
inmaxX = min(abs(in1(1,m1,z1)),abs(min(in0(1,m0,z0))))-0.005;
inminX = -inmaxX;
oumaxX =  max(max(max(out1(1,:,:))),max(max(out0(1,:,:))))+0.005;
ouminX = min(min(min(out1(1,:,:))),min(min(out0(1,:,:))));
maxY = max(max(max(in1(2,:,:))),max(max(in0(2,:,:))));
minY = min(min(min(in1(2,:,:))),min(min(in0(2,:,:))));
if maxY<abs(oumaxX-ouminX)
    maxY=abs(oumaxX-ouminX);
end
inmaxZ = inmaxX;
inminZ = inminX;
oumaxZ =  max(max(max(out1(3,:,:))),max(max(out0(3,:,:))));
ouminZ = min(min(min(out1(3,:,:))),min(min(out0(3,:,:))));

 Bz = cell(2,2,9);
s = 1/sqrt(2);
% Bz = cell(3,3,3);

% Bz{1,1,1} = [ouminX minY ouminZ 1];
% Bz{1,2,1} = [(oumaxX+ouminX)/2 minY ouminZ 1];
% Bz{1,3,1} = [oumaxX minY ouminZ 1];
% Bz{1,1,2} = [ouminX (maxY+minY)/2 ouminZ 1];
% Bz{1,2,2} = [(oumaxX+ouminX)/2 (maxY+minY)/2 ouminZ 1];
% Bz{1,3,2} = [oumaxX (maxY+minY)/2 ouminZ 1];
% Bz{1,1,3} = [ouminX maxY ouminZ 1];
% Bz{1,2,3} = [(oumaxX+ouminX)/2 maxY ouminZ 1];
% Bz{1,3,3} = [oumaxX maxY ouminZ 1];
% 
% Bz{2,1,1} = [ouminX minY (oumaxZ+ouminZ)/2 1];
% Bz{2,2,1} = [(oumaxX+ouminX)/2 minY (oumaxZ+ouminZ)/2 1];
% Bz{2,3,1} = [oumaxX minY (oumaxZ+ouminZ)/2 1];
% Bz{2,1,2} = [ouminX (maxY+minY)/2 (oumaxZ+ouminZ)/2 1];
% Bz{2,2,2} = [(oumaxX+ouminX)/2 (maxY+minY)/2 (oumaxZ+ouminZ)/2 1];
% Bz{2,3,2} = [oumaxX (maxY+minY)/2 (oumaxZ+ouminZ)/2 1];
% Bz{2,1,3} = [ouminX maxY (oumaxZ+ouminZ)/2 1];
% Bz{2,2,3} = [(oumaxX+ouminX)/2 maxY (oumaxZ+ouminZ)/2 1];
% Bz{2,3,3} = [oumaxX maxY (oumaxZ+ouminZ)/2 1];
% 
% Bz{3,1,1} = [ouminX minY oumaxZ 1];
% Bz{3,2,1} = [(oumaxX+ouminX)/2 minY oumaxZ 1];
% Bz{3,3,1} = [oumaxX minY oumaxZ 1];
% Bz{3,1,2} = [ouminX (maxY+minY)/2 oumaxZ 1];
% Bz{3,2,2} = [(oumaxX+ouminX)/2 (maxY+minY)/2 oumaxZ 1];
% Bz{3,3,2} = [oumaxX (maxY+minY)/2 oumaxZ 1];
% Bz{3,1,3} = [ouminX maxY oumaxZ 1];
% Bz{3,2,3} = [(oumaxX+ouminX)/2 maxY oumaxZ 1];
% Bz{3,3,3} = [oumaxX maxY oumaxZ 1];
% 

Bz{1,1,1} = [oumaxX minY 0 1];
Bz{1,2,1} = [oumaxX maxY 0 1];
Bz{1,1,2} = [oumaxX minY oumaxZ s];
Bz{1,2,2} = [oumaxX maxY oumaxZ s];
Bz{1,1,3} = [0 minY oumaxZ 1];
Bz{1,2,3} = [0 maxY oumaxZ 1];
Bz{1,1,4} = [ouminX minY oumaxZ s];
Bz{1,2,4} = [ouminX maxY oumaxZ s];
Bz{1,1,5} = [ouminX minY 0 1];
Bz{1,2,5} = [ouminX maxY 0 1];
Bz{1,1,6} = [ouminX minY ouminZ s];
Bz{1,2,6} = [ouminX maxY ouminZ s];
Bz{1,1,7} = [0 minY ouminZ 1];
Bz{1,2,7} = [0 maxY ouminZ 1];
Bz{1,1,8} = [oumaxX minY ouminZ s];
Bz{1,2,8} = [oumaxX maxY ouminZ s];
Bz{1,1,9} = [oumaxX minY 0 1];
Bz{1,2,9} = [oumaxX maxY 0 1];

Bz{2,1,1} = [inmaxX minY 0 1];
Bz{2,2,1} = [inmaxX maxY 0 1];
Bz{2,1,2} = [inmaxX minY inmaxZ s];
Bz{2,2,2} = [inmaxX maxY inmaxZ s];
Bz{2,1,3} = [0 minY inmaxZ 1];
Bz{2,2,3} = [0 maxY inmaxZ 1];
Bz{2,1,4} = [inminX minY inmaxZ s];
Bz{2,2,4} = [inminX maxY inmaxZ s];
Bz{2,1,5} = [inminX minY 0 1];
Bz{2,2,5} = [inminX maxY 0 1];
Bz{2,1,6} = [inminX minY inminZ s];
Bz{2,2,6} = [inminX maxY inminZ s];
Bz{2,1,7} = [0 minY inminZ 1];
Bz{2,2,7} = [0 maxY inminZ 1];
Bz{2,1,8} = [inmaxX minY inminZ s];
Bz{2,2,8} = [inmaxX maxY inminZ s];
Bz{2,1,9} = [inmaxX minY 0 1];
Bz{2,2,9} = [inmaxX maxY 0 1];

% Number of weights
n = size(Bz,1);
m1 = size(Bz,2);
l = size(Bz,3);
sj=cell2mat(Bz);

% Order of basis
p = k_-n-1;
q = l_-m1-1;
r = m_ -l -1;
%}
% 
% % %插入网格的粗网格
% [ X_ ] = getSubDivKVValues( Xi, 1); 
% % [Bz,Xi]=RefineKnotSolidXi( Bz,p,Xi,X_ );
%  [Bz,Xi,Rx]=RefineKnotSolidXi( Bz,deg.p,Xi,X_ );
% [ E_ ] = getSubDivKVValues( Eta,1);
% [Bz,Eta]=RefineKnotSolidEta( Bz,q,Eta,E_ );
% [ Z_ ] =getSubDivKVValues( Zeta,1);
% [Bz,Zeta]=RefineKnotSolidZeta( Bz,r,Zeta,Z_ );
% KVz.Xi=Xi; KVz.Eta=Eta; KVz.Zeta=Zeta;
 degz=deg;Bz11=Bz;
 
 Bz1=cell2mat(Bz);Bzz1=Bz(:,2:size(Bz,3),:);
[INNz1,IENz1,INBz1] = BLDINCIEN( degz,size(Bz,1),size(Bz,2),size(Bz,3));
% figure (3)
% plotNurbsSolidElementSimple( KVz, Bz)

n = size(Bz,1);
m = size(Bz,2);
l = size(Bz,3);

 e=0;bi=0;
% for k = 1 : (l-1)
%         for j = 1 : (m-1)
%             for i = 1 : (n-1)
%                 e=e+1;
%                for kloc = 0 : 1   %1 ,  3
%                    for jloc = 0 : 1   %2 ,  12
%                            for iloc = 0 : 1    %3 ,  55                     
%                               if (kloc+k)==(l+1)
%                                
%                                 B_j=INBz1((iloc+i), jloc+j,1); 
%                                else
%                                 B_j=INBz1((iloc+i), jloc+j, kloc+k);   % global function number
%                                end
%                               
%                               b = (1+1)*(1+1)*(1+1)-kloc*(1+1)*(1+1) - jloc*(1+1)- iloc ; % local function number
%                               Br(b,e) = B_j; % assign connectivity
%                            end
%                    end
%                end
%              
%             end
%         end
% end

for k = 1 : (l-1)
        for j = 1 : (m-1)
            for i = 1 : (n-1)
                e=e+1;
               for kloc = 0 : 1   %1 ,  3
                   for jloc = 0 : 1   %2 ,  12
                           for iloc = 0 : 1    %3 ,  55                     
                              if (jloc+j)==(m)
                               
                                B_j=INBz1((iloc+i), 1, kloc+k); 
                               else
                                B_j=INBz1((iloc+i), jloc+j, kloc+k);   % global function number
                              end
                              
                              b = (1+1)*(1+1)*(1+1)-kloc*(1+1)*(1+1) - jloc*(1+1)- iloc ; % local function number
                              Br(b,e) = B_j; % assign connectivity
                           end
                   end
               end
             
            end
        end
end

    Br1=Br;
    Br1(1,:)=Br(2,:);
    Br1(2,:)=Br(1,:);
    Br1(5,:)=Br(6,:);
    Br1(6,:)=Br(5,:);
    
% %插入网格的细网格
% [ X_ ] = getSubDivKVValues( Xi, 2); 
% % [Bz,Xi]=RefineKnotSolidXi( Bz,p,Xi,X_ );
% [Bz,Xi,Rx]=RefineKnotSolidXi( Bz,deg.p,Xi,X_ );
% [ E_ ] = getSubDivKVValues( Eta,1);
% % [Bz,Eta]=RefineKnotSolidEta( Bz,q,Eta,E_ );
% [Bz,Eta,Ry]=RefineKnotSolidEta( Bz,deg.q,Eta,E_ );
 [ E_ ] = getSubDivKVValues( Eta,1);
 [Bz,Eta]=RefineKnotSolidEta( Bz,deg.q,Eta,E_ );
[ Z_ ] =getSubDivKVValues( Zeta, 2);
[Bz,Zeta,Rz]=RefineKnotSolidZeta( Bz,deg.r,Zeta,Z_ );
KVz.Xi=Xi; KVz.Eta=Eta; KVz.Zeta=Zeta;

% i=1:size(Rx,1);j=1:size(Ry,1);k=1:size(Rz,1);
% s=1:size(Rx,2);t=1:size(Ry,2);l=1:size(Rz,2);
% [ii,jj,kk]=meshgrid(i,j,k); [ss,tt,ll]=meshgrid(s,t,l);
% 
% cn=ii+(jj-ones(size(jj)))*size(Rx,1)+(kk-ones(size(jj)))*size(Ry,1)*size(Rx,1); 
% xn=ss+(tt-ones(size(tt)))*size(Rx,2)+(ll-ones(size(ll)))*size(Ry,2)*size(Rx,2);
% zz=kron(Rx(i,s,size(Rx,3)),Ry(j,t,size(Ry,3)));
% RR=kron(zz,Rz(:,:,size(Rz,3)));
%  RR=kron(RR,eye(3));  %疏密网格间的限制矩阵
%                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                     
%  Rs=sparse(RR);
%  P1=RR';          %疏密网格间的插值矩阵
% clear cn xn zz ii jj ss tt ll kk
% numel(Rs(Rs~=0))

Bz2=cell2mat(Bz);
Bzz=Bz(:,2:size(Bz,2),:);
nodez=size(Bzz,1)*size(Bzz,2)*size(Bzz,3);
% degz.p=p;degz.q=q;degz.r=r;
% [INNz,IENz,INBz] = BLDINCIEN11( degz,size(Bz,1),size(Bz,2),size(Bz,3));
% for i=1:size(Bz,1)*size(Bz,2)*size(Bz,3)
%   h=INNz(i,:)  ;
%   gcoord_mid(i,:)=Bz{h(1),h(2),h(3)}(1:3);
% 
% end
save middle Bz2 Bz
% figure (4)
% plotNurbsSolidElementSimple( KVz, Bz)