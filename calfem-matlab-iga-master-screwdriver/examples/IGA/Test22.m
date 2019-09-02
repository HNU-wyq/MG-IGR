function [disp1,Kk1, Pp2] = Test22(S_kk1,S_ff1,P2,disp0,B1,INNz1,Br1,Bz11,Bz,nsub,ssub,idsub,Pp1)
%   profile on
% tt1 = cputime; 
%  [INBz1,INNz1,Br1,Bz11,Bz] =Midgrid(B,B01) ;

 % tt2 = cputime;
% dt(1) = tt2-tt1;
Pp2=transfer_matrix1(B1,Br1,INNz1,Bz11,Bz,nsub,ssub,idsub);
% [ Pp1, Pp2]=suoyin1(INBz1,INNz1,Br1,Bz11,B01,B1,Bz); %1初始,2参数化修改
clear Br1 Bz11 Bz INBz1 INNz1
% size(Pp1)
% tt3 = cputime;
% dt(2) = tt3-tt2;
%   S_kk1=sparse(S_kk1);
 Kk1 = Pp2'*S_kk1*Pp2;   %映射到虚拟网格M2,cond= 2.0099e+21

% s=0:size(Kk1,1)-1;
% b=s*(size(Kk1,1))+s+1;
% c= Kk1<1e-6;
%  tic
 [c1,c2]=find(Kk1);
 gf=c1((c1==c2));
 ij=setdiff([1:size(Kk1)]',gf);
 Kk1=Kk1+sparse(ij,ij,ones(size(ij)),size(Kk1,1),size(Kk1,2));
%  toc
%  logical_index=b(c(c~=0));
% logical_index=intersect(find(Kk1==0),b);
% tic
% % Kk1=sparse(Kk1);
% % Kk1([c1,c2])=1;
% % Kk1(logical_index)=1;
% % Kk1(sub2ind(size(Kk1), c1, c2))= 1;
% Kk1(b(c))= 1;
% toc
% tic
% sjj=b(c);
% X=ones(1,numel(sjj));
% kk1= sparse(sjj,sjj,X,size(Kk1,1),size(Kk1,2)) ;
% Kk1=Kk1+kk1;
% toc
Kk2=Pp1*Kk1*Pp1';

sdof1=size(Kk2,1);
% Kk2(abs(diag(Kk2))<1e-6)=1;
%  for i = 1:sdof1
%     if(abs(Kk2(i,i))<1e-6)
%         Kk2(i,i) = 1;
%     end
%  end
%  
 s=0:size(Kk2,1)-1;
b=s*(size(Kk2,1))+s+1;
c= Kk2(b)<1e-6;
logical_index=b(c(c~=0));
Kk2(logical_index)=1;
 Kk2=sparse(Kk2);
 Kk1=sparse(Kk1);
%  
r0 = zeros(sdof1,1);
sdof0=numel(disp0);
r0(1:sdof0) = disp0;
% size(r0)

ff= Pp2'*S_ff1;
% [disp0,B,B1] = PCG(1,kk,r0,ff);
% tt4 = cputime;
% dt(3) = tt4-tt3;
rr0 = AMG111(Kk1,ff,Kk2,Pp1',P2,r0);
% tt5 = cputime;
% dt(4) = tt5-tt4;

disp1 = AMG122(S_kk1,S_ff1,Kk1,Pp2,rr0);
end
% profile report
% profile close
% tt6 = cputime;
% dt(5) = tt6-tt5;