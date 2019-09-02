function [ Pp1, Pp2]=suoyin1(INBz1,INNz1,Br1,Bz11,B01,B1,Bz)

% ttt1 = cputime;
% node0 = numel(Bz);

idsub=zeros(size(INBz1,1)-1,size(INBz1,2)-1,size(INBz1,3)-1);
A =0;
for k = 1 : size(idsub,1)
        for j = 1 : size(idsub,2)
            for i = 1 : size(idsub,3)
                  A = A + 1;
                 idsub(i, j, k)=A;
            end
        end
end
% idsub=INBz1;
nsub=zeros(size(INBz1,1)-1,size(INBz1,2)-1,size(INBz1,3)-1);
% ssub = zeros(size(INBz1,1)*size(INBz1,2)*size(INBz1,3),20);
ssub = zeros(size(Br1,2),20);
% p_detected=zeros(4,numel(Bz));
% p_detect=zeros(3,numel(Bz));
 [nsub,ssub]=assmeble_grid(Br1,INNz1,Bz11,Bz,ssub,nsub,idsub);
% for i=1:size(Bz,1)*size(Bz,2)*size(Bz,3)
%  for j=1:size(Br1,2)
% %      [zf11,zf22]=dian_duobianxing_suoyin(i,j,Br1,INNz1,Bz11,Bz);
% %      e=0;
% %     for k=1:size(Br1,1)
% %      e=e+1; syin(k,:)=INNz1(Br1(k,j),:);
% %      XYZ1(:,e)=Bz11{syin(k,1),syin(k,2),syin(k,3)};  
% %      XYz(:,e)= XYZ1(1:3,e);
% %     end
%      XYZ1=[];
%      k=(j-1)*size(Br1,1)+1:j*size(Br1,1);
%      syin(1:numel(k),:)=INNz1(Br1(k),:);
%      XYZ1=[XYZ1 cell2mat(Bz11(sub2ind(size(Bz11),syin(:,1),syin(:,2),syin(:,3))))];  
%      XYz1= reshape(XYZ1,4,numel(XYZ1)/4);XYz=XYz1(1:3,:);
%     
%     p1=XYz(:,1);  
%     p4=XYz(:,4); 
%     p5=XYz(:,5);
%     p2=XYz(:,2); 
%     p7=XYz(:,7);
%     p6=XYz(:,6); 
%     p3=XYz(:,3);
%     p8=XYz(:,8);
%     a=p4-p1;  b=p5-p1;  c=p2-p1;
%     d1=cross(c,b); 
%     d2=cross(a,c);  
%     d3=cross(b,a); %di为各个面的外向法向量
%     a2=p6-p7;  b2=p3-p7;  c2=p8-p7;
%     d4=cross(b2,c2); 
%     d5=cross(c2,a2); 
%     d6=cross(a2,b2);
% 
% %      for i=1:size(Bz,1)*size(Bz,2)*size(Bz,3)
% %          p_detected(1:4,i) =Bz{i};
%            p_detected=reshape(cell2mat(Bz),4,numel(cell2mat(Bz))/4);
%            p_detect= p_detected(1:3,:);
% %         for  f=1:size(Br1,1)
%             v1= bsxfun(@minus,p_detect,p1);
%             v2= bsxfun(@minus,p_detect,p7);
% %           v1=p_detect-XYz(:,1);
% %           v2=p_detect- XYz(:,7);
%           zf1(1,:)=(sum(bsxfun(@times,v1,d1)));
%           zf1(2,:)=(sum(bsxfun(@times,v1,d2)));
%           zf1(3,:)=(sum(bsxfun(@times,v1,d3)));
%           zf2(1,:)=(sum(bsxfun(@times,v2,d4)));
%           zf2(2,:)=(sum(bsxfun(@times,v2,d5)));
%           zf2(3,:)=(sum(bsxfun(@times,v2,d6)));
%           num=intersect([find(sum(zf1<=0)==3)],[find(sum(zf2<=0)==3)]);
%           s1=syin(size(syin,1),:);
%           nsub(s1(1),s1(2),s1(3)) =numel(num);
%           ssub(idsub(s1(1),s1(2),s1(3)),1:nsub(s1(1),s1(2),s1(3)) ) = num;
% %           zf1(1:3,i)=(sum(bsxfun(@times,v1,[d1, d2, d3] )));
% %           zf2(1:3,i)=(sum(bsxfun(@times,v2,[d4, d5, d6] )));
%           
% %           zf1(1:3,i)=[ dot(v1,d1),dot(v1,d2),dot(v1,d3)];
% %           zf2(1:3,i)=[ dot(v2,d4),dot(v2,d5),dot(v2,d6)];
% %           zf11(1:3)=zf1(1:3,i); 
% %           zf22(1:3)=zf2(1:3,i);
%              
% 
% %         if all(all(zf11<=0))&&all(all(zf22<=0))    %判断点是否在空间几何体里
% %           s1=syin(size(syin,1),:);
% %           nsub(s1(1),s1(2),s1(3)) = nsub(s1(1),s1(2),s1(3))+1;
% %           ssub(idsub(s1(1),s1(2),s1(3)),nsub(s1(1),s1(2),s1(3)) ) = i; %中间网格粗细联立
% % %           break;
% %         end
% %      end
%   end


%}

 Pp1=transfer_matrix1(B01,Br1,INNz1,Bz11,Bz,nsub,ssub,idsub);%初始
%  size( Pp1)
 Pp2=transfer_matrix1(B1,Br1,INNz1,Bz11,Bz,nsub,ssub,idsub);%参数化修改

%  size( Pp2)

%  sum(sum(Pp1~=0,2))
%  sum(sum(Pp2~=0,2))
 
%  ttt2 = cputime;
%  sytime=ttt2-ttt1;
 
% P21=full(Pp2*Pp1');   %初始-》修改后
% sum(sum(P21~=0,2))


