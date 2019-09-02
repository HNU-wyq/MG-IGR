function Pp=transfer_matrix1(B01,Br1,INNz1,Bz11,Bz,nsub,ssub,idsub)
% tim1=cputime;
len = norm(Bz11{Br1(1)}(1:3) -Bz11{Br1(2)}(1:3));%=0.055为表精确，避免极小的误差，直接用数字表示

nin = zeros(numel(B01),1);
in1 = zeros(numel(B01),20);
din1 = zeros(numel(B01),20);
p_detecte=zeros(4,numel(Bz));
p_detect=zeros(3,numel(Bz));
pt1=zeros(4,numel(Bz));
pt=zeros(3,numel(Bz));
      p_detecte=reshape(cell2mat(Bz),4,numel(cell2mat(Bz))/4);
      p_detect= p_detecte(1:3,:);
      pt1=reshape(cell2mat(Bz),4,numel(cell2mat(Bz))/4);
      pt= pt1(1:3,:);
for s=1:size(Br1,2)
%      [zf11,zf22]=dian_duobianxing_suoyin(i,s,Br1,INNz1,Bz11,B01);
%      e=0;
%     for k=1:size(Br1,1)
%     e=e+1; syin(k,:)=INNz1(Br1(k,s),:);
%     XYZ1(:,e)=Bz11{syin(k,1),syin(k,2),syin(k,3)};  
%     XYz(:,e)= XYZ1(1:3,e);      %一个单元的所有坐标
%     end
     XYZ1=[];
     k=(s-1)*size(Br1,1)+1:s*size(Br1,1);
     syin(1:numel(k),:)=INNz1(Br1(k),:);

%      XYZ1(:,1:numel(k))=Bz11{syin(1:numel(k),1),syin(1:numel(k),2),syin(1:numel(k),3)};  
%      XYz(:,1:numel(k))= XYZ1(1:3,1:numel(k));
     XYZ1=[XYZ1 cell2mat(Bz11(sub2ind(size(Bz11),syin(:,1),syin(:,2),syin(:,3))))];  
     XYz1= reshape(XYZ1,4,numel(XYZ1)/4);XYz=XYz1(1:3,:);
     
    p1=XYz(:,1);
    p4=XYz(:,4); 
    p5=XYz(:,5);
    p2=XYz(:,2); 
    p7=XYz(:,7);
    p6=XYz(:,6);
    p3=XYz(:,3);
    p8=XYz(:,8);
    a=p4-p1;  b=p5-p1;  c=p2-p1;
    d1=cross(c,b);  d2=cross(a,c);  d3=cross(b,a);    %di为各个面的外向法向量
    a2=p6-p7;  b2=p3-p7;  c2=p8-p7;
    d4=cross(b2,c2);  d5=cross(c2,a2);  d6=cross(a2,b2);
    
    
%       p_detecte=reshape(cell2mat(Bz),4,numel(cell2mat(Bz))/4);
%       p_detect= p_detecte(1:3,:);
%       pt1=reshape(cell2mat(Bz),4,numel(cell2mat(Bz))/4);
%       pt= pt1(1:3,:);
%         for  f=1:size(Br1,1)
            v1= bsxfun(@minus,p_detect,p1);
            v2= bsxfun(@minus,p_detect,p7);
%           v1=p_detect-XYz(:,1);
%           v2=p_detect- XYz(:,7);
         
          zf1(1,:)=(sum(bsxfun(@times,v1,d1)));
          zf1(2,:)=(sum(bsxfun(@times,v1,d2)));
          zf1(3,:)=(sum(bsxfun(@times,v1,d3)));
          zf2(1,:)=(sum(bsxfun(@times,v2,d4)));
          zf2(2,:)=(sum(bsxfun(@times,v2,d5)));
          zf2(3,:)=(sum(bsxfun(@times,v2,d6)));
%     for i = 1:numel(B01)
%          pt = B01{i}(1:3);
%            p_detected(1:4,i) =B01{i}; p_detect= p_detected(1:3,i);
%           v1=p_detect- XYz(:,1);
%           v2=p_detect- XYz(:,7);
%           zf1(1:3,i)=(sum(bsxfun(@times,v1,[d1, d2, d3] )));
%           zf2(1:3,i)=(sum(bsxfun(@times,v2,[d4, d5, d6] )));
%           num=intersect([find(sum(zf1<=0)==3)],[find(sum(zf2<=0)==3)]);
% %           zf1(1:3,i)=[ dot(v1,d1),dot(v1,d2),dot(v1,d3)];
% %           zf2(1:3,i)=[ dot(v2,d4),dot(v2,d5),dot(v2,d6)];
%           zf11(1:3)=zf1(1:3,i); 
%           zf22(1:3)=zf2(1:3,i);
%         if all(all(zf11<=0))&&all(all(zf22<=0))   %判断点是否在空间几何体里
          num=intersect([find(sum(zf1<=0)==3)],[find(sum(zf2<=0)==3)]);
          s1=syin(size(syin,1),:);  
%     if(x>=minx&&x<=maxx&&y>=miny&&y<=maxy)
    
          ix = s1(1);iy = s1(2);iz = s1(3);
             
     for j = max(ix,1):min(ix+1,size(Bz11,1)-1)
         for k = max(iy,1):min(iy+1,size(Bz11,2)-1);
          for m = max(iz,1):min(iz+1,size(Bz11,3)-1);
            l = idsub(j,k,m);
             for k1 = 1:nsub(j,k,m)
                mi=Bz{ssub(l,k1)}(1:3);
%                 dp=bsxfun(@minus,mi,pt);
              dp=repmat(mi,1,size(pt,2))-pt;
%                 lp = norm(dp);
                lp = sum(dp.^2,1).^0.5;
                len = norm(Bz11{Br1(1,s)}(1:3) -Bz11{Br1(2,s)}(1:3));
%                  if(lp(num)<len/2)   %%R%%%%?%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %                     nin(i) = nin(i)+1;
% %                     in1(i,nin(i)) = ssub(l,k1);
% %                     din1(i,nin(i)) = lp;
% %                     num1=find(lp(num)<len/2);
% %                     nin(s)=numel(num1);
% %                     in1(s,nin(num1)) = ssub(l,k1);
% %                     din1(s,nin(num1)) = lp(num1);
%                      nin(num)=nin(num)+ones(numel(num),1);
%                     in1(sub2ind(size(din1),num',nin(num))) = ssub(l,k1);
%                     din1(sub2ind(size(din1),num',nin(num))) = lp(num);
%                  end
                 if(lp<len/2)   %%R%%%%?%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                     nin(i) = nin(i)+1;
%                     in1(i,nin(i)) = ssub(l,k1);
%                     din1(i,nin(i)) = lp;
%                     num1=find(lp(num)<len/2);
%                     nin(s)=numel(num1);
%                     in1(s,nin(num1)) = ssub(l,k1);
%                     din1(s,nin(num1)) = lp(num1);
                     nin(num)=numel(num);
                    in1(num',nin(num))= ssub(l,k1);
                    din1(num',nin(num)) = lp(num);
                 end
            end
          end 
         end
     end
%       break;
%          end

%    end
end

nd = 30;
in = zeros(numel(B01),nd);
din = zeros(numel(B01),nd);

for i = 1:numel(B01)
    pt = in1(i,1:nin(i));
    r = din1(i,1:nin(i));
    [r1,s] = sort(r);
    pt1 = pt(s);
    np = min(nin(i),nd);
    nin(i) = np;
    in(i,1:nin(i)) = pt1(1:nin(i));  %suoyin
    din(i,1:nin(i)) = r1(1:nin(i));    %daxiao
end
node1 = numel(Bz);
j = numel(Bz);

for i = 1:numel(B01)
    if(nin(i)==0)
        nin(i) = 1;
        j = j+1;
        in(i,1) = j;
    end
end
node1 = j;
sdof1 = node1*3;
sdof=numel(B01)*3;
ndof=3;
Pp = sparse(sdof,sdof1);   %修改后网格与初始网格间的传递算子
for i = 1:numel(B01)
    l = din(i,1:nin(i));
    l2 = l.^2+1e-19;    %针对同一点距离为0的情况
    dl2 = 1./l2;
    sdl2 = sum(dl2);
    w = dl2/sdl2;
    for j = 1:nin(i)
        nd = in(i,j);
        Pp((i-1)*ndof+1:(i-1)*ndof+ndof,(nd-1)*ndof+1:(nd-1)*ndof+ndof) = w(j)*eye(ndof);
    end
end
% tim2=cputime;
% time12=tim2-tim1
end