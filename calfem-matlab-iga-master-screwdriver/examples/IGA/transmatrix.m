function Rm=transmatrix(Xii,Xi,p)

s1=length(Xii)-1-p;
s2=length(Xi)-1-p;
for j=1:s1  %xi
    for i=1:s2 %cu
         if (Xi(i)<=Xii(j)&&Xii(j)<Xi(i+1))
          Rr(i,j,1)=1;
         else 
          Rr(i,j,1)=0;
         end                
    end        
end
 Rr(s2+1,:,1)=0;
size(Rr);

for j=1:s1  %xi
    for i=1:s2 %cu        
        for rr=1:p
        alfa= Judge(rr,Xi,Xii,i,j);
        alfa1=Judge(rr,Xi,Xii,i+1,j);
        Rr(i,j,rr+1)= alfa* Rr(i,j,rr) + (1-alfa1)* Rr(i+1,j,rr);    %%%%%%()
     
        end
    end        
end
%  Rr(s2,:,rr+1)=zeros(1,j,1);
 Rm= Rr(1:s2,:,2:rr+1);
size(Rm);
end

 