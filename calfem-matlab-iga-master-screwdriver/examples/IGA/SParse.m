
function b=SParse(b,bc1)

[u,v,s]=find(b);
% bc1=[2,3]';

%     c=bc1(i);
%     K(c,:)=0;
%     K(:,c)=0;
%     K(c,c)=1;
for i=1:length(bc1)
c=intersect(u,bc1);
[c1,c2]=find(u==c(i));
b(u(c1),v(c1))=0; 
end
for i=1:length(bc1)
d=intersect(v,bc1);
[d1,d2]=find(v==d(i));
b(u(d1),v(d1))=0;
end
for i=1:size(bc1)
b(bc1(i),bc1(i))=1;
end

end

  
