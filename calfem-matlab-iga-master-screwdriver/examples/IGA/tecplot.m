function  tecplot(IEN,gcoord,stress_n,a)

nel = size(IEN,2);
node= size(gcoord,1);
nodedata=fopen('Mise.plt','w');               
fprintf(nodedata,'TITLE="data"\n');

fprintf(nodedata,'VARIABLES=,"X", "Y","Z" ,"sig4","ux","uy","uz"\n');

fprintf(nodedata,'ZONE T="%d  "  ,  N=%d, E=%d, ET=BRICK, F=FEPOINT\n',1,node,nel); % 注意这里的双引号
sigm1=stress_n(:,1);sigm2=stress_n(:,2);sigm3=stress_n(:,3);sigm4=stress_n(:,4);
ux=a(1:3:size(a));uy=a(2:3:size(a));uz=a(3:3:size(a));
for i=1:node
    fprintf(nodedata,'%10.10f,%10.10f,%10.10f,%10.10f,%10.10f,%10.10f,%10.10f\n',...
      gcoord(i,1),gcoord(i,2),gcoord(i,3),sigm4(i),ux(i),uy(i),uz(i));
%      fprintf(nodedata,'%10.10f,%10.10f,%10.10f,%10.10f\n',gcoord(i,1),gcoord(i,2),gcoord(i,3),tempfull(i));
end  

for j=1:nel
    for i=1:size(IEN,1)
        fprintf(nodedata,'%d       ',IEN(i,j));  
    end
    fprintf(nodedata,'\n');
end
end
