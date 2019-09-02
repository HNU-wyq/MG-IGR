function  alfa= Judge(rr,Xi,Xii,i,j)
if abs(Xi(i)==Xi(i+rr))
            alfa=0; 
        else
            alfa = (Xii(j+rr-1) - Xi(i))/(Xi(i+rr) - Xi(i));
end
end