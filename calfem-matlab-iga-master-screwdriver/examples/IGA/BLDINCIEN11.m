function [ INC,IEN,INB] = BLDINCIEN11( deg,n,m,l )
p= deg.p ; q=deg.q; r=deg.r ;
A = 0;
    nel = (n-p) * (m-q) * (l-r); % number of elements
    nnp = n*m*l; % number of global basis functions
    nen = (p+1)*(q+1)*(r+1); % number of local basis functions
    INC = zeros([nnp],[3]); % NURBS coordinates array
    IEN = zeros([nen],[nel]); % connectivity array ( = Enod' )
    INB=zeros(n,m,l);
for k = 1 : l
        for j = 1 : m
            for i = 1 : n
                
                A = A + 1;
                INC(A,1) = i;
                INC(A,2) = j;
                INC(A,3) = k;
                INB(i, j, k)=A;
            end
        end
end
e=0;bi=0;
for k = 1 : (l-r+1)
        for j = 1 : (m-q)
            for i = 1 : (n-p)
                e=e+1;
               for kloc = 0 : r   %1 ,  3
                   for jloc = 0 : q   %2 ,  12
                           for iloc = 0 : p    %3 ,  55                     
                              if (kloc+k)==(l+1)
                                B_j=INB((iloc+i), jloc+j,1); 
                               else
                                B_j=INB((iloc+i), jloc+j, kloc+k);   % global function number
                              end
                              b = (p+1)*(q+1)*(r+1)-kloc*(p+1)*(q+1) - jloc*(p+1)- iloc ; % local function number
                              IEN(b,e) = B_j; % assign connectivity
                           end
                   end
               end
             
            end
        end
end

end
