function  solid11=h_refinement(B,knots,deg )
Xi=knots{1};Eta=knots{2};Zeta=knots{3};


[ X_ ] = getSubDivKVValues( Xi, 2); 
[B,Xi]=RefineKnotSolidXi( B,deg.p,Xi,X_ );
[ E_ ] = getSubDivKVValues( Eta,2);
[B,Eta]=RefineKnotSolidEta( B,deg.q,Eta,E_ );
[ Z_ ] =getSubDivKVValues( Zeta, 2);
[B,Zeta]=RefineKnotSolidZeta( B,deg.r,Zeta,Z_ );

solid11.knots={Xi Eta Zeta};
solid11.coefs=B;
end