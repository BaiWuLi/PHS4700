%
% Definir terrain
%
global longueur largeur hauteur hfilet lfilet dir ;
bandes=5/100;
epaisseur=10/100;
lo2=longueur/2;
xoff=0.0;
xcp=3*longueur/8;
ycp=3*largeur/8;
ba2=bandes/2;
la2=largeur/2;
yoff=0.0;
epp=3/100;
df=(lfilet-largeur)/2;
LigneX1=[xoff    xoff+longueur    xoff+longueur   xoff   xoff;
         -ba2+la2   -ba2+la2    ba2+la2    ba2+la2   -ba2+la2 ;
          hauteur     hauteur     hauteur     hauteur      hauteur ];
LigneX2=[xoff    xoff+longueur    xoff+longueur   xoff   xoff;
         yoff   yoff   yoff+bandes   yoff+bandes    yoff ;
          hauteur     hauteur     hauteur     hauteur      hauteur ];
LigneX3=[xoff    xoff+longueur    xoff+longueur   xoff   xoff;
          yoff+largeur-bandes    yoff+largeur-bandes   yoff+largeur   yoff+largeur    yoff+largeur-bandes ;
          hauteur     hauteur     hauteur     hauteur      hauteur ];
LigneY1=[xoff    xoff+bandes   xoff+bandes   xoff   xoff;
         yoff   yoff    yoff+largeur    yoff+largeur   yoff ;
          hauteur     hauteur     hauteur     hauteur      hauteur ];
LigneY2=[xoff+longueur-bandes    xoff+longueur     xoff+longueur   xoff+longueur-bandes   xoff+longueur-bandes;
         yoff   yoff    yoff+largeur    yoff+largeur   yoff ;
          hauteur     hauteur     hauteur     hauteur      hauteur ];
TableX1=[xoff+bandes    xoff+longueur-bandes    xoff+longueur-bandes   xoff+bandes   xoff+bandes;
         yoff+bandes   yoff+bandes    -ba2    -ba2    yoff+bandes ;
          hauteur     hauteur     hauteur     hauteur      hauteur ];
TableX2=[xoff+bandes    xoff+longueur-bandes    xoff+longueur-bandes   xoff+bandes   xoff+bandes;
          ba2    ba2    yoff+largeur-bandes    yoff+largeur-bandes    ba2 ;
          hauteur     hauteur     hauteur     hauteur      hauteur ];
CoteX1=[xoff    xoff+longueur    xoff+longueur   xoff   xoff;
         yoff   yoff   yoff   yoff   yoff ;
          hauteur-epaisseur     hauteur-epaisseur     hauteur     hauteur      hauteur-epaisseur ];
CoteY1=[xoff    xoff   xoff   xoff   xoff;
         yoff    yoff+largeur    yoff+largeur   yoff   yoff ;
          hauteur-epaisseur   hauteur-epaisseur     hauteur     hauteur      hauteur-epaisseur ];
Patte1F=[-xcp+epp+lo2 -xcp-epp+lo2 -xcp-epp+lo2 -xcp+epp+lo2 -xcp+epp+lo2;
         -ycp+la2-epp -ycp+la2-epp -ycp+la2-epp -ycp+la2-epp -ycp+la2-epp;
         0 0 hauteur hauteur 0];
Patte1C=[-xcp-epp+lo2 -xcp-epp+lo2 -xcp-epp+lo2 -xcp-epp+lo2 -xcp-epp+lo2;
         -ycp+la2-epp -ycp+la2+epp -ycp+la2+epp -ycp+la2-epp -ycp+la2-epp;
         0 0 hauteur hauteur 0];
Patte2F=[ xcp+epp+lo2  xcp-epp+lo2  xcp-epp+lo2  xcp+epp+lo2  xcp+epp+lo2;
         -ycp+la2-epp -ycp+la2-epp -ycp+la2-epp -ycp+la2-epp -ycp+la2-epp;
         0 0 hauteur hauteur 0];
Patte2C=[ xcp-epp+lo2  xcp-epp+lo2  xcp-epp+lo2  xcp-epp+lo2  xcp-epp+lo2;
         -ycp+la2-epp -ycp+la2+epp -ycp+la2+epp -ycp+la2-epp -ycp+la2-epp;
         0 0 hauteur hauteur 0];
Patte3F=[-xcp+epp+lo2 -xcp-epp+lo2 -xcp-epp+lo2 -xcp+epp+lo2 -xcp+epp+lo2;
          ycp+la2-epp  ycp+la2-epp  ycp+la2-epp  ycp+la2-epp  ycp+la2-epp;
         0 0 hauteur hauteur 0];
Patte3C=[-xcp-epp+lo2 -xcp-epp+lo2 -xcp-epp+lo2 -xcp-epp+lo2 -xcp-epp+lo2;
          ycp+la2-epp  ycp+la2+epp  ycp+la2+epp  ycp+la2-epp  ycp+la2-epp;
         0 0 hauteur hauteur 0];
Patte4F=[ xcp+epp+lo2  xcp-epp+lo2  xcp-epp+lo2  xcp+epp+lo2  xcp+epp+lo2;
          ycp+la2-epp  ycp+la2-epp  ycp+la2-epp  ycp+la2-epp  ycp+la2-epp;
         0 0 hauteur hauteur 0];
Patte4C=[ xcp-epp+lo2  xcp-epp+lo2  xcp-epp+lo2  xcp-epp+lo2  xcp-epp+lo2;
          ycp+la2-epp  ycp+la2+epp  ycp+la2+epp  ycp+la2-epp  ycp+la2-epp;
         0 0 hauteur hauteur 0];
Filet=[0+lo2 0+lo2 0+lo2 0+lo2 0+lo2;
       yoff-df yoff-df yoff+largeur+df yoff+largeur+df yoff-df;
        hauteur hauteur+hfilet hauteur+hfilet hauteur hauteur];
fill3(LigneX1(1,:),LigneX1(2,:),LigneX1(3,:),[1 1 1]);
fill3(LigneX2(1,:),LigneX2(2,:),LigneX2(3,:),[1 1 1]);
fill3(LigneX3(1,:),LigneX3(2,:),LigneX3(3,:),[1 1 1]);
fill3(LigneY1(1,:),LigneY1(2,:),LigneY1(3,:),[1 1 1]);
fill3(LigneY2(1,:),LigneY2(2,:),LigneY2(3,:),[1 1 1]);
fill3(TableX1(1,:),TableX1(2,:),TableX1(3,:),[0.5 1 0.5]);
fill3(TableX2(1,:),TableX2(2,:),TableX2(3,:),[0.5 1 0.5]);
fill3(CoteX1(1,:),CoteX1(2,:),CoteX1(3,:),[0.5 1 0.5]);
fill3(CoteY1(1,:),CoteY1(2,:),CoteY1(3,:),[0.5 1 0.5]);
fill3(Patte1F(1,:),Patte1F(2,:),Patte1F(3,:),[0 0 0]);
fill3(Patte1C(1,:),Patte1C(2,:),Patte1C(3,:),[0 0 0]);
fill3(Patte2F(1,:),Patte2F(2,:),Patte2F(3,:),[0 0 0]);
fill3(Patte2C(1,:),Patte2C(2,:),Patte2C(3,:),[0 0 0]);
fill3(Patte3F(1,:),Patte3F(2,:),Patte3F(3,:),[0 0 0]);
fill3(Patte3C(1,:),Patte3C(2,:),Patte3C(3,:),[0 0 0]);
fill3(Patte4F(1,:),Patte4F(2,:),Patte4F(3,:),[0 0 0]);
fill3(Patte4C(1,:),Patte4C(2,:),Patte4C(3,:),[0 0 0]);
fill3(Filet(1,:),Filet(2,:),Filet(3,:),[1 1 0]);
