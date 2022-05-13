function [G,g]=gravitation(kmap,grav,gravface,gravelem)
global inedge bedge elem centelem coord normals
Klef=zeros(2,2);
Krel=zeros(2,2);
G=zeros(size(elem,1),1);
R=[0 1 ;-1 0 ];
K1=zeros(3,3);
K2=zeros(3,3);
K=zeros(3,3);
for ifacont=1:size(bedge,1)
    lef=bedge(ifacont,3);
    ve1=coord(bedge(ifacont,2),1:2)-coord(bedge(ifacont,1),1:2);
    
    % tensor de permeabilidade do elemento a esquerda
    Klef(1,1)=kmap(elem(lef,5),2);
    Klef(1,2)=kmap(elem(lef,5),3);
    Klef(2,1)=kmap(elem(lef,5),4);
    Klef(2,2)=kmap(elem(lef,5),5);
   
    g(ifacont,1)=dot((R*ve1')',(Klef*grav(lef,1:2)')');
   
    G(lef,1)=G(lef,1)+g(ifacont,1);

end
for iface=1:size(inedge,1)
   
         lef=inedge(iface,3);
         rel=inedge(iface,4);
         
         % calculo do ponto meio da face
         %vm=(coord(inedge(ii,1),1:2)+coord(inedge(ii,2),1:2))*0.5;
         vd1aux=coord(inedge(iface,2),:)-coord(inedge(iface,1),:);
         vd1=coord(inedge(iface,2),1:2)-coord(inedge(iface,1),1:2);
         % calculo da distancia do centro ao ponto meio da face
         %dj1=norm(centelem(lef,1:2)-vm);
         %dj2=norm(centelem(rel,1:2)-vm);
         
         %Do ponto do início da aresta até o centro da célula da direita
         vd2=centelem(rel,:)-coord(inedge(iface,1),:);
         cd=cross(vd1aux,vd2);
         dj1=norm(cd)/norm(vd1aux); % altura a direita
         
         %Do ponto do início da aresta até o centro da célula da direita
         ve2=centelem(lef,:)-coord(inedge(iface,1),:);
         ce=cross(vd1aux,ve2);
         dj2=norm(ce)/norm(vd1aux); % altura a esquerda
         
         % tensor do elemento a esquerda
         
         Klef(1,1)=kmap(elem(lef,5),2);
         Klef(1,2)=kmap(elem(lef,5),3);
         Klef(2,1)=kmap(elem(lef,5),4);
         Klef(2,2)=kmap(elem(lef,5),5);
         
         % tensor do elemento a direita
         
         Krel(1,1)=kmap(elem(rel,5),2);
         Krel(1,2)=kmap(elem(rel,5),3);
         Krel(2,1)=kmap(elem(rel,5),4);
         Krel(2,2)=kmap(elem(rel,5),5);
         
         Keq=inv(inv(Klef/dj1)+inv(Krel/dj2)); % equation 21
         graveq=((dj1*grav(lef,:)+dj2*grav(rel,:))'); % equation 22
         g(iface+size(bedge,1),1)=dot((R*vd1')',(Keq*graveq)'); % equation 20
         
         G(lef,1)=G(lef,1)+g(iface+size(bedge,1),1);
         G(rel,1)=G(rel,1)-g(iface+size(bedge,1),1);
     
end
    
    
end
