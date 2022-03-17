function [G,g]=gravitation(kmap,grav)
global inedge bedge elem centelem coord normals
Klef=zeros(2,2);
Krel=zeros(2,2);
G=zeros(size(elem,1),1);
for i=1:size(bedge,1)
    
   
        
        lef=bedge(i,3);  
        % tensor de permeabilidade do elemento a esquerda
        Klef(1,1)=kmap(elem(lef,5),2);
        Klef(1,2)=kmap(elem(lef,5),3);
        Klef(2,1)=kmap(elem(lef,5),4);
        Klef(2,2)=kmap(elem(lef,5),5);
        g(i,1)=dot(normals(i,1:2),(Klef*grav(lef,:)')');
        
        G(lef,1)=G(lef,1)+g(i,1);
end

for ii=1:size(inedge,1)
    
        
        lef=inedge(ii,3);
        rel=inedge(ii,4);
        % calculo do ponto meio da face        
        vm=(coord(inedge(ii,1),1:2)+coord(inedge(ii,2),1:2))*0.5;
        % calculo da distancia do centro ao ponto meio da face 
        dj1=norm(centelem(lef,1:2)-vm); 
        dj2=norm(centelem(rel,1:2)-vm);
        
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
        
        Keq=inv(dj1*inv(Klef)+dj2*inv(Krel)); % equation 21
        graveq=(dj1*grav(lef,:)+dj2*grav(rel,:))'; % equation 22
        g(ii+size(bedge,1),1)=dot(normals(ii+size(bedge,1),1:2),(Keq*graveq)'); % equation 20
        
        G(lef,1)=G(lef,1)+g(ii+size(bedge,1),1);
        G(rel,1)=G(rel,1)-g(ii+size(bedge,1),1);
end


end