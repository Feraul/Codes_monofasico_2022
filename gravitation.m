function [G,g]=gravitation(kmap,grav,gravface)
global inedge bedge elem centelem coord normals
Klef=zeros(2,2);
Krel=zeros(2,2);
G=zeros(size(elem,1),1);
R=[0 1 ;-1 0 ];

for i=1:size(inedge,1)+size(bedge,1)
    
    if i<size(bedge,1) || i==size(bedge,1)
        
        ve1=coord(bedge(i,2),1:2)-coord(bedge(i,1),1:2);
        lef=bedge(i,3);
        % tensor de permeabilidade do elemento a esquerda
        Klef(1,1)=kmap(elem(lef,5),2);
        Klef(1,2)=kmap(elem(lef,5),3);
        Klef(2,1)=kmap(elem(lef,5),4);
        Klef(2,2)=kmap(elem(lef,5),5);
        x=(coord(bedge(i,1),1)+coord(bedge(i,2),1))*0.5;
        y=(coord(bedge(i,1),2)+coord(bedge(i,2),2))*0.5;
        gravface1=[cos(x)*cos(y) -sin(x)*sin(y)];
        
        %gravface1=[0 sind(x)*cosd(y)/y];
        g(i,1)=dot((R*ve1')',(Klef*grav(lef,1:2)')');
        %g(i,1)=dot(normals(i,1:2),(Klef*gravface1')');
         
        %g(i,1)=gravface(i,1);
        G(lef,1)=G(lef,1)+g(i,1);
    else
        
        
        ii=i-size(bedge,1);
        
        lef=inedge(ii,3);
        rel=inedge(ii,4);
        
        % calculo do ponto meio da face
        vm=(coord(inedge(ii,1),1:2)+coord(inedge(ii,2),1:2))*0.5;
        vd1=coord(inedge(ii,2),1:2)-coord(inedge(ii,1),1:2);
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
        g(i,1)=dot((R*vd1')',(Keq*graveq)'); % equation 20
        
        G(lef,1)=G(lef,1)+g(i,1);
        G(rel,1)=G(rel,1)-g(i,1);
    end
    
    
end
end