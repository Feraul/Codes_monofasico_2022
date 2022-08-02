function [gravresult, gravrate]=gravitation(kmap,vec_gravelem, vec_gravface)

global inedge bedge elem centelem  coord normals

Klef=zeros(2,2);
Krel=zeros(2,2);
gravrate=zeros(size(bedge,1)+size(inedge,1),1);
gravresult=zeros(size(elem,1),1);
RR=[0 1 ; -1 0];
K1=zeros(3,3);
K2=zeros(3,3);
K=zeros(3,3);
i=1;
for ifacont=1 : size(bedge,1)
    
    lef = bedge(ifacont,3);
    
    
    % calculo do ponto meio da face
    % Distância do centro da face ao centro do elemento
    %ve1 = coord(bedge(ifacont,2),1:2)-coord(bedge(ifacont,1),1:2);
    %ve2 = ve1 - centelem(lef,1:2);
    %dj1 = norm(ve2);
    centface=0.5*(coord(bedge(ifacont,2),1:2)+coord(bedge(ifacont,1),1:2));
    % tensor de permeabilidade do elemento a esquerda
    Klef(1,1)=kmap(elem(lef,5),2);
    Klef(1,2)=kmap(elem(lef,5),3);
    Klef(2,1)=kmap(elem(lef,5),4);
    Klef(2,2)=kmap(elem(lef,5),5);
    
    
    % for edge faces, the permeability distance weighted harmonic mean
    % is equal to the permeability tensor
    iK = Klef;
    
    % normal da face
    n = abs(normals(ifacont,1:2));
    
    % for edge faces, the weighted arithmetic mean of element center
    % gravities is equal to the gravities
    fg = vec_gravelem(lef,1:2);
    
    gravrate(ifacont,1) = n*iK*fg';
    
    if(bedge(ifacont,5) > 200)
        kg(i,1) = Klef(1,1)*vec_gravface(ifacont,1)+Klef(1,2)*vec_gravface(ifacont,2);
        kg(i,2) = Klef(2,1)*vec_gravface(ifacont,1)+Klef(2,2)*vec_gravface(ifacont,2);
        
        
        gravrate(ifacont,1) = abs(normals(ifacont,1))*kg(i,1) + abs(normals(ifacont,2))*kg(i,2);
        
        
        %if(centface(ifacont,2)==0)
        %    gravrate(ifacont,1)=gravrate(ifacont,1)*(-1);
        %end
        %i=i+1;
        
        %disp(gravrate(ifacont,1))
    end
    
    %gravrate(ifacont,1)=dot((RR*ve1')',(Klef*grav(lef,1:2)')');
    
    if (centface(1,1)==0 | centface(1,2)==0)
        gravresult(lef,1)=gravresult(lef,1)+gravrate(ifacont,1);
    else
        gravresult(lef,1)=gravresult(lef,1)-gravrate(ifacont,1);
    end
    
end

for ifacont=1 : size(inedge,1)
    
    iaux = ifacont + size(bedge,1);
    
    lef=inedge(ifacont,3);
    rel=inedge(ifacont,4);
    centface=0.5*(coord(inedge(ifacont,2),:)+coord(inedge(ifacont,1),:));
    % Distance from edge center to left element center
    dj1_vec = centface(1,:)-centelem(lef,:); %vetor distância
    dj1 = norm(dj1_vec); %módulo da distância
    
    % Distance from edge center to right element center
    dj2_vec = centface(1,:)-centelem(rel,:); %vetor distância
    dj2 = norm(dj2_vec); %módulo da distância
    
    
    % calculo do ponto meio da face
    %vd1aux = coord(inedge(iface,2),:)-coord(inedge(iface,1),:);
    %vd1=coord(inedge(iface,2),1:2)-coord(inedge(iface,1),1:2);
    
    % Do ponto do início da aresta até o centro da célula da direita
    %vd2=centelem(rel,:)-coord(inedge(iface,1),:);
    %cd=cross(vd1aux,vd2);
    %dj1=norm(cd)/norm(vd1aux); % altura a direita
    
    %Do ponto do início da aresta até o centro da célula da esquerda
    %ve2=centelem(lef,:)-coord(inedge(iface,1),:);
    %ce=cross(vd1aux,ve2);
    %dj2=norm(ce)/norm(vd1aux); % altura a esquerda
    
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
    
    
    % invert element-centers permeability tensor
    iKlef = inv(Klef);
    iKrel = inv(Krel);
    
    % take distance weighted harmonic mean of permeability
    
    iK = inv((dj1*iKlef + dj2*iKrel)/(dj1+dj2));
    %inv(inv(Klef/dj1)+inv(Krel/dj2)); % equation 21
    
    % normal da face
    n = abs(normals(iaux,1:2));
    
    % take the weighted arithmetic mean of element center gravities
    fg = (dj1*vec_gravelem(lef,1:2)+dj2*vec_gravelem(rel,1:2))/(dj1+dj2);
    
    gravrate(iaux,1) = n*iK*fg';
    
    %kg(iaux,1) = Klef(1,1)*vec_gravface(iaux,1)+Klef(1,2)*vec_gravface(iaux,2);
    %kg(iaux,2) = Klef(2,1)*vec_gravface(iaux,1)+Klef(2,2)*vec_gravface(iaux,2);
    
    
    %graveq = ((dj1*grav(lef,:)+dj2*grav(rel,:))'); %equation 22
    %gravrate(iface+size(bedge,1),1)=dot((RR*vd1')',(Keq*graveq(1:2))');
    
    if (coord(inedge(ifacont,1),1)==1 | coord(inedge(ifacont,1),2)==0)
        gravresult(lef,1)=gravresult(lef,1) - gravrate(iaux,1);
        gravresult(rel,1)=gravresult(rel,1) + gravrate(iaux,1);
    else
        gravresult(lef,1)=gravresult(lef,1) + gravrate(iaux,1);
        gravresult(rel,1)=gravresult(rel,1) - gravrate(iaux,1);
    end
    
end

end