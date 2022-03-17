function [E]=xicalculate(no,kmap,epsilon)
global nsurn2 nsurn1 esurn2 esurn1 elem centelem coord


O=zeros(esurn2(no+1)-esurn2(no),3); % vetor de baricentro na vizinhança do nó "ni".
nec=size(O,1);
K=zeros(3);
R=[0 1 0; -1 0 0; 0 0 0];

for k=1:nec,
    %Verifica se o elemento é um quadrilátero ou um triângulo.
    ielem= esurn1(esurn2(no)+k);
    xilem=centelem(ielem,:);
    
    vertices=elem(ielem,find(elem(ielem,1:4)~=0));
    
    arroundvertices=nsurn1(nsurn2(no) + 1:nsurn2(no + 1));
    
    [lia1,locb1]=ismember(arroundvertices, vertices);
    
    j=1;
    v=zeros(2,1);
    for i=1:length(locb1)
        if locb1(i)~=0 && i==1 && locb1(length(locb1))~=0
            v(2,1)=arroundvertices(i,1);
            j=1;
        elseif locb1(i)~=0
            v(j,1)=arroundvertices(i,1);
            j=j+1;
        end
    end
    
    v1=v(1,1);
    v2=v(2,1);
    % t 
    tsigma=coord(v1,:)-coord(no,:);
    ttau=coord(v2,:)-coord(no,:);
    % tensor de permeabilidade
    K(1,1)=kmap(elem(ielem,5),2);
    K(1,2)=kmap(elem(ielem,5),3);
    K(2,1)=kmap(elem(ielem,5),4);
    K(2,2)=kmap(elem(ielem,5),5);
    % calculo da distancia ortogonal
    dKsigma=norm(cross(xilem-coord(no,:),tsigma))/norm(tsigma);
    dKtau=norm(cross(xilem-coord(no,:),ttau))/norm(ttau);
    
    E(k,1)=dot(K*(R*-tsigma'),R*-tsigma')/dKsigma;
    E(k,2)=dot(K*(R*ttau'),R*ttau')/dKtau;
    
end

end