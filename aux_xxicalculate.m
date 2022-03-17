function [EE]=aux_xxicalculate(no,kmap,epsilon)
global nsurn2 nsurn1 esurn2 esurn1 elem centelem coord


O=zeros(esurn2(no+1)-esurn2(no),3); % vetor de baricentro na vizinhança do nó "ni".
P=zeros(nsurn2(no+1)-nsurn2(no),3);
nec=size(O,1);
K=zeros(3);
R=[0 1 0; -1 0 0; 0 0 0];

for k=1:nec,
    
    ielem= esurn1(esurn2(no)+k);
    xilem=centelem(ielem,:);
    
    %======================================
    vertices=elem(ielem,find(elem(ielem,1:4)~=0));% obetendo os vértices do ielem
    arroundvertices=nsurn1(nsurn2(no) + 1:nsurn2(no + 1)); % vertices em geral
    
    [v]=reordenatevertice(arroundvertices,vertices);
    v1=v(1,1);
    v2=v(2,1);
    % calculo do talfa
    tsigma=coord(v1,:)-coord(no,:);
    ttau=coord(v2,:)-coord(no,:);
    % calculo dos xalfa
    xsigma=coord(no,:)+epsilon*tsigma;
    xtau=coord(no,:)+epsilon*ttau;
    
    %==================================
    
    K(1,1)=kmap(elem(ielem,5),2);
    K(1,2)=kmap(elem(ielem,5),3);
    K(2,1)=kmap(elem(ielem,5),4);
    K(2,2)=kmap(elem(ielem,5),5);
    
    d2=norm(cross(xilem-coord(no,:),xsigma-coord(no,:)))/norm(xsigma-coord(no,:));
    d3=norm(cross(xilem-coord(no,:),xtau-coord(no,:)))/norm(xtau-coord(no,:));
    
    EE(k,1)=dot(K*(R*tsigma'),R*(xilem-coord(no,:))')/d2;
    
    EE(k,2)=dot(K*(R*ttau'),R*(xilem-coord(no,:))')/d3;
    
end

end

function[v]=reordenatevertice(arroundvertices,vertices)
       
    [lia1,locb1]=ismember(arroundvertices,vertices);
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
end