function [M,I] = globalmatrixtpfa(Kde, Kn, nflag, Hesq,gravresult,gravrate)

global inedge bedge elem coord centelem bcflag

% Constrói a matriz global.
% prealocação da matriz global e do vetor termo de fonte
M=zeros(size(elem,1),size(elem,1));
I=zeros(size(elem,1),1);

% Loop de faces de contorno

for ifacont=1:size(bedge,1)
    v0=coord(bedge(ifacont,2),:)-coord(bedge(ifacont,1),:);
    normcont=norm(v0);
    % calculo das constantes nas faces internas
    A=-Kn(ifacont)/(Hesq(ifacont)*norm(v0));
    
    if bedge(ifacont,5)<200
        
        c1=nflag(bedge(ifacont,1),2);
        c2=nflag(bedge(ifacont,2),2);
        
        %Preenchimento
        
        M(bedge(ifacont,3),bedge(ifacont,3))=M(bedge(ifacont,3),bedge(ifacont,3))- A*(norm(v0)^2);
        
        I(bedge(ifacont,3))=I(bedge(ifacont,3))-c1*A*(norm(v0)^2)+gravrate(ifacont,1);
        
    else
        x=bcflag(:,1)==bedge(ifacont,5);
        r=find(x==1);
        I(bedge(ifacont,3))=I(bedge(ifacont,3))- normcont*bcflag(r,2);
        
        
    end
end

for iface=1:size(inedge,1),
    norma=norm(coord(inedge(iface,1),:)-coord(inedge(iface,2),:));
    %Contabiliza as contribuições do fluxo numa faces  para os elementos %
    %a direita e a esquerda dela.                                        %
    M(inedge(iface,3), inedge(iface,3))=M(inedge(iface,3), inedge(iface,3))-Kde(iface,1);
    M(inedge(iface,3), inedge(iface,4))=M(inedge(iface,3), inedge(iface,4))+Kde(iface,1);
    M(inedge(iface,4), inedge(iface,4))=M(inedge(iface,4), inedge(iface,4))-Kde(iface,1);
    M(inedge(iface,4), inedge(iface,3))=M(inedge(iface,4), inedge(iface,3))+Kde(iface,1);
    I(inedge(iface,3))=I(inedge(iface,3))+gravrate(iface+size(bedge,1),1);
    I(inedge(iface,4))=I(inedge(iface,4))-gravrate(iface+size(bedge,1),1);
end
end
