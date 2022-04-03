function [M,I] = globalmatrixtpfa(Kde, Kn, nflagface, Hesq,gravresult,gravrate,gravno,gravelem,gravitational)

global inedge bedge elem coord centelem bcflag

% Constrói a matriz global.
% prealocação da matriz global e do vetor termo de fonte
M=zeros(size(elem,1),size(elem,1));
I=zeros(size(elem,1),1);

% Loop de faces de contorno
m=0;
m3=0;
for ifacont=1:size(bedge,1)
    v0=coord(bedge(ifacont,2),:)-coord(bedge(ifacont,1),:);
    normcont=norm(v0);
    lef=bedge(ifacont,3);
    % calculo das constantes nas faces internas
    A=-Kn(ifacont)/(Hesq(ifacont)*norm(v0));
    
    if bedge(ifacont,5)<200
        
        c1=nflagface(ifacont,2);
        
        if strcmp(gravitational,'yes')
            m1=gravno(bedge(ifacont,1),1);
            
            m=A*(norm(v0)^2*m1-norm(v0)^2*gravelem(lef));
        end
        
        %Preenchimento
        
        M(bedge(ifacont,3),bedge(ifacont,3))=M(bedge(ifacont,3),bedge(ifacont,3))- A*(norm(v0)^2);
        
        I(bedge(ifacont,3))=I(bedge(ifacont,3))-c1*A*(norm(v0)^2)-m;
        
    else
        x=bcflag(:,1)==bedge(ifacont,5);
        r=find(x==1);
        I(bedge(ifacont,3))=I(bedge(ifacont,3))- normcont*bcflag(r,2);
        
        
    end
    
end

for iface=1:size(inedge,1),
     lef=inedge(iface,3);
    rel=inedge(iface,4);
    %Contabiliza as contribuições do fluxo numa faces  para os elementos %
    %a direita e a esquerda dela.                                        %
    M(inedge(iface,3), inedge(iface,3))=M(inedge(iface,3), inedge(iface,3))-Kde(iface,1);
    M(inedge(iface,3), inedge(iface,4))=M(inedge(iface,3), inedge(iface,4))+Kde(iface,1);
    M(inedge(iface,4), inedge(iface,4))=M(inedge(iface,4), inedge(iface,4))-Kde(iface,1);
    M(inedge(iface,4), inedge(iface,3))=M(inedge(iface,4), inedge(iface,3))+Kde(iface,1);
    if strcmp(gravitational,'yes')
        m3= Kde(iface)*(gravelem(rel,1)-gravelem(lef,1));
        I(inedge(iface,3))=I(inedge(iface,3))-m3;
        I(inedge(iface,4))=I(inedge(iface,4))+m3;
        
    end
end
end
