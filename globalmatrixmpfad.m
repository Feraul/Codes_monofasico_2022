function [ M, I ] = globalmatrixmpfad( w,s, Kde, Ded, Kn, Kt, nflagno, ...
    Hesq,fonte,gravresult,gravrate,gravno,gravelem)

global coord elem esurn1 esurn2  bedge inedge  centelem bcflag gravitational

%-----------------------inicio da rutina ----------------------------------%
%Constr�i a matriz global.

M=sparse(size(elem,1),size(elem,1)); %Prealoca��o de M.
I=sparse(size(elem,1),1);
% fonte
I=I+fonte;
% contribui��o dos po�os
m=0;

for ifacont=1:size(bedge,1)
    lef=bedge(ifacont,3);
    
    v0=coord(bedge(ifacont,2),:)-coord(bedge(ifacont,1),:); %fase.
    v1=centelem(bedge(ifacont,3),:)-coord(bedge(ifacont,1),:);
    v2=centelem(bedge(ifacont,3),:)-coord(bedge(ifacont,2),:);
    normcont=norm(v0);
    
    % Tratamento do n� nos v�rtices 2 e 4%
      A=-Kn(ifacont)/(Hesq(ifacont)*norm(v0));  
    if bedge(ifacont,5)<200
        c1=nflagno(bedge(ifacont,1),2);
        c2=nflagno(bedge(ifacont,2),2);
       
        %Preenchimento
        
        if strcmp(gravitational,'yes')
            m1=gravno(bedge(ifacont,1),1);
            m2=gravno(bedge(ifacont,2),1);
            m=A*(dot(v2,-v0)*m1+dot(v1,v0)*m2-norm(v0)^2*gravelem(lef))-(m2-m1)*Kt(ifacont);
           %m=gravrate(ifacont);
        else
            m=0;
        end
        
        M(bedge(ifacont,3),bedge(ifacont,3))=M(bedge(ifacont,3),bedge(ifacont,3))-A*(norm(v0)^2);
        
        I(bedge(ifacont,3))=I(bedge(ifacont,3))-A*(dot(v2,-v0)*c1+dot ...
            (v1,v0)*c2)+(c2-c1)*Kt(ifacont)-m;
        
    else
       
%                 if strcmp(gravitational,'yes')
%                     no1=bedge(ifacont,1);
%                     no2=bedge(ifacont,2);
%                     nec1=esurn2(no1+1)-esurn2(no1);
%                     nec2=esurn2(no2+1)-esurn2(no2);
%                     g1=0;
%                     if nflagno(no1,1)<200
%                         g1=gravno(no1,1);
%                     else
%                         for j=1:nec1
%                             element1=esurn1(esurn2(no1)+j);
%                             g1=g1+w(esurn2(no1)+j)*gravelem(element1);
%                         end
%                     end
%                     g2=0;
%                     if  nflagno(no2,1)<200
%                         g2=gravno(no2,1);
%                     else
%                         for j=1:nec2
%                             element2=esurn1(esurn2(no2)+j);
%                             g2=g2+w(esurn2(no2)+j)*gravelem(element2);
%                         end
%                     end
%         
%                     m=A*(dot(v2,-v0)*g1+dot(v1,v0)*g2-norm(v0)^2*gravelem(lef))-(g2-g1)*Kt(ifacont);
%                     % m=gravrate(ifacont);
%                 else
%                     m=0;
%                 end
       
        % contorno de Neumann
        x=bcflag(:,1)==bedge(ifacont,5);
        r=find(x==1);
        I(bedge(ifacont,3))=I(bedge(ifacont,3)) -normcont*bcflag(r,2);%+m;
    end
   
end


% contribui��o nas faces internas
for iface=1:size(inedge,1)
    lef=inedge(iface,3);
    rel=inedge(iface,4);
    %Contabiliza as contribui��es do fluxo numa aresta para os elementos %
    %a direita e a esquerda dela.                                        %
    
    M(inedge(iface,3), inedge(iface,3))=M(inedge(iface,3), inedge(iface,3))- Kde(iface);
    M(inedge(iface,3), inedge(iface,4))=M(inedge(iface,3), inedge(iface,4))+ Kde(iface);
    M(inedge(iface,4), inedge(iface,4))=M(inedge(iface,4), inedge(iface,4))- Kde(iface);
    M(inedge(iface,4), inedge(iface,3))=M(inedge(iface,4), inedge(iface,3))+ Kde(iface);
    
    %Se os n�s das arestas estiverem em fronteiras de Dirichlet, suas
    %contribui��es ser�o contabilizadas logo abaixo.
    
    if nflagno(inedge(iface,1),1)<200
        I(inedge(iface,3))=I(inedge(iface,3))-Kde(iface)*Ded(iface)*nflagno(inedge(iface,1),2);
        I(inedge(iface,4))=I(inedge(iface,4))+Kde(iface)*Ded(iface)*nflagno(inedge(iface,1),2);
    end
    if nflagno(inedge(iface,2),1)<200
        I(inedge(iface,3))=I(inedge(iface,3))+Kde(iface)*Ded(iface)*nflagno(inedge(iface,2),2);
        I(inedge(iface,4))=I(inedge(iface,4))-Kde(iface)*Ded(iface)*nflagno(inedge(iface,2),2);
    end
    % quando o n� pertece ao contorno de Neumann
    if nflagno(inedge(iface,1),1)==202
        
        I(inedge(iface,3))=I(inedge(iface,3))-Kde(iface)*Ded(iface)*s(inedge(iface,1)); %ok
        
        I(inedge(iface,4))=I(inedge(iface,4))+Kde(iface)*Ded(iface)*s(inedge(iface,1)); %ok
    end
    if nflagno(inedge(iface,2),1)==202
        
        I(inedge(iface,3))=I(inedge(iface,3))+Kde(iface)*Ded(iface)*s(inedge(iface,2)); %ok
        
        I(inedge(iface,4))=I(inedge(iface,4))-Kde(iface)*Ded(iface)*s(inedge(iface,2)); %ok
        
    end
    
    %Contabiliza��o das contribui��es dos n�s que n�o est�o na
    %fronteiras de Dirichlet.
    
    if nflagno(inedge(iface,1),1)>200
        for j=1:(esurn2(inedge(iface,1)+1)-esurn2(inedge(iface,1)))
            
            post_cont=esurn2(inedge(iface,1))+j;
            
            M(inedge(iface,3), esurn1(post_cont))=M(inedge(iface,3),esurn1(post_cont)) + Kde(iface)*Ded(iface)*w(post_cont);
            
            M(inedge(iface,4), esurn1(post_cont))=M(inedge(iface,4),esurn1(post_cont)) - Kde(iface)*Ded(iface)*w(post_cont);
            
        end
    end
    if nflagno(inedge(iface,2),1)>200
        for j=1:(esurn2(inedge(iface,2)+1)-esurn2(inedge(iface,2))),
            
            post_cont=esurn2(inedge(iface,2))+j;
            
            M(inedge(iface,3), esurn1(post_cont))=M(inedge(iface,3),esurn1(post_cont)) - Kde(iface)*Ded(iface)*w(post_cont);
            
            M(inedge(iface,4), esurn1(post_cont))=M(inedge(iface,4),esurn1(post_cont)) + Kde(iface)*Ded(iface)*w(post_cont);
        end
    end
    if strcmp(gravitational,'yes')
        no1=inedge(iface,1);
        no2=inedge(iface,2);
        nec1=esurn2(no1+1)-esurn2(no1);
        nec2=esurn2(no2+1)-esurn2(no2);
        g1=0;
        if nflagno(no1,1)<200
         g1=gravno(no1,1);   
        else
            for j=1:nec1
                element1=esurn1(esurn2(no1)+j);
                g1=g1+w(esurn2(no1)+j)*gravelem(element1);
            end
        end
        g2=0;
        if nflagno(no2,1)<200
          g2=gravno(no2,1);  
        else
            for jj=1:nec2
                element2=esurn1(esurn2(no2)+jj);
                g2=g2+w(esurn2(no2)+jj)*gravelem(element2);
            end
        end
        
        m= Kde(iface)*(gravelem(rel,1)-gravelem(lef,1)-Ded(iface)*(g2-g1));
        
        %m=gravrate(size(bedge,1)+iface,1);
        I(inedge(iface,3))=I(inedge(iface,3))-m;
        I(inedge(iface,4))=I(inedge(iface,4))+m;
        
    end
end

