function [EE]=xxicalculate(no,kmap,epsilon)
global nsurn2 nsurn1 esurn2 esurn1 elem centelem coord


O=zeros(esurn2(no+1)-esurn2(no),3); % vetor de baricentro na vizinhança do nó "ni".
P=zeros(nsurn2(no+1)-nsurn2(no),3);
nec=size(O,1);
K=zeros(3);
Ksigma=zeros(3);
Ktau=zeros(3);
R=[0 1 0; -1 0 0; 0 0 0];

if size(P,1)==size(O,1) %Se for um nó interno.
    for k=1:nec,
        
        ielem= esurn1(esurn2(no)+k);
        xilem=centelem(ielem,:);
        % primeiro elemento, por tanto precisa do proximo
        if (k==1)
            
            ielemk=esurn1(esurn2(no)+size(O,1)); % elemento anterior
            xilemk=centelem(ielemk,:);
            
            ielemkk=esurn1(esurn2(no)+k+1); % elemento posterior
            xilemkk=centelem(ielemkk,:);
            
        elseif k==nec
            
            ielemk=esurn1(esurn2(no)+k-1); % elemento anterior
            xilemk=centelem(ielemk,:);
            
            ielemkk=esurn1(esurn2(no)+1); % elemento posterior
            xilemkk=centelem(ielemkk,:);
            
        else
            ielemk=esurn1(esurn2(no)+k-1);% elemento anterior
            xilemk=centelem(ielemk,:);
            
            ielemkk=esurn1(esurn2(no)+k+1);% elemento posterior
            xilemkk=centelem(ielemkk,:);
        end
        %======================================
        vertices=elem(ielem,find(elem(ielem,1:4)~=0));% obetendo os vértices do ielem
        arroundvertices=nsurn1(nsurn2(no) + 1:nsurn2(no + 1)); % vertices em geral
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
        
        Ksigma(1,1)=kmap(elem(ielemk,5),2);
        Ksigma(1,2)=kmap(elem(ielemk,5),3);
        Ksigma(2,1)=kmap(elem(ielemk,5),4);
        Ksigma(2,2)=kmap(elem(ielemk,5),5);
        
        Ktau(1,1)=kmap(elem(ielemkk,5),2);
        Ktau(1,2)=kmap(elem(ielemkk,5),3);
        Ktau(2,1)=kmap(elem(ielemkk,5),4);
        Ktau(2,2)=kmap(elem(ielemkk,5),5);
        
        d1=norm(cross(xilemk-coord(no,:),xsigma-coord(no,:)))/norm(xsigma-coord(no,:));
        d2=norm(cross(xilem-coord(no,:),xsigma-coord(no,:)))/norm(xsigma-coord(no,:));
        d3=norm(cross(xilem-coord(no,:),xtau-coord(no,:)))/norm(xtau-coord(no,:));
        d4=norm(cross(xilemkk-coord(no,:),xtau-coord(no,:)))/norm(xtau-coord(no,:));
        
        EE(k,1)=dot(Ksigma*(R*tsigma'),R*-(xilemk-coord(no,:))')/d1;
        EE(k,2)=dot(K*(R*-tsigma'),R*(xilem-coord(no,:))')/d2;
        
        EE(k,3)=dot(K*(R*ttau'),R*-(xilem-coord(no,:))')/d3;
        EE(k,4)=dot(Ktau*(R*-ttau'),R*(xilemkk-coord(no,:))')/d4;
        
        
    end
else %Se for um nó do contorno.
    for k=1:nec,
        ielem= esurn1(esurn2(no)+k); %elemento em questão
        xilem=centelem(ielem,:); % baricentro
        if (k==1)&&(size(P,1)~=size(O,1))&& k~=nec
            
            ielemkk=esurn1(esurn2(no)+k+1);% elemento posterior
            xilemkk=centelem(ielemkk,:); % baricentro
            
            %=========================================
            vertices=elem(ielem,find(elem(ielem,1:4)~=0));
            arroundvertices=nsurn1(nsurn2(no) + 1:nsurn2(no + 1));
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
            
            v1=v(1,1);
            v2=v(2,1);
            % calculo do talfa
            tsigma=coord(v1,:)-coord(no,:);
            ttau=coord(v2,:)-coord(no,:);
            % calculo dos xalfa
            xsigma=coord(no,:)+epsilon*tsigma;
            xtau=coord(no,:)+epsilon*ttau;
            
            %======================================
            K(1,1)=kmap(elem(ielem,5),2);
            K(1,2)=kmap(elem(ielem,5),3);
            K(2,1)=kmap(elem(ielem,5),4);
            K(2,2)=kmap(elem(ielem,5),5);
            
            Ktau(1,1)=kmap(elem(ielemkk,5),2);
            Ktau(1,2)=kmap(elem(ielemkk,5),3);
            Ktau(2,1)=kmap(elem(ielemkk,5),4);
            Ktau(2,2)=kmap(elem(ielemkk,5),5);
            
            d2=norm(cross(xilem-coord(no,:),xsigma-coord(no,:)))/norm(xsigma-coord(no,:));
            d3=norm(cross(xilem-coord(no,:),xtau-coord(no,:)))/norm(xtau-coord(no,:));
            d4=norm(cross(xilemkk-coord(no,:),xtau-coord(no,:)))/norm(xtau-coord(no,:));
            
            
            EE(k,1)=0;
            EE(k,2)=dot(K*(R*-tsigma'),R*(xilem-coord(no,:))')/d2;
            
            EE(k,3)=dot(K*(R*ttau'),R*-(xilem-coord(no,:))')/d3;
            EE(k,4)=dot(Ktau*(R*-ttau'),R*(xilemkk-coord(no,:))')/d4;
            
            
        elseif (k==1)&&(size(P,1)~=size(O,1))&& k==nec
            % unico elemento no contorno
            vertices=elem(ielem,find(elem(ielem,1:4)~=0));
            
            arroundvertices=nsurn1(nsurn2(no) + 1:nsurn2(no + 1));
            
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
            
            v1=v(1,1);
            v2=v(2,1);
            % calculo do talfa
            tsigma=coord(v1,:)-coord(no,:);
            ttau=coord(v2,:)-coord(no,:);
            % calculo dos xalfa
            xsigma=coord(no,:)+epsilon*tsigma;
            xtau=coord(no,:)+epsilon*ttau;
            
            K(1,1)=kmap(elem(ielem,5),2);
            K(1,2)=kmap(elem(ielem,5),3);
            K(2,1)=kmap(elem(ielem,5),4);
            K(2,2)=kmap(elem(ielem,5),5);
            
            d1=2*norm(cross(xilem-coord(no,:),xsigma-coord(no,:)))/norm(xsigma-coord(no,:));
            
            d3=2*norm(cross(xilem-coord(no,:),xtau-coord(no,:)))/norm(xtau-coord(no,:));
            
            EE(k,1)=0;
            EE(k,2)=dot(K*(R*-tsigma'),R*(xilem-coord(no,:))')/d1;
            
            EE(k,3)=dot(K*(R*ttau'),R*-(xilem-coord(no,:))')/d3;
            EE(k,4)=0;
            
            
        elseif (k==nec)&&(size(P,1)~=size(O,1)) && k~=1
            
            ielemk=esurn1(esurn2(no)+k-1); % elemento anterior
            xilemk=centelem(ielemk,:); % baricentro
            %=======================================================
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
            % calculo do talfa
            tsigma=coord(v1,:)-coord(no,:);
            ttau=coord(v2,:)-coord(no,:);
            % calculo dos xalfa
            xsigma=coord(no,:)+epsilon*tsigma;
            xtau=coord(no,:)+epsilon*ttau;
            
            %===============================================
            K(1,1)=kmap(elem(ielem,5),2);
            K(1,2)=kmap(elem(ielem,5),3);
            K(2,1)=kmap(elem(ielem,5),4);
            K(2,2)=kmap(elem(ielem,5),5);
            
            Ksigma(1,1)=kmap(elem(ielemk,5),2);
            Ksigma(1,2)=kmap(elem(ielemk,5),3);
            Ksigma(2,1)=kmap(elem(ielemk,5),4);
            Ksigma(2,2)=kmap(elem(ielemk,5),5);
            
            d1=norm(cross(xilemk-coord(no,:),xsigma-coord(no,:)))/norm(xsigma-coord(no,:));
            d2=norm(cross(xilem-coord(no,:),xsigma-coord(no,:)))/norm(xsigma-coord(no,:));
            d3=norm(cross(xilem-coord(no,:),xtau-coord(no,:)))/norm(xtau-coord(no,:));
            
            
            EE(k,1)=dot(Ksigma*(R*tsigma'),R*-(xilemk-coord(no,:))')/d1;
            EE(k,2)=dot(K*(R*-tsigma'),R*(xilem-coord(no,:))')/d2;
            
            EE(k,3)=dot(K*(R*ttau'),R*-(xilem-coord(no,:))')/d3;
            EE(k,4)=0;
            
        else
            ielemk=esurn1(esurn2(no)+k-1); % elemento anterior
            xilemk=centelem(ielemk,:);% baricentro
            
            ielemkk=esurn1(esurn2(no)+k+1); % elemento posterior
            xilemkk=centelem(ielemkk,:); % baricentro
            %===========================================
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
            
            tsigma=coord(v1,:)-coord(no,:);
            ttau=coord(v2,:)-coord(no,:);
            
            xsigma=coord(no,:)+epsilon*tsigma;
            xtau=coord(no,:)+epsilon*ttau;
            
            %==================================
            
            K(1,1)=kmap(elem(ielem,5),2);
            K(1,2)=kmap(elem(ielem,5),3);
            K(2,1)=kmap(elem(ielem,5),4);
            K(2,2)=kmap(elem(ielem,5),5);
            
            Ksigma(1,1)=kmap(elem(ielemk,5),2);
            Ksigma(1,2)=kmap(elem(ielemk,5),3);
            Ksigma(2,1)=kmap(elem(ielemk,5),4);
            Ksigma(2,2)=kmap(elem(ielemk,5),5);
            
            Ktau(1,1)=kmap(elem(ielemkk,5),2);
            Ktau(1,2)=kmap(elem(ielemkk,5),3);
            Ktau(2,1)=kmap(elem(ielemkk,5),4);
            Ktau(2,2)=kmap(elem(ielemkk,5),5);
            
            d1=norm(cross(xilemk-coord(no,:),xsigma-coord(no,:)))/norm(xsigma-coord(no,:));
            d2=norm(cross(xilem-coord(no,:),xsigma-coord(no,:)))/norm(xsigma-coord(no,:));
            d3=norm(cross(xilem-coord(no,:),xtau-coord(no,:)))/norm(xtau-coord(no,:));
            d4=norm(cross(xilemkk-coord(no,:),xtau-coord(no,:)))/norm(xtau-coord(no,:));
            
            EE(k,1)=dot(Ksigma*(R*tsigma'),R*-(xilemk-coord(no,:))')/d1;
            EE(k,2)=dot(K*(R*-tsigma'),R*(xilem-coord(no,:))')/d2;
            
            EE(k,3)=dot(K*(R*ttau'),R*-(xilem-coord(no,:))')/d3;
            EE(k,4)=dot(Ktau*(R*-ttau'),R*(xilemkk-coord(no,:))')/d4;
            
        end
        
    end
end
end