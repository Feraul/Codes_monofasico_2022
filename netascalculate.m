function [netas,E]=netascalculate(no,kmap,epsilon)
global nsurn2 nsurn1 esurn2 esurn1 elem coord centelem


O=zeros(esurn2(no+1)-esurn2(no),3); % vetor de baricentro na vizinhança do nó "ni".
P=zeros(nsurn2(no+1)-nsurn2(no),3);
nec=size(O,1);

K=zeros(3);
Kigma=zeros(3);
Ktau=zeros(3);
R=[0 1 0; -1 0 0; 0 0 0];

if size(P,1)==size(O,1) %Se for um nó interno.
    for k=1:nec,
        ielem= esurn1(esurn2(no)+k); % elemento em questão 
        
        if (k==1)&&(size(P,1)==size(O,1))
            
            ielemk=esurn1(esurn2(no)+size(O,1)); % elemento anterior
            ielemkk=esurn1(esurn2(no)+k+1); % elemento posterior
            
        elseif (k==nec)&&(size(P,1)==size(O,1))
            
            ielemk=esurn1(esurn2(no)+k-1); % elemento anterior
            ielemkk=esurn1(esurn2(no)+1); % elemento posterior
            
        else
            ielemk=esurn1(esurn2(no)+k-1); % elemento anterior
            ielemkk=esurn1(esurn2(no)+k+1); % elemento posterior
        end
        
        
        % calculando vertices para ielem
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
        %% ================================================================
        xilem=centelem(ielem,:);
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
        
        %% ================================================================
        tf1=coord(v1,:)-coord(no,:);
        tf2=coord(v2,:)-coord(no,:);
        
        xm1=coord(no,:)+epsilon*tf1;
        xm2=coord(no,:)+epsilon*tf2;
        
        vec1=xm1-coord(no,:);
        vec2=xm2-coord(no,:);
        
        normal1=(R*(xm2-xm1)');
        
        normsigma= norm(coord(no,:)-coord(v1,:));
        normtau= norm(coord(no,:)-coord(v2,:));
        
        K(1,1)=kmap(elem(ielem,5),2);
        K(1,2)=kmap(elem(ielem,5),3);
        K(2,1)=kmap(elem(ielem,5),4);
        K(2,2)=kmap(elem(ielem,5),5);
        
        sielem=norm(cross(xm1-coord(no,:),xm2-coord(no,:)))*0.5;
        % para elemento anterior===========================================
        
        vertices1a=elem(ielemk,find(elem(ielemk,1:4)~=0));
        
        
        [lia1,locb1]=ismember(arroundvertices,vertices1a);
        j=1;
        vaa=zeros(2,1);
        for i=1:length(locb1)
            if locb1(i)~=0 && i==1 && locb1(length(locb1))~=0
                vaa(2,1)=arroundvertices(i,1);
                j=1;
            elseif locb1(i)~=0
                vaa(j,1)=arroundvertices(i,1);
                j=j+1;
            end
        end
        v111=vaa(1,1);
        v211=vaa(2,1);
        
        tf1222=coord(v111,:)-coord(no,:);
        tf2222=coord(v211,:)-coord(no,:);
        
        xm111=coord(no,:)+epsilon*tf1222;
        xm211=coord(no,:)+epsilon*tf2222;
        
        
        normal22=(R*(xm211-xm111)');
        
        sielemsigma=norm(cross(xm111-coord(no,:),xm211-coord(no,:)))*0.5;
        
        Kigma(1,1)=kmap(elem(ielemk,5),2);
        Kigma(1,2)=kmap(elem(ielemk,5),3);
        Kigma(2,1)=kmap(elem(ielemk,5),4);
        Kigma(2,2)=kmap(elem(ielemk,5),5);
        
        % para elemento posterior =========================================
        
        vertices1=elem(ielemkk,find(elem(ielemkk,1:4)~=0));
        
        [lia1,locb1]=ismember(arroundvertices,vertices1);
        j=1;
        va=zeros(2,1);
        for i=1:length(locb1)
            if locb1(i)~=0 && i==1 && locb1(length(locb1))~=0
                va(2,1)=arroundvertices(i,1);
                j=1;
            elseif locb1(i)~=0
                va(j,1)=arroundvertices(i,1);
                j=j+1;
            end
        end
        v11=va(1,1);
        v21=va(2,1);
        
        tf12=coord(v11,:)-coord(no,:);
        tf22=coord(v21,:)-coord(no,:);
        
        xm11=coord(no,:)+epsilon*tf12;
        xm21=coord(no,:)+epsilon*tf22;
        
        normal2=(R*(xm21-xm11)');
        
        sielemtau=norm(cross(xm11-coord(no,:),xm21-coord(no,:)))*0.5;
        
        Ktau(1,1)=kmap(elem(ielemkk,5),2);
        Ktau(1,2)=kmap(elem(ielemkk,5),3);
        Ktau(2,1)=kmap(elem(ielemkk,5),4);
        Ktau(2,2)=kmap(elem(ielemkk,5),5);
        
        %==============================================================
        
        netas(k,1)=dot(Kigma*normal22,(normal22+(R*(vec1)')))/(2*sielemsigma);
        
        netas(k,2)=dot(K*normal1,(normal1+(R*(-vec1)')))/(2*sielem);
        
        netas(k,3)=dot(K*normal1,(normal1+(R*vec2')))/(2*sielem);
        
        netas(k,4)=dot(Ktau*normal2,(normal2+(R*(-vec2)')))/(2*sielemtau);
        
        
    end
else %Se for um nó do contorno.
    for k=1:nec,
        ielem= esurn1(esurn2(no)+k); % elemento em questão
        
        if (k==1)&&(size(P,1)~=size(O,1))&& k~=nec
            
            ielemkk=esurn1(esurn2(no)+k+1); % elemento posterior
            
            % para ielem
            vertices=elem(ielem,find(elem(ielem,1:4)~=0));
            
            arroundvertices=nsurn1(nsurn2(no) + 1:nsurn2(no + 1));
            
            [lia1,locb1]=ismember( arroundvertices,vertices);
            
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
            
            tf1=coord(v1,:)-coord(no,:);
            tf2=coord(v2,:)-coord(no,:);
            
            xm1=0.5*(coord(no,:)+coord(v1,:));%epsilon*tf1;
            xm2=0.5*(coord(no,:)+coord(v2,:));%epsilon*tf2;
            
            vec1=xm1-coord(no,:);
            vec2=xm2-coord(no,:);
            
            normal1=(R*(xm2-xm1)')/norm(xm2-xm1);
            
            normsigma= norm(coord(no,:)-coord(v1,:));
            normtau= norm(coord(no,:)-coord(v2,:));
            
            K(1,1)=kmap(elem(ielem,5),2);
            K(1,2)=kmap(elem(ielem,5),3);
            K(2,1)=kmap(elem(ielem,5),4);
            K(2,2)=kmap(elem(ielem,5),5);
            
            sielem=norm(cross(xm1-coord(no,:),xm2-coord(no,:)))*0.5;
            
            % para elemento posterior =====================================
            
            vertices1=elem(ielemkk,find(elem(ielemkk,1:4)~=0));
            
            [lia1,locb1]=ismember(arroundvertices,vertices1);
            
            j=1;
            va=zeros(2,1);
            for i=1:length(locb1)
                if locb1(i)~=0 && i==1 && locb1(length(locb1))~=0
                    va(2,1)=arroundvertices(i,1);
                    j=1;
                elseif locb1(i)~=0
                    va(j,1)=arroundvertices(i,1);
                    j=j+1;
                end
            end
            v11=va(1,1);
            v21=va(2,1);
            
            tf1aux=coord(v11,:)-coord(no,:);
            tf2aux=coord(v21,:)-coord(no,:);
            
            %xm11=coord(no,:)+epsilon*tf1aux;
            %xm21=coord(no,:)+epsilon*tf2aux;
            xm11=0.5*(coord(no,:)+coord(v11,:));%epsilon*tf1;
            xm21=0.5*(coord(no,:)+coord(v21,:));%epsilon*tf2;
            normal2=(R*(xm21-xm11)')/norm(xm21-xm11);
            
            sielemtau=norm(cross(xm11-coord(no,:),xm21-coord(no,:)))*0.5;
            
            Ktau(1,1)=kmap(elem(ielemkk,5),2);
            Ktau(1,2)=kmap(elem(ielemkk,5),3);
            Ktau(2,1)=kmap(elem(ielemkk,5),4);
            Ktau(2,2)=kmap(elem(ielemkk,5),5);
            
            %===============================================================
            netas(k,1)=0;
            
            netas(k,2)=dot(K*normal1,(normal1+(R*-vec1')))/(2*sielem);
            
            netas(k,3)=dot(K*normal1,(normal1+(R*vec2')))/(2*sielem);
            
            netas(k,4)=dot(Ktau*normal2,normal2+(R*(-vec2)'))/(2*sielemtau);
            
        elseif (k==1)&&(size(P,1)~=size(O,1))&& k==nec
            
            % elemento único
            vertices=elem(ielem,find(elem(ielem,1:4)~=0)); % para ielem
            
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
            
            tf1=coord(v1,:)-coord(no,:);
            tf2=coord(v2,:)-coord(no,:);
            %xm1=coord(no,:)+epsilon*tf1;
            %xm2=coord(no,:)+epsilon*tf2;
            
            xm1=0.5*(coord(no,:)+coord(v1,:));%epsilon*tf1;
            xm2=0.5*(coord(no,:)+coord(v2,:));%epsilon*tf2;
            
            vec1=xm1-coord(no,:);
            vec2=xm2-coord(no,:);
            
            normal1=(R*(xm2-xm1)')/norm(xm2-xm1);
            
            normsigma= norm(coord(no,:)-coord(v1,:));
            normtau= norm(coord(no,:)-coord(v2,:));
            
            K(1,1)=kmap(elem(ielem,5),2);
            K(1,2)=kmap(elem(ielem,5),3);
            K(2,1)=kmap(elem(ielem,5),4);
            K(2,2)=kmap(elem(ielem,5),5);
            
            sielem=norm(cross(xm1-coord(no,:),xm2-coord(no,:)))*0.5;
            %==============================================================
            netas(k,1)=0;
            
            netas(k,2)=dot(K*normal1,(normal1+(R*-vec1')))/(2*sielem);
            
            netas(k,3)=dot(K*normal1,(normal1+(R*vec2')))/(2*sielem);
            
            netas(k,4)=0;
            
        elseif (k==nec)&&(size(P,1)~=size(O,1))
            
            ielemk=esurn1(esurn2(no)+k-1); % elemento anterior
            
            
            % para ielem
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
            
            tf1=coord(v1,:)-coord(no,:);
            tf2=coord(v2,:)-coord(no,:);
            %xm1=coord(no,:)+epsilon*tf1;
            %xm2=coord(no,:)+epsilon*tf2;
            
            xm1=0.5*(coord(no,:)+coord(v1,:));%epsilon*tf1;
            xm2=0.5*(coord(no,:)+coord(v2,:));%epsilon*tf2;
            
            vec1=xm1-coord(no,:);
            vec2=xm2-coord(no,:);
            
            normal1=(R*(xm2-xm1)')/norm(xm2-xm1);
            
            normsigma= norm(coord(no,:)-coord(v1,:));
            normtau= norm(coord(no,:)-coord(v2,:));
            
            K(1,1)=kmap(elem(ielem,5),2);
            K(1,2)=kmap(elem(ielem,5),3);
            K(2,1)=kmap(elem(ielem,5),4);
            K(2,2)=kmap(elem(ielem,5),5);
            
            Kigma(1,1)=kmap(elem(ielemk,5),2);
            Kigma(1,2)=kmap(elem(ielemk,5),3);
            Kigma(2,1)=kmap(elem(ielemk,5),4);
            Kigma(2,2)=kmap(elem(ielemk,5),5);
            
            sielem=norm(cross(xm1-coord(no,:),xm2-coord(no,:)))*0.5;
            
            % para elemento anterior ======================================
            
            vertices1=elem(ielemk,find(elem(ielemk,1:4)~=0));
            
            [lia1,locb1]=ismember(arroundvertices,vertices1);
            
            j=1;
            va=zeros(2,1);
            for i=1:length(locb1)
                if locb1(i)~=0 && i==1 && locb1(length(locb1))~=0
                    va(2,1)=arroundvertices(i,1);
                    j=1;
                elseif locb1(i)~=0
                    va(j,1)=arroundvertices(i,1);
                    j=j+1;
                end
            end
            v11=va(1,1);
            v21=va(2,1);
            
            
            tf1maux=coord(v11,:)-coord(no,:);
            tf2maux=coord(v21,:)-coord(no,:);
            %xm11=coord(no,:)+epsilon*tf1maux;
            %xm21=coord(no,:)+epsilon*tf2maux;
            
            xm11=0.5*(coord(no,:)+coord(v11,:));%epsilon*tf1;
            xm21=0.5*(coord(no,:)+coord(v21,:));%epsilon*tf2;
            
            normal2=(R*(xm21-xm11)')/norm(xm21-xm11);
            
            
            sielemsigma=norm(cross(xm11-coord(no,:),xm21-coord(no,:)))*0.5;
            %===============================================================
            netas(k,1)=dot(Kigma*normal2,(normal2+(R*(vec1)')))/(2*sielemsigma);
            
            netas(k,2)=dot(K*normal1,(normal1+(R*-vec1')))/(2*sielem);
            
            netas(k,3)=dot(K*normal1,(normal1+(R*vec2')))/(2*sielem);
            
            netas(k,4)=0;
            
        else
            
            ielemk=esurn1(esurn2(no)+k-1); % elemento anterior
            
            ielemkk=esurn1(esurn2(no)+k+1); % elemento posterior
            
            % para ielem
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
            tf1=coord(v1,:)-coord(no,:);
            tf2=coord(v2,:)-coord(no,:);
            %xm1=coord(no,:)+epsilon*tf1;
            %xm2=coord(no,:)+epsilon*tf2;
            
            xm1=0.5*(coord(no,:)+coord(v1,:));%epsilon*tf1;
            xm2=0.5*(coord(no,:)+coord(v2,:));%epsilon*tf2;
            vec1=xm1-coord(no,:);
            vec2=xm2-coord(no,:);
            
            normal1=(R*(xm2-xm1)')/norm(xm2-xm1);
            
            normsigma= norm(coord(no,:)-coord(v1,:));
            normtau= norm(coord(no,:)-coord(v2,:));
            
            K(1,1)=kmap(elem(ielem,5),2);
            K(1,2)=kmap(elem(ielem,5),3);
            K(2,1)=kmap(elem(ielem,5),4);
            K(2,2)=kmap(elem(ielem,5),5);
            
            sielem=norm(cross(xm1-coord(no,:),xm2-coord(no,:)))*0.5;
            
            % para ielemento anterior =====================================
            
            vertices1a=elem(ielemk,find(elem(ielemk,1:4)~=0));
            
            [lia1,locb1]=ismember(arroundvertices, vertices1a );
            
            j=1;
            vaa=zeros(2,1);
            for i=1:length(locb1)
                if locb1(i)~=0 && i==1 && locb1(length(locb1))~=0
                    vaa(2,1)=arroundvertices(i,1);
                    j=1;
                elseif locb1(i)~=0
                    vaa(j,1)=arroundvertices(i,1);
                    j=j+1;
                end
            end
            
            v111=vaa(1,1);
            v211=vaa(2,1);
            
            tf1baux=coord(v111,:)-coord(no,:);
            tf2baux=coord(v211,:)-coord(no,:);
            %xm111=coord(no,:)+epsilon*tf1baux;
            %xm211=coord(no,:)+epsilon*tf2baux;
            
            xm111=0.5*(coord(no,:)+coord(v111,:));%epsilon*tf1;
            xm211=0.5*(coord(no,:)+coord(v211,:));%epsilon*tf2;
            normal22=(R*(xm211-xm111)')/norm(xm211-xm111);
            
            
            sielemsigma=norm(cross(xm111-coord(no,:),xm211-coord(no,:)))*0.5;
            
            Kigma(1,1)=kmap(elem(ielemk,5),2);
            Kigma(1,2)=kmap(elem(ielemk,5),3);
            Kigma(2,1)=kmap(elem(ielemk,5),4);
            Kigma(2,2)=kmap(elem(ielemk,5),5);
            % para elemento posterior ====================================
            
            vertices1=elem(ielemkk,find(elem(ielemkk,1:4)~=0));
            
            [lia1,locb1]=ismember(arroundvertices, vertices1);
            
            j=1;
            va=zeros(2,1);
            for i=1:length(locb1)
                if locb1(i)~=0 && i==1 && locb1(length(locb1))~=0
                    va(2,1)=arroundvertices(i,1);
                    j=1;
                elseif locb1(i)~=0
                    va(j,1)=arroundvertices(i,1);
                    j=j+1;
                end
            end
            
            v11=va(1,1);
            v21=va(2,1);
            
            tf1caux=coord(v11,:)-coord(no,:);
            tf2caux=coord(v21,:)-coord(no,:);
            %xm11=coord(no,:)+epsilon*tf1caux;
            %xm21=coord(no,:)+epsilon*tf2caux;
            
            xm11=0.5*(coord(no,:)+coord(v11,:));%epsilon*tf1;
            xm21=0.5*(coord(no,:)+coord(v21,:));%epsilon*tf2;
            
            normal2=(R*(xm21-xm11)')/norm(xm21-xm11);
            
            sielemtau=norm(cross(xm11-coord(no,:),xm21-coord(no,:)))*0.5;
            
            Ktau(1,1)=kmap(elem(ielemkk,5),2);
            Ktau(1,2)=kmap(elem(ielemkk,5),3);
            Ktau(2,1)=kmap(elem(ielemkk,5),4);
            Ktau(2,2)=kmap(elem(ielemkk,5),5);
            %==============================================================
            netas(k,1)=dot(Kigma*normal22,(normal22+(R*(vec1)')))/(2*sielemsigma);
            
            netas(k,2)=dot(K*normal1,(normal1+(R*-vec1')))/(2*sielem);
            
            netas(k,3)=dot(K*normal1,(normal1+(R*vec2')))/(2*sielem);
            
            netas(k,4)=dot(Ktau*normal2,(normal2+(R*(-vec2)')))/(2*sielemtau);
            
        end
        %% ================================================================
        % calculo do E
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
        % ================================================================
    end
end


end