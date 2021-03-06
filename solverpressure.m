function [pressure,errorelativo,flowrate,flowresult,tabletol,coercividade]=...
    solverpressure(kmap,nflagface,nflagno,fonte,...
    tol, nit,p_old,mobility,gamma,wells,parameter,...
    Hesq, Kde, Kn, Kt, Ded,weightDMP,auxface,...
    calnormface,gravresult,gravrate,w,s,gravno,gravelem,gravface)
global iteration pmetodo
errorelativo=0;
tabletol=0;
coercividade=0;

switch pmetodo
    
    case {'nlfvLPEW', 'nlfvDMPSY','nlfvDMPV1','nlfvHP', 'nlfvPPS','nlfvLPS'}
        % intepolacao dos pontos auxiliarea com o antigo campo de pressao
        [pinterp]=pressureinterp(p_old,nflagface,nflagno,w,s,parameter,...
            weightDMP,mobility);
        % calculo da matrizes globlais iniciais
        [M_old,RHS_old]=globalmatrix(p_old,pinterp,gamma,nflagface,nflagno,...
            parameter,kmap,fonte,w,s,weightDMP,auxface,wells,...
            mobility,Hesq, Kde, Kn, Kt, Ded,calnormface,gravresult,gravrate,...
            gravno,gravelem,gravface);
        if strcmp(iteration,'AA')
            
            tic
            [pressure,tabletol,iter,ciclos]=picardAA(M_old,RHS_old,nit,tol,kmap,...
                parameter,w,s,nflagface,fonte,p_old,gamma,...
                nflagno,weightDMP,auxface,wells,mobility,Hesq, ...
                Kde, Kn, Kt, Ded,calnormface,gravresult,gravrate,...
            gravno,gravelem,gravface);
            toc
        elseif strcmp(iteration,'fullpicard')
            
            tic
            [pressure,tabletol,iter,ciclos]=fullpicard(M_old,RHS_old,...
                nit,tol,kmap,...
                parameter,w,s,nflagface,fonte,p_old,gamma,...
                nflagno,weightDMP,auxface,wells,mobility,Hesq, ...
                Kde, Kn, Kt, Ded,calnormface,gravresult,gravrate,...
            gravno,gravelem,gravface);
            toc
        elseif strcmp(iteration,'iterbroyden')
            
            p_old1=M_old\RHS_old;
            
            % solver de press?o pelo m?todo Broyden
            %             [pressure,step,errorelativo,flowrate,flowresult]=iterbroyden(M_old,RHS_old,p_old,tol,kmap,parameter,...
            %                 metodoP,auxflag,w,s,nflagface,fonte,gamma,nflagno,benchmark,R0_old,p_old1,...
            %                 weightDMP,auxface,wells,mobility,Hesq, Kde, Kn, Kt, Ded);
            
            [pressure, iter,ciclos,tolerancia]=broyden(M_old,RHS_old,...
                p_old,tol,kmap,parameter,...
                w,s,nflagface,fonte,gamma,nflagno,p_old1,...
                weightDMP,auxface,calnormface,wells,mobility,gravresult);
        elseif strcmp(iteration,'RRE')
            
            tic
            [pressure,tabletol,iter,ciclos]=picardRRE(M_old,RHS_old,nit,tol,kmap,...
                parameter,w,s,nflagface,fonte,p_old,gamma,...
                nflagno,weightDMP,auxface,wells,mobility,Hesq, ...
                Kde, Kn, Kt, Ded,calnormface);
            toc
        elseif strcmp(iteration,'MPE')
            
            tic
            [pressure,tabletol,iter,ciclos]=picardMPE(M_old,RHS_old,nit,tol,kmap,...
                parameter,w,s,nflagface,fonte,p_old,gamma,...
                nflagno,weightDMP,auxface,wells,mobility,Hesq, ...
                Kde, Kn, Kt, Ded,calnormface);
            toc
        elseif strcmp(iteration, 'iterdiscretnewton')
            
            p_old1=M_old\RHS_old;
            % interpola??o nos n?s ou faces
            [pinterp1]=pressureinterp(p_old1,nflagface,nflagno,w,s,...
                parameter,weightDMP,mobility);
            % calculo da matriz globlal inicial
            [M_old1,RHS_old1]=globalmatrix(p_old1,pinterp1,gamma,nflagface,nflagno,...
                parameter,kmap,fonte,w,s,weightDMP,auxface,wells,...
                mobility,Hesq, Kde, Kn, Kt, Ded,calnormface);
            
            % resolvedor de press?o pelo m?todo de Newton-Discreto
            [pressure,iter,ciclos,tolerancia]=iterdiscretnewton(M_old1,...
                RHS_old1,M_old,RHS_old,nit,tol,kmap,...
                parameter,w,s,nflagface,fonte,p_old,gamma,nflagno,benchmark,...
                weightDMP,auxface,wells,mobility,Hesq, Kde, Kn, Kt,...
                Ded,calnormface,p_old1);
            
            
        elseif strcmp(iteration, 'iterhybrid')
            
            p_old1=M_old\RHS_old;
            
            % interpola??o nos n?s ou faces
            [pinterp1]=pressureinterp(p_old1,nflagface,w,s,parameter,weightDMP,mobility);
            
            % calculo da matriz globlal inicial
            [M_old1,RHS_old1]=globalmatrix(p_old1,pinterp1,gamma,nflagface,nflagno,...
                parameter,kmap,fonte,w,s,weightDMP,auxface,wells,mobility,Hesq, Kde, Kn, Kt, Ded);
            
            % solver pressure pelo m?todo hybrido
            [pressure,iter,ciclos,tolerancia]=iterhybrid(M_old1,RHS_old1,tol,kmap,...
                parameter,w,s,nflagface,fonte,p_old,gamma,...
                nflagno,p_old1,weightDMP,auxface,wells,mobility,Hesq, ...
                Kde, Kn, Kt, Ded);
            
        elseif strcmp(iteration, 'JFNK')
            
            [pinterp]=pressureinterp(p_old,nflagface,nflagno,w,s,...
                parameter,weightDMP,mobility);
            % calculo da matriz globlal inicial
            [M_old,RHS_old]=globalmatrix(p_old,pinterp,gamma,nflagface,nflagno,...
                parameter,kmap,fonte,w,s,weightDMP,auxface,wells,...
                mobility,Hesq, Kde, Kn, Kt, Ded,calnormface);
            
            p_old1=M_old\RHS_old;
            % calculo do residuo
            R0=norm(M_old*p_old-RHS_old);
            
            % interpola??o nos n?s ou faces
            [pinterp1]=pressureinterp(p_old1,nflagface,nflagno,w,s,...
                parameter,weightDMP,mobility);
            
            % calculo da matriz globlal inicial
            [M_old1,RHS_old1]=globalmatrix(p_old1,pinterp1,gamma,...
                nflagface,nflagno,...
                parameter,kmap,fonte,w,s,weightDMP,...
                auxface,wells,mobility,Hesq, Kde, Kn, Kt, Ded,calnormface);
            
            % calculo da press?o
            [pressure,iter,ciclos,tolerancia]= JFNK1(tol,kmap,parameter,w,...
                s,nflagface,fonte,gamma,...
                nflagno,M_old1,RHS_old1,p_old1,R0,weightDMP,...
                auxface,wells,mobility,Hesq, Kde, Kn, Kt, Ded,calnormface);
            
        end
        
        % fornecimento das informacoes
        niteracoes=iter*ciclos;
        
        name = pmetodo;
        X = sprintf('>> calculo o campo de pressao pelo metodo: %s ',name);
        disp(X)
        if strcmp(iteration,'iterbroyden')
            x=['>> tolerancia alcancado:',num2str(tolerancia)];
            disp(x);
        else
            x=['>> tolerancia alcancado:',num2str(tabletol(size(tabletol,1),2))];
            disp(x);
        end
        name1=iteration;
        z = sprintf('>> metodo iterativo utilizado: %s ',name1);
        disp(z)
        y=['>> numero de iteracoes:',num2str(niteracoes)];
        disp(y);
        
        
    case {'lfvHP','lfvLPEW','mpfad','tpfa'}
        % incializando variaveis
        
        pinterp=0;
        % calculo da matriz globlal inicial
        [M_old,RHS_old]=globalmatrix(p_old,pinterp,gamma,nflagface,nflagno,...
            parameter,kmap,fonte,w,s,weightDMP,auxface,wells,...
            mobility,Hesq, Kde, Kn, Kt, Ded,calnormface,gravresult,gravrate,...
            gravno,gravelem,gravface);
        pressure=M_old\RHS_old;
        
        
        
        tabletol=0;
        name = pmetodo;
        X = sprintf('>> calculo o campo de pressao pelo metodo: %s ',name);
        disp(X)
        
end

% intepolacao dos pontos auxiliare a com o novo campo de pressao
pinterp=pressureinterp(pressure,nflagface,nflagno,w,s,parameter,weightDMP);
if strcmp(pmetodo,'nlfvDMPSY')
    %implementa??o do fluxo NLFV-DMP
    [flowrate,flowresult]=flowrateNLFVDMP(pressure, pinterp, parameter,...
        nflagface,kmap,gamma,weightDMP,mobility);
elseif strcmp(pmetodo,'nlfvHP')
    [flowrate,flowresult]=flowrateNLFVHP(pressure, pinterp, parameter);
    coercividade=0;
elseif strcmp(pmetodo,'nlfvLPEW')
    %implementa??o do fluxo NLFV
    [flowrate,flowresult,coercividade]=flowrateNLFV(pressure, pinterp,...
        parameter,mobility);
elseif strcmp(pmetodo, 'nlfvLPS')
    %implementa??o do fluxo NLFV
    [flowrate,flowresult,coercividade]=flowrateNLFVPP(pressure, pinterp,...
        parameter,mobility);
    
elseif strcmp(pmetodo, 'nlfvPPS')
    %implementa??o do fluxo NLFV
    [flowrate,flowresult,coercividade]=flowrateNLFVPP(pressure, pinterp,...
        parameter,mobility);
    
elseif strcmp(pmetodo, 'lfvHP')
    
    %calculo das vaz?es
    [flowrate,flowresult]=flowratelfvHP(parameter,weightDMP,mobility,...
        pinterp,pressure);
elseif strcmp(pmetodo, 'lfvLPEW')
    %calculo das vaz?es
    [flowrate,flowresult]=flowratelfvLPEW(parameter,weightDMP,mobility,...
        pinterp,pressure);
    
elseif  strcmp(pmetodo, 'tpfa')
    [flowrate, flowresult]=flowrateTPFA(pressure,Kde,Kn,Hesq,nflagface,...
        mobility,gravresult,gravrate,pinterp);
else
    %calculo das vaz?es
    [flowrate,flowresult]=calflowrateMPFAD(pressure,w,s,Kde,Ded,Kn,Kt,...
        Hesq,nflagno,1,gravresult,gravrate,pinterp,gravno,gravelem);
end


end