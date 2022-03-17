function [pressure,errorelativo,flowrate,flowresult,tabletol,coercividade]=solverpressure(kmap,nflagface,nflagno,fonte,...
    tol, nit,p_old,mobility,gamma,wells,parameter,metodoP,...
    auxflag,interpol,Hesq, Kde, Kn, Kt, Ded,weightDMP,auxface,benchmark,...
    iteration,nflag,calnormface,N,gravresult,gravrate)

errorelativo=0;
tabletol=0;
coercividade=0;
% calculo dos pesos que correspondem aos metodos de interpolacao
if strcmp(interpol,'LPEW1')
    [w,s] = Pre_LPEW_1(kmap,N);
elseif strcmp(interpol,'eLPEW2')
    [w,s] = Pre_ELPEW_2(kmap,N,gravrate);
else
    [w,s] = Pre_LPEW_2(kmap,N,gravrate);
end
switch metodoP
    
    case {'nlfvLPEW', 'nlfvDMPSY','nlfvDMPV1','nlfvHP', 'nlfvPPS'}
        
        if strcmp(iteration,'AA')
            
            [pinterp]=pressureinterp(p_old,nflagface,nflagno,w,s,auxflag,metodoP,parameter,weightDMP,mobility);
            % calculo da matriz globlal inicial
            [M_old,RHS_old]=globalmatrix(p_old,pinterp,gamma,nflagface,nflagno,...
                parameter,kmap,fonte,metodoP,w,s,benchmark,weightDMP,auxface,wells,...
                mobility,Hesq, Kde, Kn, Kt, Ded,calnormface,gravresult);
            tic
            [pressure,tabletol,iter,ciclos]=picardAA(M_old,RHS_old,nit,tol,kmap,...
                parameter,metodoP,auxflag,w,s,nflagface,fonte,p_old,gamma,...
                nflagno,benchmark,weightDMP,auxface,wells,mobility,Hesq, ...
                Kde, Kn, Kt, Ded,calnormface);
            toc
        elseif strcmp(iteration,'fullpicard')
            
            % interpolação nos nós ou faces
            [pinterp]=pressureinterp(p_old,nflagface,nflagno,w,s,auxflag,metodoP,parameter,weightDMP,mobility);
            % calculo da matriz globlal inicial
            [M_old,RHS_old]=globalmatrix(p_old,pinterp,gamma,nflagface,nflagno,...
                parameter,kmap,fonte,metodoP,w,s,benchmark,weightDMP,auxface,wells,...
                mobility,Hesq, Kde, Kn, Kt, Ded,calnormface,gravresult);
            tic
            [pressure,tabletol,iter,ciclos]=fullpicard(M_old,RHS_old,nit,tol,kmap,...
                parameter,metodoP,auxflag,w,s,nflagface,fonte,p_old,gamma,...
                nflagno,benchmark,weightDMP,auxface,wells,mobility,Hesq, ...
                Kde, Kn, Kt, Ded,calnormface,gravresult);
            toc
         elseif strcmp(iteration,'iterbroyden')
            [pinterp]=pressureinterp(p_old,nflagface,nflagno,w,s,auxflag,metodoP,parameter,weightDMP,mobility);
            % calculo da matriz globlal inicial
            [M_old,RHS_old]=globalmatrix(p_old,pinterp,gamma,nflagface,nflagno,...
                parameter,kmap,fonte,metodoP,w,s,benchmark,weightDMP,auxface,wells,...
                mobility,Hesq, Kde, Kn, Kt, Ded,calnormface,0);
            p_old1=M_old\RHS_old;
            
            % solver de pressão pelo método Broyden
            %             [pressure,step,errorelativo,flowrate,flowresult]=iterbroyden(M_old,RHS_old,p_old,tol,kmap,parameter,...
            %                 metodoP,auxflag,w,s,nflagface,fonte,gamma,nflagno,benchmark,R0_old,p_old1,...
            %                 weightDMP,auxface,wells,mobility,Hesq, Kde, Kn, Kt, Ded);
            
            [pressure, iter,ciclos,tolerancia]=broyden(M_old,RHS_old,p_old,tol,kmap,parameter,...
                metodoP,auxflag,w,s,nflagface,fonte,gamma,nflagno,benchmark,p_old1,...
                weightDMP,auxface,calnormface,wells,mobility,gravresult);
        elseif strcmp(iteration,'RRE')
            %             [M_old,RHS_old]=assemblematrixlfvLPEW(parameter,fonte,w,s,nflagno,weightDMP,wells,mobility);
            %             p_old=M_old\RHS_old;
            %             p_old=p_old-min(p_old,0);
            % interpolação nos nós ou faces
            [pinterp]=pressureinterp(p_old,nflagface,nflagno,w,s,auxflag,metodoP,parameter,weightDMP,mobility);
            % calculo da matriz globlal inicial
            [M_old,RHS_old]=globalmatrix(p_old,pinterp,gamma,nflagface,nflagno,...
                parameter,kmap,fonte,metodoP,w,s,benchmark,weightDMP,auxface,wells,...
                mobility,Hesq, Kde, Kn, Kt, Ded,calnormface);
            tic
            [pressure,tabletol,iter,ciclos]=picardRRE(M_old,RHS_old,nit,tol,kmap,...
                parameter,metodoP,auxflag,w,s,nflagface,fonte,p_old,gamma,...
                nflagno,benchmark,weightDMP,auxface,wells,mobility,Hesq, ...
                Kde, Kn, Kt, Ded,calnormface);
            toc
        elseif strcmp(iteration,'MPE')
            %             [M_old,RHS_old]=assemblematrixlfvLPEW(parameter,fonte,w,s,nflagno,weightDMP,wells,mobility);
            %             p_old=M_old\RHS_old;
            %             p_old=p_old-min(p_old,0);
            % interpolação nos nós ou faces
            [pinterp]=pressureinterp(p_old,nflagface,nflagno,w,s,auxflag,metodoP,parameter,weightDMP,mobility);
            % calculo da matriz globlal inicial
            [M_old,RHS_old]=globalmatrix(p_old,pinterp,gamma,nflagface,nflagno,...
                parameter,kmap,fonte,metodoP,w,s,benchmark,weightDMP,auxface,wells,...
                mobility,Hesq, Kde, Kn, Kt, Ded,calnormface);
            tic
            [pressure,tabletol,iter,ciclos]=picardMPE(M_old,RHS_old,nit,tol,kmap,...
                parameter,metodoP,auxflag,w,s,nflagface,fonte,p_old,gamma,...
                nflagno,benchmark,weightDMP,auxface,wells,mobility,Hesq, ...
                Kde, Kn, Kt, Ded,calnormface);
            toc
        elseif strcmp(iteration, 'iterdiscretnewton')
            
            p_old1=M_old\RHS_old;
            % interpolação nos nós ou faces
            [pinterp1]=pressureinterp(p_old1,nflagface,nflagno,w,s,auxflag,metodoP,parameter,weightDMP,mobility);
            % calculo da matriz globlal inicial
            [M_old1,RHS_old1]=globalmatrix(p_old1,pinterp1,gamma,nflagface,nflagno,...
                parameter,kmap,fonte,metodoP,w,s,benchmark,weightDMP,auxface,wells,...
                mobility,Hesq, Kde, Kn, Kt, Ded,calnormface);
            
            % resolvedor de pressão pelo método de Newton-Discreto
            [pressure,step,errorelativo,flowrate,flowresult]=iterdiscretnewton(M_old1,RHS_old1,M_old,RHS_old,nit,tol,kmap,...
                parameter,metodoP,auxflag,w,s,nflagface,fonte,p_old,gamma,nflagno,benchmark,...
                weightDMP,auxface,wells,mobility,Hesq, Kde, Kn, Kt, Ded,calnormface,p_old1);
            
            
        elseif strcmp(iteration, 'iterhybrid')
            
            p_old1=M_old\RHS_old;
            
            % interpolação nos nós ou faces
            [pinterp1]=pressureinterp(p_old1,nflagface,w,s,auxflag,metodoP,parameter,weightDMP,mobility);
            
            % calculo da matriz globlal inicial
            [M_old1,RHS_old1]=globalmatrix(p_old1,pinterp1,gamma,nflagface,nflagno,...
                parameter,kmap,fonte,metodoP,w,s,benchmark,weightDMP,auxface,wells,mobility,Hesq, Kde, Kn, Kt, Ded);
            
            % solver pressure pelo método hybrido
            [pressure,step,errorelativo,flowrate,flowresult]=iterhybrid(M_old1,RHS_old1,tol,kmap,...
                parameter,metodoP,auxflag,w,s,nflagface,fonte,p_old,gamma,...
                nflagno,benchmark,p_old1,weightDMP,auxface,wells,mobility,Hesq, Kde, Kn, Kt, Ded);
            
        elseif strcmp(iteration, 'JFNK')
            
            [pinterp]=pressureinterp(p_old,nflagface,nflagno,w,s,auxflag,metodoP,parameter,weightDMP,mobility);
            % calculo da matriz globlal inicial
            [M_old,RHS_old]=globalmatrix(p_old,pinterp,gamma,nflagface,nflagno,...
                parameter,kmap,fonte,metodoP,w,s,benchmark,weightDMP,auxface,wells,...
                mobility,Hesq, Kde, Kn, Kt, Ded,calnormface);
            
            p_old1=M_old\RHS_old;
            % calculo do residuo
            R0=norm(M_old*p_old-RHS_old);
            
            % interpolação nos nós ou faces
            [pinterp1]=pressureinterp(p_old1,nflagface,nflagno,w,s,auxflag,metodoP,parameter,weightDMP,mobility);
            
            % calculo da matriz globlal inicial
            [M_old1,RHS_old1]=globalmatrix(p_old1,pinterp1,gamma,nflagface,nflagno,...
                parameter,kmap,fonte,metodoP,w,s,benchmark,weightDMP,...
                auxface,wells,mobility,Hesq, Kde, Kn, Kt, Ded,calnormface);
            
            % calculo da pressão
            [pressure,iter,ciclos,tolerancia]= JFNK1(tol,kmap,parameter,metodoP,auxflag,w,s,nflagface,fonte,gamma,...
                nflagno,benchmark,M_old1,RHS_old1,p_old1,R0,weightDMP,...
                auxface,wells,mobility,Hesq, Kde, Kn, Kt, Ded,calnormface);
            
        end
        
        % fornecimento das informacoes
        niteracoes=iter*ciclos;
        
        name = metodoP;
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
        
        % interpolação nos nós ou faces
        [pinterp]=pressureinterp(p_old,nflagface,nflagno,w,s,auxflag,...
            metodoP,parameter,weightDMP,mobility,gravresult,gravrate);
        
        % calculo da matriz globlal inicial
        [M_old,RHS_old]=globalmatrix(p_old,pinterp,gamma,nflagface,nflagno,...
            parameter,kmap,fonte,metodoP,w,s,benchmark,weightDMP,auxface,wells,...
            mobility,Hesq, Kde, Kn, Kt, Ded,calnormface,gravresult,gravrate);
        pressure=M_old\RHS_old;
        
        tabletol=0;
        iter=[1,1];
        ciclos=1;
end



if strcmp(metodoP,'nlfvDMPSY')
    pinterp=pressureinterp(pressure,nflagface,nflagno,w,s,auxflag,metodoP,parameter,weightDMP);
    
    %implementação do fluxo NLFV-DMP
    [flowrate,flowresult]=flowrateNLFVDMP(pressure, pinterp, parameter,...
        nflagface,kmap,gamma,weightDMP,mobility);
 
    
elseif strcmp(metodoP,'nlfvHP')
    pinterp=pressureinterp(pressure,nflagface,nflagno,w,s,auxflag,metodoP,parameter,weightDMP);
    
    [flowrate,flowresult]=flowrateNLFVHP(pressure, pinterp, parameter);
    coercividade=0;
elseif strcmp(metodoP,'nlfvLPEW')
    pinterp=pressureinterp(pressure,nflagface,nflagno,w,s,auxflag,metodoP,parameter,weightDMP);
    
    %implementação do fluxo NLFV
    [flowrate,flowresult,coercividade]=flowrateNLFV(pressure, pinterp, parameter,mobility);
    
elseif strcmp(metodoP, 'nlfvPPS')
    
    pinterp=pressureinterp(pressure,nflagface,nflagno,w,s,auxflag,metodoP,parameter,weightDMP);
    
    %implementação do fluxo NLFV
    [flowrate,flowresult,coercividade]=flowrateNLFVPP(pressure, pinterp, parameter,mobility);
    
elseif strcmp(metodoP, 'lfvHP')
    %interpolação nos nós ou faces
    [pinterp]=pressureinterp(pressure,nflagface,nflagno,w,s,auxflag,metodoP,parameter,weightDMP,mobility,g);
    %calculo das vazões
    [flowrate,flowresult]=flowratelfvHP(parameter,weightDMP,mobility,pinterp,pressure,g);
elseif strcmp(metodoP, 'lfvLPEW')
    [pinterp]=pressureinterp(pressure,nflagface,nflagno,w,s,auxflag,metodoP,parameter,weightDMP,mobility);
    %calculo das vazões
    [flowrate,flowresult]=flowratelfvLPEW(parameter,weightDMP,mobility,pinterp,pressure);
    
elseif  strcmp(metodoP, 'tpfa')
    [flowrate, flowresult]=flowrateTPFA(pressure,Kde,Kn,Hesq,nflag,mobility);
else
    %calculo das vazões
    [flowrate,flowresult]=calflowrateMPFAD(pressure,w,s,Kde,Ded,Kn,Kt,Hesq,nflag,auxflag,1,gravresult,gravrate);
end


end