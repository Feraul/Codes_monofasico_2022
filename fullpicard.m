function [p_new,tabletol,step,ciclos]=fullpicard(M_old,RHS_old,nitpicard,tolpicard,kmap,...
    parameter,metodoP,auxflag,w,s,nflagface,fonte,p_old,gamma,nflagno,benchmark,...
    weightDMP,auxface,wells,mobility,Hesq, Kde, Kn, Kt, Ded,calnormface,gravresult)
ciclos=1;
%% calculo do residuo Inicial
R0=norm(M_old*p_old-RHS_old);

%% inicializando dados para iteração Picard
step=0;
er=1;
contador=0;
while (tolpicard<er || tolpicard==er) && (step<nitpicard)
    %% atualiza iterações
    step=step+1
    
%    p_new=M_old\RHS_old;  % inversão sem pivotamento
%     if rcond(full(M_old))<1e-5
        [L,U] = ilu(M_old,struct('type','ilutp','droptol',1e-6));
%        [p_new,fl1,rr1,it1,rv1]=bicgstab(M_old,RHS_old,10^-12,150);%,100,L,U,p_old);
%     else
%        [p_new,fl1,rr1,it1,rv1]=bicgstab(M_old,RHS_old,min(10^-5,10^-1*er),1000);
%    end
    [p_new]=gmres(M_old,RHS_old,10,1e-9,1000,L,U);
    %% plotagem no visit
         S=ones(size(p_new,1),1);
         postprocessor(p_new,step)
%         p_max=max(p_new)
%         p_min=min(p_new)
    %% Interpolação das pressões na arestas (faces)
    [pinterp_new]=pressureinterp(p_new,nflagface,nflagno,w,s,auxflag,metodoP,parameter,weightDMP);
    
    %% Calculo da matriz global
    [M_new,RHS_new]=globalmatrix(p_new,pinterp_new,gamma,nflagface,nflagno...
        ,parameter,kmap,fonte,metodoP,w,s,benchmark,weightDMP,auxface,...
        wells,mobility,Hesq, Kde, Kn, Kt, Ded,calnormface,gravresult);
    min(RHS_new)
    % calculo do residuo
    R = norm(M_new*p_new - RHS_new);
    % calculo do erro
    if (R0 ~= 0.0)
        er = abs(R/R0)
    else
        er = 0.0; %exact
    end

    M_old=M_new;
    RHS_old=RHS_new;
    p_old=p_new;
    %if rem(step,5)==0 || step==1
    tabletol(contador+1,1:2)=[contador, er];
    contador=contador+1;
    %end
end
end