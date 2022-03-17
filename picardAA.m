function [p_new,tabletol,iter,ciclos]=picardAA(M_old,RHS_old,nitpicard,tolpicard,kmap,...
    parameter,metodoP,auxflag,w,s,nflagface,fonte,p_old,gamma,nflagno,benchmark,...
    weightDMP,auxface,wells,mobility,Hesq, Kde, Kn, Kt, Ded,calnormface)

%% calculo do residuo Inicial
R0=norm(M_old*p_old-RHS_old);
ciclos=1;
%% inicializando dados para iteração Picard

% if rcond(full(M_old))<1e-5
    [L,U] = ilu(M_old,struct('type','ilutp','droptol',1e-6));
    
%    [p_old,]=bicgstab(M_old,RHS_old,1e-11,1000,L,U,p_old);
% else
%    [p_old,]=bicgstab(M_old,RHS_old,1e-11,1000);
%end
[p_old,fl1,rr1,it1,rv1]=gmres(M_old,RHS_old,10,1e-9,1000,L,U);
[p_new,iter,res_hist,tabletol]=AndAcc(p_old,1e-8,kmap,parameter,metodoP,auxflag,w,s,...
    nflagface,fonte,gamma,nflagno,benchmark,weightDMP,auxface,...
    wells,mobility,Hesq, Kde, Kn, Kt, Ded,calnormface,R0,tolpicard);
end