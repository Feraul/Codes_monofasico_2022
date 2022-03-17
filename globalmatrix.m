function [M,I]=globalmatrix(p,pinterp,gamma,nflagface,nflagno,parameter,kmap,...
    fonte,metodoP,w,s,benchmark,weightDMP,auxface,wells,mobility,Hesq, Kde, Kn, Kt, Ded,calnormface,gravresult,gravrate)


if strcmp(metodoP,'nlfvDMP1')
    
    [M,I]=assemblematrixDMP(pinterp,gamma,nflagface,parameter,kmap,fonte);
elseif strcmp(metodoP,'nlfvDMP2')
    
    [M,I]=assemblematrixDMPv1(pinterp,gamma,nflagface,parameter,kmap,fonte);
elseif strcmp(metodoP,'nlfvLPS') || strcmp(metodoP,'nlfvPPS')
    
    [M,I]=assemblematrixLPSPPS(p,pinterp,parameter,fonte);
elseif strcmp(metodoP,'nlfvLPEW')
    
    [M,I]=assemblematrixGYZS(pinterp,parameter,fonte,wells,mobility,calnormface);
elseif strcmp(metodoP,'nlfvHP')
    [M,I]=assemblematrixNLFVHP(pinterp,parameter,fonte,wells,Hesq,Kn,Kt,nflagno);
elseif strcmp(metodoP,'nlfvDMPSY')
    
    [M,I]=assemblematrixDMPSY(p,pinterp,gamma,nflagface,parameter,kmap,fonte,...
        benchmark,weightDMP,auxface,wells,mobility);
elseif strcmp(metodoP,'nlfvDMPV1')
    
    [M,I]=assemblematrixNLFVDMP(p,pinterp,gamma,nflagface,parameter,kmap,...
    fonte,benchmark,weightDMP,auxface,wells,mobility);
elseif strcmp(metodoP,'lfvLPEW')
    
    [M,I]=assemblematrixlfvLPEW(parameter,fonte,w,s,nflagno,weightDMP,wells,mobility);
elseif strcmp(metodoP,'lfvHP')
    
    [M,I]=assemblematrixlfvHPv3(parameter,fonte,nflagface,weightDMP,wells,mobility,g);
elseif strcmp(metodoP,'mpfad')
    [ M, I ] = globalmatrixmpfad( w,s, Kde, Ded, Kn, Kt, nflagno, Hesq,wells,1,fonte,gravresult,gravrate);
elseif strcmp(metodoP,'tpfa')
    [ M, I ] = globalmatrixtpfa( Kde, Kn, nflag, Hesq,wells,mobility);
   
end
end