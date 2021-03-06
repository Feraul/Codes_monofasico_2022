function [pointarmonic,parameter,gamma,p_old,tol,nit,er,nflagface,...
    nflagno,weightDMP,Hesq,Kde,Kn,Kt,Ded,auxface,calnormface,...
    gravresult,gravrate,weight,contrcontor]=preprocessorlocal(kmap,...
    N,grav,gravface,gravelem)
global elem gravitational pmetodo interpol correction bedge strategy
% inicializando as variaveis
nflagno=0;
nflagface=0;
pointarmonic=0;
parameter=0;
auxface=0;
weightDMP=0;
Hesq=0;
Kde=0;
Kn=0;
Kt=0;
Ded=0;
calnormface=0;
gravresult=0;
gravrate=0;
weight=0;
contrcontor=0;

% gravitational term
if strcmp(gravitational,'yes')
    if strcmp(strategy,'starnoni')
        
        %[gravresult,gravrate]=gravitation(kmap,grav,gravface);
        [gravresult,gravrate]=gravitation_aux(kmap,grav,gravface);
    end
end

if strcmp(pmetodo,'nlfvLPEW')
    %% calculo dos parametros ou constantes (ksi)
    % esta rutina estamos usando de 7/2/2016
    %[parameter]=coefficientPPS(kmap); % urgente revisar
    %temos usado para muitos estes o seguinte rutina
    [parameter,calnormface]=coefficientLPSangle(kmap);
    % adequa??o dos flags de contorno
    nflagno= contflagno;
    
    % calculo dos pesos que correspondem aos metodos de interpolacao
    if strcmp(interpol,'LPEW1')
        % interpolaca LPEW1 proposto por Gao e Wu 2010
        [weight,contrcontor] = Pre_LPEW_1(kmap,N);
    elseif strcmp(interpol,'eLPEW2')
        % interpolaca LPEW2 modificado por proposto por Miao e Wu 2021
        [weight,contrcontor] = Pre_ELPEW_2(kmap,N,gravrate);
    elseif strcmp(interpol,'LS')
        [ weight,contrcontor] = LS(kmap);
    elseif strcmp(interpol,'eLS')
        disp('>> falta implementar!')
    else
        % interpolaca LPEW1 proposto por Gao e Wu 2010
        [weight,contrcontor] = Pre_LPEW_2(kmap,N);
    end
    
elseif strcmp(pmetodo,'nlfvLPS')
    %% calculo dos parametros ou constantes (ksi)
    % esta rutina estamos usando de 7/2/2016
    %[parameter]=coefficientPPS(kmap); % urgente revisar
    %temos usado para muitos estes o seguinte rutina
    [parameter,calnormface]=coefficientLPSangle(kmap);
    % adequa??o dos flags de contorno
    nflagno= contflagno;
    
    % calculo dos pesos que correspondem aos metodos de interpolacao
    if strcmp(interpol,'LPEW1')
        % interpolaca LPEW1 proposto por Gao e Wu 2010
        [weight,contrcontor] = Pre_LPEW_1(kmap,N);
    elseif strcmp(interpol,'eLPEW2')
        % interpolaca LPEW2 modificado por proposto por Miao e Wu 2021
        [weight,contrcontor] = Pre_ELPEW_2(kmap,N,gravrate);
    elseif strcmp(interpol,'LS')
        [ weight,contrcontor] = LS(kmap);
    elseif strcmp(interpol,'eLS')
        disp('>> falta implementar!')
    else
        % interpolaca LPEW1 proposto por Gao e Wu 2010
        [weight,contrcontor] = Pre_LPEW_2(kmap,N,gravrate);
    end
    
elseif strcmp(pmetodo,'interpfree')
    [parameter]=coeffinterpfree(kmap,F);
    
    
elseif strcmp(pmetodo,'nlfvPPS')
    %% calculo dos parametros ou constantes (ksi)
    % esta rutina estamos usando de 7/2/2016
    %[parameter]=coefficientPPS(kmap); % urgente revisar
    %temos usado para muitos estes o seguinte rutina
    [parameter,calnormface]=coefficientLPSangle(kmap);
    % adequa??o dos flags de contorno
    nflagno= contflagno;
    
    if strcmp(correction,'yes')
        % calculo dos pontos harmonicos com correcao
        [pointarmonic,weightDMP,raioaux]=harmonicopointcorrection(kmap);
    else
        % calculoa dos pontos harmonicos sem correcao
        [pointarmonic,weightDMP,raioaux]=harmonicopoint(kmap,N,benchmark);
    end
elseif strcmp(pmetodo,'nlfvHP')
    %% faces alrededor de um elemento
    [facelement]=element_face;
    
    if strcmp(correction,'yes')
        % calculo dos pontos harmonicos com correcao
        [pointarmonic,weightDMP,raioaux]=harmonicopointcorrection(kmap);
    else
        % calculoa dos pontos harmonicos sem correcao
        [pointarmonic,weightDMP,raioaux]=harmonicopoint(kmap,N);
    end
    %% calculo dos parametros ou constantes (ksi)
    % temos usado este parametro durante muito tempo em muitos testes
    [parameter,auxface]=coefficientPPSharmonicpoint(facelement,pointarmonic,kmap,raioaux);
    % adequa??o dos flag de face de contorno
    nflagface= contflagface;
    % adequa??o dos nos flags de contorno
    nflagno= contflagno;
    %calculo de parametros
    [Hesq, Kde, Kn, Kt, Ded]=preMPFAD(kmap);
    
elseif strcmp(pmetodo,'nlfvDMPSY')|| strcmp(pmetodo,'lfvHP') || strcmp(pmetodo,'nlfvDMPV1')
    %% faces alrededor de um elemento
    [facelement]=element_face;
    
    if strcmp(correction,'yes')
        % calculo dos pontos harmonicos com correcao
        [pointarmonic,weightDMP,raioaux]=harmonicopointcorrection(kmap);
    else
        % calculoa dos pontos harmonicos sem correcao
        [pointarmonic,weightDMP,raioaux]=harmonicopoint(kmap,N,benchmark);
    end
    %% calculo dos parametros ou constantes (ksi)
    % temos usado este parametro durante muito tempo em muitos testes
    [parameter,auxface]=coefficientPPSharmonicpoint(facelement,pointarmonic,kmap,raioaux);
    % esta rutina estamos usando de 7/2/2016
    %[parameter]=coefficientPPSusingHP(kmap,facelement,pointarmonic); %para lfvHP
    % adequa??o dos flag de face de contorno
    nflagface= contflagface;
    % adequa??o dos nos flags de contorno
    nflagno= contflagno;
    
elseif strcmp(pmetodo,'lfvLPEW')
    %% calculo dos parametros ou constantes (ksi)
    
    %[parameter]=coefficientLPS(kmap);
    [parameter]=coefficientLPSangle(kmap);
    % adequa??o dos flags de contorno
    nflagno= contflagno;
    % calculo dos pesos DMP
    [weightDMP]=weightnlfvDMP(kmap);
    
    
    % calculo dos pesos que correspondem aos metodos de interpolacao
    if strcmp(interpol,'LPEW1')
        % interpolaca LPEW1 proposto por Gao e Wu 2010
        [weight,contrcontor] = Pre_LPEW_1(kmap,N);
    elseif strcmp(interpol,'eLPEW2')
        % interpolaca LPEW2 modificado por proposto por Miao e Wu 2021
        [weight,contrcontor] = Pre_ELPEW_2(kmap,N,gravrate);
    elseif strcmp(interpol,'LS')
        [ weight,contrcontor] = LS(kmap);
    elseif strcmp(interpol,'eLS')
        disp('>> falta implementar!')
    else
        % interpolaca LPEW1 proposto por Gao e Wu 2010
        [weight,contrcontor] = Pre_LPEW_2(kmap,N,gravrate);
    end
    
    % outra maneira de calcular os pesos proposto no artigo
    %[weightDMP]=weightlfv(parameter);
elseif strcmp(pmetodo,'mpfad')
    
    %% calculo das constantes fisicos-geometrico
    [Hesq, Kde, Kn, Kt, Ded]=preMPFAD(kmap);
    % adequa??o dos flags de contorno
    nflagno= contflagno;
    
    % calculo dos pesos que correspondem aos metodos de interpolacao
    if strcmp(interpol,'LPEW1')
        % interpolaca LPEW1 proposto por Gao e Wu 2010
        [weight,contrcontor] = Pre_LPEW_1(kmap,N);
    elseif strcmp(interpol,'eLPEW2')
        % interpolaca LPEW2 modificado por proposto por Miao e Wu 2021
        [weight,contrcontor] = Pre_ELPEW_2(kmap);
    elseif strcmp(interpol,'LS')
        [ weight,contrcontor] = LS(kmap);
        
    elseif strcmp(interpol,'eLS')
        disp('>> falta implementar!')
        
    else
        % caso contrario utiliza interpolacao LPEW2
        % interpolaca LPEW1 proposto por Gao e Wu 2010
        [weight,contrcontor] = Pre_LPEW_2(kmap,N);
    end
    
else
    % calculo das constantes fisicos-geometrico para o TPFA
    [Hesq, Kde, Kn, Kt, Ded]=preMPFAD(kmap);
    % adequa??o dos flags de contorno
    nflagface= contflagface;
end
%% dados inicializa??o m?todos dos volumes finitos n?o linear
gamma=0.0;                     % este parametro esta no intervalo [0,1] pode ser utilizado para o m?todo nao linear MPFA
p_old=1*ones(size(elem,1),1);  % inicializando a presao
tol=1e-10;                      % tolerancia para metodos n?o lineares
nit=2000;                      % numero de iteracoes de Picard
er=1;                          % inicializacao do erro
end