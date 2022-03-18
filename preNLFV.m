function [pointarmonic,parameter,gamma,p_old,tol,nit,er,nflagface,...
    nflagno,weightDMP,Hesq,Kde,Kn,Kt,Ded,auxface,calnormface,...
    gravresult,gravrate,weight,contrcontor]=preNLFV(kmap,...
    N,metodoP,benchmark,bedge,grav,gravitational,correction,interpol)
global elem
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
if strcmp(metodoP,'nlfvLPEW')
    %% calculo dos parametros ou constantes (ksi)
    % esta rutina estamos usando de 7/2/2016
    %[parameter]=coefficientPPS(kmap); % urgente revisar
    %temos usado para muitos estes o seguinte rutina
    [parameter,calnormface]=coefficientLPSangle(kmap);
    % adequação dos flags de contorno
    nflagno= contflagno(benchmark,bedge);
    
    % calculo dos pesos que correspondem aos metodos de interpolacao
    if strcmp(interpol,'LPEW1')
        % interpolaca LPEW1 proposto por Gao e Wu 2010
        [weight,contrcontor] = Pre_LPEW_1(kmap,N);
    elseif strcmp(interpol,'eLPEW2')
        % interpolaca LPEW2 modificado por proposto por Miao e Wu 2021
        [weight,contrcontor] = Pre_ELPEW_2(kmap,N,gravrate);
    else
        % interpolaca LPEW1 proposto por Gao e Wu 2010
        [weight,contrcontor] = Pre_LPEW_2(kmap,N,gravrate);
    end
    % calculo do termo gravitacional
    if strcmp(gravitational,'yes')
        [gravresult,gravrate]=gravitation(kmap,grav);
    end
elseif strcmp(metodoP,'nlfvPPS')
    %% calculo dos parametros ou constantes (ksi)
    % esta rutina estamos usando de 7/2/2016
    %[parameter]=coefficientPPS(kmap); % urgente revisar
    %temos usado para muitos estes o seguinte rutina
    [parameter,calnormface]=coefficientLPSangle(kmap);
    % adequação dos flags de contorno
    nflagno= contflagno(benchmark,bedge);
    
    if strcmp(correction,'yes')
        % calculo dos pontos harmonicos com correcao
        [pointarmonic,weightDMP,raioaux]=harmonicopointcorrection(kmap);
    else
        % calculoa dos pontos harmonicos sem correcao
        [pointarmonic,weightDMP,raioaux]=harmonicopoint(kmap,N,benchmark);
    end
    
    % calculo do termo gravitacional
    if strcmp(gravitational,'yes')
        [gravresult,gravrate]=gravitation(kmap,grav);
    end
elseif strcmp(metodoP,'nlfvHP')
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
    % adequação dos flag de face de contorno
    nflagface= contflagface(benchmark,bedge);
    % adequação dos nos flags de contorno
    nflagno= contflagno(benchmark,bedge);
    %calculo de parametros
    [Hesq, Kde, Kn, Kt, Ded]=preMPFAD(kmap);
    
elseif strcmp(metodoP,'nlfvDMPSY')|| strcmp(metodoP,'lfvHP') || strcmp(metodoP,'nlfvDMPV1')
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
    % adequação dos flag de face de contorno
    nflagface= contflagface(benchmark,bedge);
    % adequação dos nos flags de contorno
    nflagno= contflagno(benchmark,bedge);
    % gravitational term
    if strcmp(gravitational,'yes')
        [gravresult,gravrate]=gravitation(kmap,grav);
    end
elseif strcmp(metodoP,'lfvLPEW')
    %% calculo dos parametros ou constantes (ksi)
    
    %[parameter]=coefficientLPS(kmap);
    [parameter]=coefficientLPSangle(kmap);
    % adequação dos flags de contorno
    nflagno= contflagno(benchmark,bedge);
    % calculo dos pesos DMP
    [weightDMP]=weightnlfvDMP(kmap);
    
    
    % calculo dos pesos que correspondem aos metodos de interpolacao
    if strcmp(interpol,'LPEW1')
        % interpolaca LPEW1 proposto por Gao e Wu 2010
        [weight,contrcontor] = Pre_LPEW_1(kmap,N);
    elseif strcmp(interpol,'eLPEW2')
        % interpolaca LPEW2 modificado por proposto por Miao e Wu 2021
        [weight,contrcontor] = Pre_ELPEW_2(kmap,N,gravrate);
    else
        % interpolaca LPEW1 proposto por Gao e Wu 2010
        [weight,contrcontor] = Pre_LPEW_2(kmap,N,gravrate);
    end
    
    if strcmp(gravitational,'yes')
        [gravresult,gravrate]=gravitation(kmap,grav);
    end
    % outra maneira de calcular os pesos proposto no artigo
    %[weightDMP]=weightlfv(parameter);
elseif strcmp(metodoP,'mpfad')
    
    %% calculo das constantes fisicos-geometrico
    [Hesq, Kde, Kn, Kt, Ded]=preMPFAD(kmap);
    % adequação dos flags de contorno
    nflagno= contflagno(benchmark,bedge);
    
    % calculo dos pesos que correspondem aos metodos de interpolacao
    if strcmp(interpol,'LPEW1')
        % interpolaca LPEW1 proposto por Gao e Wu 2010
        [weight,contrcontor] = Pre_LPEW_1(kmap,N);
    elseif strcmp(interpol,'eLPEW2')
        % interpolaca LPEW2 modificado por proposto por Miao e Wu 2021
        [weight,contrcontor] = Pre_ELPEW_2(kmap,N,gravrate);
    else
        % interpolaca LPEW1 proposto por Gao e Wu 2010
        [weight,contrcontor] = Pre_LPEW_2(kmap,N,gravrate);
    end
    % gravitational term
    if strcmp(gravitational,'yes')
        [gravresult,gravrate]=gravitation(kmap,grav);
    end
else
    %% calculo das constantes fisicos-geometrico para o TPFA
    [Hesq, Kde, Kn, Kt, Ded]=preMPFAD(kmap);
    % adequação dos flags de contorno
    nflagno= contflagno(benchmark);
    
end
%% dados inicialização métodos dos volumes finitos não linear
gamma=0.0;                     % este parametro esta no intervalo [0,1] pode ser utilizado para o método nao linear MPFA
p_old=1*ones(size(elem,1),1);  % inicializando a presao
tol=1e-9;                      % tolerancia para metodos não lineares
nit=2000;                      % numero de iteracoes de Picard
er=1;                          % inicializacao do erro
end