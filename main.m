% Simulador para resolver a equacao de eliptica (2-D) 
% Desenvolvedor: Prof. Fernando R.L. Contreras

%% Este codigo somente roda MONOFASICO
clear all
clc
format short
global coord centelem elem esurn1 esurn2 nsurn1 nsurn2 bedge inedge ...
    normals esureface1 esureface2 esurefull1 esurefull2 elemarea dens ...
    visc satlimit pormap bcflag courant totaltime filepath foldername;
%%========================================================================%

[coord,centelem,elem,esurn1,esurn2,nsurn1,nsurn2,bedge,inedge,normals,...
    esureface1,esureface2,esurefull1,esurefull2,elemarea,dens,visc,...
    satlimit,pormap,bcflag,courant,totaltime,filepath,foldername,kmap,...
    wells] = preprocessor;
%% NOTAS 
% 1. Para o interpolador com correcao de pontos harmonicos precisa ainda
% implementar o caso artigo Zhang Kobaise figura 12. 
% 2. Ainda falta investir no termo gravitacional
% 3. Verificar o codigo eLPW2

%%
% bedge(217,4)=102;
% bedge(223,4)=102;
% bedge(229,4)=102;
% bedge(235,4)=102;
%% só em malha com furo 1
% x=bedge(217:240,1);
% y=bedge(217:240,2);
% bedge(217:240,1)=y;
% bedge(217:240,2)=x;
%% só em malha com furo 2
% x=bedge(129:144,1);
% y=bedge(129:144,2);
% bedge(129:144,1)=y;
% bedge(129:144,2)=x;
%% distorcao de malhas estruturadas
% esta funcao pode ser ativado se deseja distorcer alguma malha estruturada
%[auxcoord]=distortedramd;

%% somente use nas malhas "Tipo1malha1", "Tipo1malha2", "Tipo1malha3" e "Tipo1malha4"

% Historial para malha "Tipo1malha0" não habilite nada, el ya viene
% ordenado en sentido antihorario tanto en elemento como no contorno
% para malha "Tipo1malha1" "Tipo1malha2 " "Tipo1malha3" e "Tipo1malha4"
% vamos habilitar na linha 803-806 do preprocessador,
% e na linha, caso contrario vai dar erro no calculo do esurn1 esurn2 e no
% nsurn1 e no nsurn2

%   x=bedge(:,1);
%   y=bedge(:,2);
%   bedge(:,1)=y;
%   bedge(:,2)=x;
%   x1=elem(:,1);
%   x2=elem(:,3);
%   elem(:,1)=x2;
%   elem(:,3)=x1;

%% Modificação Malha Kershaw
%bedge(:,4:5)=101;
%----------------------------
%% tratamento malha Hermeline
%bedge(:,4:5)=101;
% malha 16x16
%x=bedge(16:24,1);
%y=bedge(16:24,2);
%bedge(16:24,1)=y;
%bedge(16:24,2)=x;
% malha 32x32
% x=bedge(33:48,1);
% y=bedge(33:48,2);
% bedge(33:48,1)=y;
% bedge(33:48,2)=x;
% malha 64x64
% x=bedge(64:96,1);
% y=bedge(64:96,2);
% bedge(64:96,1)=y;
% bedge(64:96,2)=x;
% malha 128x128
%x=bedge(128:192,1);
%y=bedge(128:192,2);
%bedge(128:192,1)=y;
%bedge(128:192,2)=x;
%% calculo o flag do elemento que deseja
%   a=6287;
%   b=445;
%   c=5740;
%   d=0;
%   [elemento]=searchelement(a,b,c,d)
%% escolha o tipo de erro discreto que deseja usar
% erromethod1 ---> erro utilizado por Gao e Wu 2010
% erromethod2 --->  ''     ''     por Lipnikov et al 2010
% erromethod3 --->  ''     ''     por Eigestad et al 2005
% erromethod4 --->  ''     ''     por Shen e Yuan 2015
erromethod='erromethod1';
%% Defina o tipo de solver 
% tpfa      --> método Linear dos volumes finito TPFA
% mpfad     --> (MPFA-D)
% lfvLPEW   --> método linear basedo no método não linear usando LPEW (MPFA-HD), ter cuidado linha 52 e 54 do preNLFV
% lfvHP     --> (MPFA-H)
% lfvEB     --> método completamente baseado na face (MPFA-BE).
% nlfvLPEW  --> (NLFV-PP)
% nlfvDMPSY --> método não linear que preserva DMP baseado no artigo (Gao e Wu, 2013) e (Sheng e Yuan, 20...)
% nlfvHP
% nlfvPPS
pmetodo='nlfvLPEW';
%% metodo de interação: picard, newton, broyden, secant,
% método de itereção proprio de métodos não lineares iterfreejacobian,iterdiscretnewton, JFNK
% iteration='iterdiscretnewton';
% iteration='iterbroyden';
% iteration='JFNK';
% iteration='fullpicard';
% iteration='MPE'; 
% iteration='RRE'; % picard com acelerador rank reduced extrapolation
iteration='AA'; % picard com aceleracao de Anderson
%iteration='iterhybrid';
%% defina o ponto de interpolacao
interpol='LPEW2';
%interpol='LPEW1';
%interpol='eLS';
%interpol='eLPEW2';
%% correcao dos pontos harmonicos
% digite 'yes' ou 'no'
correction= 'no';
%% digite segundo o benchmark
% procure no "adequapermeab.m" o caso que deseja rodar e logo digite o nome
% do caso
benchmark='edqueiroz'; 
%% com termo gravitacional
% com termo gravitacional 'yes' ou 'no'
gravitational='no';
%% adequação dos flags
%nflag= calflag(pmetodo);
%% este flag só use quando o problema é Buckley-Leverett com fluxo imposto
% na face com flag 202, porém a saturação será imposto na face
auxflag=202;

%% adequação das permeabilidades, otros entes fisico-geometricos segundo o bechmark
[elem,kmap,normKmap,solanal,bedge,fonte,velanal,grav]=benchmarks(benchmark,...
    kmap,elem,bedge);

% F faces na vizinhanca de um elemento
% V 
% N
[F,V,N]=elementface;
%% pre-processador local
[pointarmonic,parameter,gamma,p_old,tol,nit,er,nflagface,nflagno,...
    weightDMP,Hesq,Kde,Kn,Kt,Ded,auxface,calnormface,gravresult,gravrate]=...
    preNLFV(kmap,N,pmetodo,benchmark,bedge,grav,gravitational,correction);
nflag=nflagno;
% não habilite
%[aroundface]=aroundfacelement(F,pointarmonic);
%% calculo das mobilidades
% mobilidade sera utilizado unitario porque este simulador "in house" somente
% trata problemas de escoamento monofasico
mobility=zeros(size(bedge,1)+size(inedge,1),1);
mobility(:)=1;
%%  Calculo da pressao pelos método lineares e nao-lineares
[p,errorelativo,flowrate,flowresult,tabletol,coercividade]=solverpressure(...
    kmap,nflagface,nflagno,fonte, tol,...
    nit,p_old,mobility,gamma,wells,parameter,pmetodo,auxflag,interpol,...
    Hesq, Kde, Kn, Kt, Ded,weightDMP,auxface,...
    benchmark,iteration,nflag,calnormface,N,gravresult,gravrate);
% pos-processador no visit
% postprocessor(full(abs(p-solanal)),1)
% postprocessor(full(p),2)
% postprocessor(solanal,3)
%% calculo do erro, pressao maxima e minima
[erropressure,errovelocity,pressuremax,pressuremin]=errorateconv(solanal, p, velanal,flowrate,...
    erromethod,benchmark)

%hold on
%grid on
%% plota os erros
%plot(tabletol(:,1),log(tabletol(:,2)),'s-m','LineWidth',2)
%ylabel('log(Error)')
%xlabel('Number iterations')

