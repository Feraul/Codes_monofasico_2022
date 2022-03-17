function [ w,s] = Pre_ELPEW_2(kmap,N,g)
global coord  esurn2
% devolve os pesos "w" cujos elementos s�o organizados em um vetor
%Retorna todos os par�metros necess�rios �s express�es dos fluxos.
apw=ones(size(coord,1),1);
r=zeros(size(coord,1),2);
s=0;
epsilon=0.5;

for No=1:size(coord,1)
    % calcula
    
   O=zeros(esurn2(No+1)-esurn2(No),3); 
   
   %[EE]=xxicalculate(No,kmap,epsilon);
   %[netas,E]=netascalculate(No,kmap,epsilon);
   %[omegas]=omegacalculate(E,EE,netas,No);
   
   [netas,E,EE]=parameters_weight(No,kmap,epsilon);
   [omegas]=aux_omegacalculate(E,EE,netas,No);
    for k=0:size(O,1)-1
        w(apw(No)+k,1)=omegas(k+1)/sum(omegas); %Os pesos fazem sentido
    end

    apw(No+1)=apw(No)+size(O,1);
    % calculando os pesos relativos a condi��o de contorno de Neumann
    % interpola�ao das press�es nos contornos de Neumann
%     vetor = nsurn1(nsurn2(No) + 1:nsurn2(No + 1));
%     comp1 = N(No,1);
%     comp2 = N(No,length(vetor));
%     if comp1 > size(inedge,1) && comp2 > size(inedge,1)
%         a = bcflag(:,1) == bedge(comp1 - size(inedge,1),5);
%         s1 = find(a == 1);
%         b = bcflag(:,1) == bedge(comp2 - size(inedge,1),5);
%         s2 = find(b == 1);
%         
%         s(No,1) = -(1/sum(lambda))*(r(No,1)*(bcflag(s1,2)+g(comp1,1)) + ...
%             r(No,2)*(bcflag(s2,2)+g(comp1,1)));
%     end  %End of IF
end
end

