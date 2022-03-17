function [omega]=omegacalculate(E,EE,netas,no)

global  esurn2 

O=zeros(esurn2(no+1)-esurn2(no),3); % vetor de baricentro na vizinhança do nó "ni".
nec=size(O,1);

for k=1:nec,

    omega(k,1)=((netas(k,1)+netas(k,2))*(E(k,1))/(EE(k,1)+EE(k,2)))+ ((netas(k,3)+netas(k,4))*(E(k,2))/(EE(k,3)+EE(k,4)));
    
end


end