function [C] = deCasteljau( n,V,t )
%algoritmo di deCasteljau
%n=grado della curva di Bezier
%V vertici di controllo V0,...,Vn
%t parametro dove calcolare la curva

for i=1:n+1
        Q(i)=V(i);
end
for k=1:n    %passo n+1 non esiste perche risulta i=1:0
    for i=1:n-k+1
        Q(i)=(1.0-t)*Q(i)+t*Q(i+1);
    end
end
C=Q(1);
end

