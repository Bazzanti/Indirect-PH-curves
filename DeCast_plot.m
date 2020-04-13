function [] = DeCast_plot(n,Vx,Vy)
%date le coordinate dei vertici del poligono di controllo Vx e Vy 
%usa DeCasteljau per visualizzare la curva di Bézier

%inizializzazione variabili
t=linspace(0,1,100);
Cx=zeros(1,100);
Cy=zeros(1,100);
%calcolo curva
for i=1:100 
    Cx(i)=deCasteljau(n,Vx,t(i));
    Cy(i)=deCasteljau(n,Vy,t(i));
end
%disegno
plot(Cx,Cy,Vx,Vy,'-*') 

end