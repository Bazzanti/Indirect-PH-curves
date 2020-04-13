function [P1,P2] = indPH (P0,P3,T0,T1)
%indPH: calcola i punti P1 e P2 per una curva indirect PH
%
%Con l'algoritmo trovo i punti intermedi P1 e P2 tra i punti di controllo
%[P0,P1,P2,P3] per una curva di Bezier cubica
%
% input: P0,P3, primo e ultimo punto della ph cubica
%         T0,T1 vettori direzione del primo e ultimo punto   
%output: P1,P2

%uso P0=[0,0] e P3=[L,0] già riparametrizzati sotto trasformazione rigida
L0 = P3(1);

%I= intersezione tra R0 e R1 semirette
%R0=@(u) P0+abs(u)*T0; u>0
%R1=@(v) P3-abs(v)*T1; v>0
v = (+P3(1)*T0(2))/(-T0(1)*T1(2)+T1(1)*T0(2));
u = (-T1(2)*v)/T0(2);

I = u*T0;
%nella costruzione del teorema viene usato I(2)>0, per cui lo imposto
%in tale modo e modifico il risultato finale con segno inverso se  I(2)<0
%suppI=1; %lo imposto in precedenza positivo
%if I(2)<0
%    suppI=I(2);
%    I(2)=-I(2);
%end

%angoli alfa e beta tra P0,T0 e P3,T1
Alfa = atan2(T0(2),T0(1)); %0,pi angolo tra T0 e P3-P0
Beta = atan2(-T1(2),T1(1)); %0,pi angolo tra T1 e P0-P3

%dato che ho usato I(2)>0, se inizialmente era negativo allora
%agli angoli Alfa e Beta diventano positivi
%if suppI < 0
%   Alfa = -Alfa;
%   Beta = -Beta;
%end   
  
%calcolo l'h ottimale usando F
m = I(1)/L0;
n = I(2)/L0;

F = @(h) 2*(6*m^2-3*m+6*n^2+1)*h^4 ...
     +(12*m^2-27*m+12*n^2+11)*h^3 ...
     +18*(-2*m+1)*h^2+(-12*m^2-3*m+12*n^2+4)*h ...
     + 2*(-6*m^2+9*m-6*n^2-4);

x0 = [0.49 , 2.04]; % intervallo iniziale, teorema
h = fzero(F,x0); %trovo lo 0 di F
h=5;
%trova L1=P1-P0 e L3=P3-P2 tramite shape control
L1 = (2*L0*sin(Beta)) / ((h+2)*(sin(Alfa+Beta))) ;
L3 = (2*h*L0*sin(Alfa)) / ((2*h+1)*(sin(Alfa+Beta))) ;

%calcolo P1 e P2
P1 = P0+L1*T0/norm(T0);
P2 = P3-L3*T1/norm(T1) ;
    
%dato che ho usato I(2)>0, se inizialmente era negativo allora
%i punti P1 e P2 sono venuti nell'ordinata negativa
%if suppI < 0
%   P1(2) = -P1(2)
%   P2(2) = -P2(2)
%end   
    
end

