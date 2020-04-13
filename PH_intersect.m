function [Q1,Q2,Pm] = PH_intersect (P3,T0,T1)
% input: P0=[0,0],P3=[L,0] riparametrizzati, primo e ultimo punto della ph cubica
%         T0,T1 vettori direzione del primo e ultimo punto   
%output: Q1, Q2 punti iniziali trovati con il metodo del rettangolo
%        Pm punto di intersezione tra Q1Q2 e P0P3
%spezzo la curva in 2 e trovo 3 nuovi punti per creare 2 indirect-Ph

%Q1(u)=P0+uT0    Q2(v)=P3-vT1    u,v>0
%angoli alfa e beta tra P0,T0 e P3,-T1  tra [-pi,pi]

l=P3(1);
Alfa = atan2(T0(2),T0(1));  %atan2(Y,X) [-pi,pi]
Beta = atan2(-T1(2),-T1(1)); %-T1         
    
%viene usato arctan(2) perche il metodo usato usa il rettangolo con lato
%uno il doppio dell altro

%punto dipendente da Alfa 
if -atan(2) <=Alfa && Alfa< atan(2) 
    Q1(1) = l/4;
    Q1(2) = l*T0(2) / (4*T0(1));
elseif atan(2) <=Alfa && Alfa< pi-atan(2)
    Q1(2) = l/2;
    Q1(1) = l*T0(1) / (2*T0(2));
elseif (pi-atan(2) <=Alfa && Alfa< pi)  ||  (-pi <=Alfa && Alfa< -pi+atan(2))
    Q1(1) = -l/4;
    Q1(2) = -l*T0(2) / (4*T0(1));
elseif -pi+atan(2) <=Alfa && Alfa< -atan(2)
    Q1(2) = -l/2;
    Q1(1) = -l*T0(1) / (2*T0(2));    
end

%punto dipendente da Beta
if -atan(2) <=Beta && Beta< atan(2) 
    Q2(1) = 5*l/4;
    Q2(2) = l*T1(2) / (4*T1(1));
 
elseif atan(2) <=Beta && Beta< pi-atan(2)
    Q2(2) = l/2;
    Q2(1) = l + l*T1(1) / (2*T1(2));
    
elseif (pi-atan(2) <=Beta && Beta<= pi)  ||  (-pi <=Beta  && Beta< -pi+atan(2))  
    Q2(1) = 3*l/4;
    Q2(2) = -l*T1(2) / (4*T1(1));
    
elseif -pi+atan(2) <=Beta && Beta< -atan(2)
    Q2(2) = -l/2;
    Q2(1) = l-l*T1(1) / (2*T1(2));
end

%punto Pm intersezione tra Q1Q2 e bisettore di P0P3
Pm = [l/2 , ( (l/2-Q1(1)) * (Q2(2)-Q1(2)) / (Q2(1)-Q1(1))) + Q1(2)];   
    
end