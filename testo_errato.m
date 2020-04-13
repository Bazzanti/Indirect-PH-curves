%%
%punto Pm intersezione tra Q1Q2 e P0P3
Pm = [(-Q1(2) * (Q2(1)-Q1(1)) / (Q2(2)-Q1(2))) + Q1(1), 0];


 %%
 else%u,v non sono corretti
    %spezzo la curva in 2 creando 3 nuovi punti
    [Q1,Q2,Pm] = PH_intersect (P3,T0,T1);
    Q1Q2=Q2-Q1;
    %PRIMA CURVA
    %prima ind-ph usando il set (P0,T0;Pm,Q1Q2)
    %P0=[0,0] e Pm=[l,0] quindi non li devo riparametrizzare
    [P1_1,P2_1] = indPH(P0,Pm,T0,Q1Q2);  
    %trasformazione rigida inversa
    [P0_1] = trasf_rigida_inversa (P0,ang,supP0);
    [P1_1] = trasf_rigida_inversa (P1_1,ang,supP0);
    [P2_1] = trasf_rigida_inversa (P2_1,ang,supP0);
    [P3_1] = trasf_rigida_inversa (Pm,ang,supP0);    
    
    %primo plot
    Vx1=[P0_1(1),P1_1(1),P2_1(1),P3_1(1)];
    Vy1=[P0_1(2),P1_1(2),P2_1(2),P3_1(2)];
    hold on
    DeCast_plot(3,Vx1,Vy1)
  
    %SECONDA CURVA
    %seconda ind-PH
    %riparametrizzazione di Q1 Q2 Pm P3
    %Pm e P3 sono nell'asse delle ascisse quindi devo solo traslare
    if P3(1)>=Pm(1)
        [P0_2] = trasf_rigida (Pm,0,Pm);
        [P3_2] = trasf_rigida (P3,0,Pm);
        [P1_2,P2_2] = indPH (P0_2,P3_2,Q1Q2,T1);
        %ora faccio la traslazione inversa di Pm
        [P0_2] = trasf_rigida_inversa (P0_2,0,Pm);    
        [P1_2] = trasf_rigida_inversa (P1_2,0,Pm);
        [P2_2] = trasf_rigida_inversa (P2_2,0,Pm);
        [P3_2] = trasf_rigida_inversa (P3_2,0,Pm);
    
  
     elseif Pm(1)>P3(1) %rotazione di 180
        [P0_2] = trasf_rigida (Pm,0,Pm); %Pm va a 0
        P3_2 = [-(P3(1)-Pm(1)) , P3(2)]; %P3 va ribaltato
        [Q1Q2_2] = trasf_rigida (Q1Q2,pi,[0,0]);
        [T1_2] = trasf_rigida (T1,pi,[0,0]);
        [P1_2,P2_2] = indPH (P0_2,P3_2,Q1Q2_2,T1_2); %vettori invertiti
        
        %in questo caso particolare devo prima effettuare la rotazione
        %ed in seguito la traslazione perché ho la possibilità di un 
        %cambio di quadrante  che provoca errori di rotazione
        
        %ora faccio la traslazione inversa di +180 gradi
        [P0_2] = trasf_rigida_inversa (P0_2,pi,[0,0]);    
        [P1_2] = trasf_rigida_inversa (P1_2,pi,[0,0]);
        [P2_2] = trasf_rigida_inversa (P2_2,pi,[0,0]);
        [P3_2] = trasf_rigida_inversa (P3_2,pi,[0,0]);
        %ora faccio la traslazione inversa di Pm
        [P0_2] = trasf_rigida_inversa (P0_2,0,Pm);    
        [P1_2] = trasf_rigida_inversa (P1_2,0,Pm);
        [P2_2] = trasf_rigida_inversa (P2_2,0,Pm);
        [P3_2] = trasf_rigida_inversa (P3_2,0,Pm);
        end   
            
        
        
    %inversa trasformazione rigida di P0+angolo:
    [P0_2] = trasf_rigida_inversa (P0_2,ang,supP0);
    [P1_2] = trasf_rigida_inversa (P1_2,ang,supP0);
    [P2_2] = trasf_rigida_inversa (P2_2,ang,supP0);
    [P3_2] = trasf_rigida_inversa (P3_2,ang,supP0);
    
    Vx2=[P0_2(1),P1_2(1),P2_2(1),P3_2(1)];
    Vy2=[P0_2(2),P1_2(2),P2_2(2),P3_2(2)];
    DeCast_plot(3,Vx2,Vy2)
 

                 
%%
   [Q1] = trasf_rigida_inversa (Q1,ang,supP0);
   [Q2] = trasf_rigida_inversa (Q2,ang,supP0);   
   plot(Q1(1),Q1(2),'-s','MarkerSize',10)
   plot(Q2(1),Q2(2),'-s','MarkerSize',10)


 %%
    
        