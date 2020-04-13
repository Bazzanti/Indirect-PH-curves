function [] = Main_indPH (P0,P3,T0,T1) 
%indPH: calcola  una curva indirect PH
%
%Con l'algoritmo trovo i punti intermedi P1 e P2 tra i punti di controllo
%[P0,P1,P2,P3] per una curva di Bezier cubica
%
% input: P0,P3, primo e ultimo punto della ph cubica
%         T0,T1 vettori direzione del primo e ultimo punto   
%output:curva di Bezier cubica


%(1) Trasformazione rigida 
supP0=P0;
ang = atan2(P3(2)-P0(2),P3(1)-P0(1));  %atan2(Y,X) [-pi,pi]
[P0] = trasf_rigida (P0,ang,supP0);
[P3] = trasf_rigida (P3,ang,supP0);
[T0] = trasf_rigida (T0,ang,[0,0]);
[T1] = trasf_rigida (T1,ang,[0,0]);

%uso P0 e P3 riparametrizzati sotto trasformazione rigida
%P0 = [0,0]; P3 = [L,0]; 
%CALCOLO RAGGI 
v = (-P0(1)*T0(2)-T0(1)*(P3(2)-P0(2))+P3(1)*T0(2))/(-T0(1)*T1(2)+T1(1)*T0(2));
u = (-T1(2)*v + P3(2)-P0(2))/T0(2);

%CALCOLO CURVA
if u>0 && v>0 && u~=inf && v~=inf % u,v sono corretti quindi posso calcolare la curva indPH

    [P1,P2] = indPH (P0,P3,T0,T1);
    %(1i) inversa trasformazione rigida :
    [P0] = trasf_rigida_inversa (P0,ang,supP0);
    [P1] = trasf_rigida_inversa (P1,ang,supP0);
    [P2] = trasf_rigida_inversa (P2,ang,supP0);
    [P3] = trasf_rigida_inversa (P3,ang,supP0);
    
    Vx = [P0(1),P1(1),P2(1),P3(1)];
    Vy = [P0(2),P1(2),P2(2),P3(2)];
    %Plot con decasteljau 
    
    DeCast_plot(3,Vx,Vy)
    axis equal
    
    %RHO0 RHO1
    %R0= 3*(P1(1)-P0(1))/T0(1)
    %R1= 3*(P3(1)-P2(1))/T1(1)
    
    %CURVATURA
    
    %derivate funz handle
    [derX,dersecX] = Der_bezier3 (P0(1),P1(1),P2(1),P3(1));
    [derY,dersecY] = Der_bezier3 (P0(2),P1(2),P2(2),P3(2));
     
    %%%%%%%%%%%%
    %LENGTH
    %ds=@(x) sqrt((derX(x)).^2 + (derY(x)).^2);
    %Length = integral(ds,0,1)
     
    %%%%%%%%%%%%%%
    %OFFSET
    %
    hold on
    [Cx,Cy]=DeCast_curve(3,Vx,Vy);
    j=linspace(0,1);
    Offx=zeros(1,100);
    Offy=zeros(1,100);
    d=-0.1; %distanza offset
    for i=1:100
        Offx(i)=Cx(i) + d* (-derY(j(i)))./sqrt(derX(j(i)).^2+derY(j(i)).^2); 
        Offy(i)=Cy(i) + d* (derX(j(i)))./sqrt(derX(j(i)).^2+derY(j(i)).^2); 
    end
    plot(Offx,Offy)
    %
    %%%%%%%%%%%%%%
    
    %calcolo curvatura
    t=linspace(0,1);
    K=zeros(1,100);
    for i=1:100
        K(i)= abs((derX(t(i)).*dersecY(t(i))-derY(t(i)).*dersecX(t(i))))./(derX(t(i)).^2+derY(t(i)).^2).^(3/2); 
    end

    
    %plot curvatura
    figure()
    plot(t,K)
 
else%u,v non sono corretti
    %spezzo la curva in 2 creando 3 nuovi punti
    [Q1,Q2,Pm] = PH_intersect (P3,T0,T1);
    Q1Q2=Q2-Q1;
    alfa = atan2(Pm(2),Pm(1)); %angolo di Pm
    beta = atan2(P3(2)-Pm(2),P3(1)-Pm(1)); % angolo P3-Pm
   %no beta = atan2(P3(2),P3(1));
    
    %PRIMA CURVA
    %prima ind-ph usando il set (P0,T0;Pm,Q1Q2)
    
    %(2) Trasformazione P0=[0,0] e Pm=[l/2,y] 
    [Pm_1] = trasf_rigida (Pm,alfa,[0,0]);
    [T0_1] = trasf_rigida (T0,alfa,[0,0]);
    [Q1Q2_1] = trasf_rigida (Q1Q2,alfa,[0,0]);
    
    %Trovo P1 e P2 della prima curva
    [P1_1,P2_1] = indPH(P0,Pm_1,T0_1,Q1Q2_1);
    %(2i) Trasformazione inversa 
    [P1_1] = trasf_rigida_inversa (P1_1,alfa,[0,0]);
    [P2_1] = trasf_rigida_inversa (P2_1,alfa,[0,0]);
    %(1i) Trasformazione rigida inversa finale
    [P0_1] = trasf_rigida_inversa (P0,ang,supP0);
    [P1_1] = trasf_rigida_inversa (P1_1,ang,supP0);
    [P2_1] = trasf_rigida_inversa (P2_1,ang,supP0);
    [P3_1] = trasf_rigida_inversa (Pm,ang,supP0);

    %primo plot
    Vx1=[P0_1(1),P1_1(1),P2_1(1),P3_1(1)];
    Vy1=[P0_1(2),P1_1(2),P2_1(2),P3_1(2)];
    hold on
    
    %Rho1 Rho2
    %R0= 3*(P1_1(1)-P0_1(1))/T0_1(1)
    %R1= 3*(P3_1(1)-P2_1(1))/Q1Q2_1(1)
    
    
    DeCast_plot(3,Vx1,Vy1)
    axis equal
    
    %SECONDA CURVA
    %seconda ind-PH usando il set (Pm,Q2Q1;P3,T1)
   
    %(3) Trasformazione Pm=[l/2,y] e P3=[l,0] 
    [P0_2] = trasf_rigida (Pm,0,Pm); % va a 0,0
    [P3_2] = trasf_rigida (P3,beta,Pm);
    [Q1Q2_2] = trasf_rigida (Q1Q2,beta,[0,0]);
    [T1_2] = trasf_rigida (T1,beta,[0,0]);
    
    %Trovo P1 e P2 della seconda curva
    [P1_2,P2_2] = indPH (P0_2,P3_2,Q1Q2_2,T1_2);
    
    %%(3i) Trasformazione inversa 
    [P1_2] = trasf_rigida_inversa (P1_2,beta,Pm); 
    [P2_2] = trasf_rigida_inversa (P2_2,beta,Pm);

    %(1i) inversa trasformazione rigida di P0+angolo:
    [P0_2] = trasf_rigida_inversa (Pm,ang,supP0);
    [P1_2] = trasf_rigida_inversa (P1_2,ang,supP0);
    [P2_2] = trasf_rigida_inversa (P2_2,ang,supP0);
    [P3_2] = trasf_rigida_inversa (P3,ang,supP0);
    
    Vx2=[P0_2(1),P1_2(1),P2_2(1),P3_2(1)];
    Vy2=[P0_2(2),P1_2(2),P2_2(2),P3_2(2)];
    
    %Rho1 Rho2
    %R0= 3*(P1_2(1)-P0_2(1))/Q1Q2_2(1)
    %R1= 3*(P3_2(2)-P2_2(2))/T1_2(2)
    
    DeCast_plot(3,Vx2,Vy2) 
    axis equal
    hold off
    
    %CURVATURA
    %curva 1
    [derX_1,dersecX_1] = Der_bezier3 (P0_1(1),P1_1(1),P2_1(1),P3_1(1));
    [derY_1,dersecY_1] = Der_bezier3 (P0_1(2),P1_1(2),P2_1(2),P3_1(2));
    figure()
    t=linspace(0,1);
    K1=zeros(1,100);
    for i=1:100
        K1(i)= abs((derX_1(t(i)).*dersecY_1(t(i))-derY_1(t(i)).*dersecX_1(t(i))))./(derX_1(t(i)).^2+derY_1(t(i)).^2).^(3/2); 
    end
    plot(t,K1) 
    hold on
     %ds=@(x) sqrt((derX_1(x)).^2 + (derY_1(x)).^2);
     %Length = integral(ds,0,1)
    %curva 2
    [derX_2,dersecX_2] = Der_bezier3 (P0_2(1),P1_2(1),P2_2(1),P3_2(1));
    [derY_2,dersecY_2] = Der_bezier3 (P0_2(2),P1_2(2),P2_2(2),P3_2(2));
    t2=linspace(0,1);
    K2=zeros(1,100);
    for i=1:100
        K2(i)= abs((derX_2(t2(i)).*dersecY_2(t2(i))-derY_2(t2(i)).*dersecX_2(t2(i))))./(derX_2(t2(i)).^2+derY_2(t2(i)).^2).^(3/2) ; 
    end 
    plot(1+t2,K2)
    hold off
    % ds=@(x) sqrt((derX_2(x)).^2 + (derY_2(x)).^2);
     %Length = integral(ds,0,1)
end

end