clear
%init

%P0=[0,0];
%P3=[1,0];
%T0=[0.13,0.26]; %T0=[1,2];
%T1=[1.27,-1.27]; %T1=[1,-1];

%P0=[1,3];
%P3=[4,2];
%T0=[0.95,-0.5]; %T0=[2,-1];
%T1=[0.39,0.58]; %T1=[2,3];
%primapt 2curve
%P0=[1,1];
%P3=[2,0.667];
%T0=[0.44,-0.29];
%T1=[1.06, 1.28];

P3=[3,1];
P0=[2,0.667];
%T0=[0.25,0.33];
%T1=[0, -1];
T0=[0.5,-1];
T1=[0.5,-1];

%drawArrow = @(x,y) quiver( x(1),y(1),x(2)-x(1),y(2)-y(1),0 ) ;   
%drawArrow([0,1],[0,1]); hold on


%PH5
CP = Main_PH5 (P0,P3,T0,T1);

Vx=[CP(1,1),CP(1,2),CP(1,3),CP(1,4),CP(1,5),CP(1,6)];
Vy=[CP(2,1),CP(2,2),CP(2,3),CP(2,4),CP(2,5),CP(2,6)];

hold on
DeCast_plot(5,Vx,Vy);
%Main_indPH (P0,P3,T0,T1);
axis equal

%Curvatura
[derX,dersecX] = Der_bezier5 (CP(1,1),CP(1,2),CP(1,3),CP(1,4),CP(1,5),CP(1,6));
[derY,dersecY] = Der_bezier5 (CP(2,1),CP(2,2),CP(2,3),CP(2,4),CP(2,5),CP(2,6));

figure()
title('Curvatura')
t=linspace(0,1);
K=zeros(1,100);

for i=1:100
    K(i)= abs((derX(t(i)).*dersecY(t(i))-derY(t(i)).*dersecX(t(i))))./(derX(t(i)).^2+derY(t(i)).^2).^(3/2); 
 
end
plot(t,K)

 ds=@(x) sqrt((derX(x)).^2 + (derY(x)).^2);
 Length5 = integral(ds,0,1)
     
%ind PH

%figure()
%title('Curvatura')
%t=linspace(0,1);
%K=zeros(1,100);
%K1=zeros(1,100);
%K2=zeros(1,100);
%K3=zeros(1,100);
%K4=zeros(1,100);
%for i=1:100
%    K(i)= abs((derX(t(i)).*dersecY(t(i))-derY(t(i)).*dersecX(t(i))))./(derX(t(i)).^2+derY(t(i)).^2).^(3/2); 
%    K1(i)= dersecX(t(i));
%     K2(i)= derX(t(i));
%      K3(i)= dersecY(t(i));
%       K4(i)= derY(t(i));
%end
%plot(t,K)
%figure()
%plot(t,K1)
%figure()
%plot(t,K2)
%figure()
%plot(t,K3)
%figure()
%plot(t,K4)

%ind PH



