function [CP] = Main_PH5 (pI,pF,T0,T1) 

%Ph5 con hermite interpolation
d0 = complex(T0(1), T0(2));
d1 = complex(T1(1), T1(2));

pI = complex(pI(1), pI(2));
pF = complex(pF(1), pF(2));

% -- pag. 530 (25.18)
w0 = sqrt(d0);
w2 = sqrt(d1);

coeff_1 = 2/3;
coeff_2 = w0 + w2;
coeff_3 = w0^2 + w2^2 + (w2*w0)/3 - 5*(pF-pI);

w1 = roots([coeff_1 coeff_2 coeff_3]);

useSol_1    = false;
 if(useSol_1)
                w1 = w1(1);
            else
                w1 = w1(2);
 end
            


p1 = pI + 1/5 * w0^2;
p2 = p1 + 1/5 * w0*w1;
p3 = p2 + 1/15 * (2*w1^2 + w0*w2);
p4 = p3 + 1/5 * w1*w2;
p5 = p4 + 1/5 * w2^2;
CPcmplx = [pI p1 p2 p3 p4 p5];
CP = [real(CPcmplx) ;  imag(CPcmplx)];

% -- save w,u,v coefficents (building polynomials W complex, and coeffs not complex)
W = [w0 w1 w2];
U = [real(w0) real(w1) real(w2)];
V = [imag(w0) imag(w1) imag(w2)];

% -- save sigmas
sigma(1) = U(1)^2 + V(1)^2;
sigma(2) = U(1)*U(2) + V(1)*V(2);
sigma(3) = 2/3* (U(2)^2 + V(2)^2) + 1/3*(U(1)*U(3) + V(1)*V(3));
sigma(4) = U(2)*U(3) + V(2)*V(3);
sigma(5) = U(3)^2 + V(3)^2;

% -- save arc length coefficients s
s(1) = 0;
s(2) = 1/5 * (sigma(1));
s(3) = 1/5 * (sigma(1) + sigma(2));
s(4) = 1/5 * (sigma(1) + sigma(2) + sigma(3));
s(5) = 1/5 * (sigma(1) + sigma(2) + sigma(3) + sigma(4));
s(6) = 1/5 * (sigma(1) + sigma(2) + sigma(3) + sigma(4) + sigma(5));

end






%{
%Trasformazione rigida a formato standard
L=P5-P0;
ang = atan2(P5(2)-P0(2),P5(1)-P0(1)); 
[P0] = trasf_rigida (P0,ang,supP0);
[P5] = trasf_rigida (P5,ang,supP0);
[T0] = trasf_rigida (T0,ang,[0,0]);
[T1] = trasf_rigida (T1,ang,[0,0]);
P5(1)=P5(1)/L; %P3=[1,0]

%CALCOLO punti tramite derivate e coordinate
%derivate
d0=T0(2)/T0(1);
d1=T1(2)/T1(1);

R0=sqrt(d0/d1);
R1=-sqrt(d0/d1);

%4 soluzioni
alfa1 = (3*(1-R0) + sqrt( 3*(1-R0)^2 - 4*(6*R0^2+2*R0+6-30/d1)))/2;
alfa2 = (3*(1-R0) - sqrt( 3*(1-R0)^2 - 4*(6*R0^2+2*R0+6-30/d1)))/2;
alfa3 = (3*(1-R1) + sqrt( 3*(1-R1)^2 - 4*(6*R1^2+2*R1+6-30/d1)))/2;
alfa4 = (3*(1-R1) - sqrt( 3*(1-R1)^2 - 4*(6*R1^2+2*R1+6-30/d1)))/2;

%%%%%%%%%%%%%%%%%%%
u(1) = (alfa1+sqrt(alfa1^2-4*R0))/2;
u(2) = (alfa1-sqrt(alfa1^2-4*R0))/2;
u(3) = (alfa2+sqrt(alfa2^2-4*R0))/2;
u(4) = (alfa2-sqrt(alfa2^2-4*R0))/2;

for i=1:4
a(i)=u()


end
%}