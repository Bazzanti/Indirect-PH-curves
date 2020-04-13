function [der,der2] = Der_bezier5 (V0,V1,V2,V3,V4,V5)
%calcolo derivata di una curva di Bezier di 5 grado dati i vertici
%input: V0,V1,V2,V3,V4,V5 vertici Vx o Vy di una curva di Bezier
%output: funzione della derivata prima (der) 
%        e della derivata seconda (dev2)

%derivata prima
der = @(u) 5*( (u.^4)*(-V0+5*V1-10*V2+10*V3-5*V4+V5) + (u.^3)*(4*V0+-16*V1+24*V2-16*V3+4*V4)+ (u.^2)*(-6*V0+18*V1-18*V2+6*V3) + (u)*(4*V0-8*V1+4*V2) -V0 + V1);

%derivata seconda
der2 = @(u) 5*( 4*(u.^3)*(-V0+5*V1-10*V2+10*V3-5*V4+V5) + 3*(u.^2)*(4*V0+-16*V1+24*V2-16*V3+4*V4)+ 2*(u)*(-6*V0+18*V1-18*V2+6*V3) + (4*V0-8*V1+4*V2) );
end