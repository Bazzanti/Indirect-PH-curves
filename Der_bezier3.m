function [der,der2] = Der_bezier3 (V0,V1,V2,V3)
%calcolo derivata di una curva di Bezier di 3 grado dati i vertici
%input: V0,V1,V2,V3 vertici Vx o Vy di una curva di Bezier
%output: funzione della derivata prima (der) 
%        e della derivata seconda (dev2)

%derivata prima
der = @(u) 3*( (u.^2)*(-V0+3*V1-3*V2+V3) + u.*(2*V0-4*V1+2*V2)+V1-V0);

%derivata seconda
der2 = @(u) 6*u.*(-V0+3*V1-3*V2+V3) + 6*(V0-2*V1+V2);

end