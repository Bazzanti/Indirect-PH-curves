function [I]=inters2rette(m1,k1,m2,k2)
% [X, Y]=inters2retta(m1,k1,m2,k2)
% Calcola le coordinate del puntio di intersezione di due rette definite
% come:
%   retta 1 ==> y=m1*x+k1
%   retta 2 ==> y=m2*x+k2
%
% se le rette non hanno un punto di intersezione ==> X=NaN, Y=NaN
%
%
X=NaN;
Y=NaN;
D=(-m1+m2);
if(D ~= 0)
   X=(k1-k2)/D;
   Y=(-m1*k2+m2*k1)/D;
end
I=[X,Y];
end