function [PP] = trasf_rigida (P,Alfa,P0)
%trasf_rigida: Dato P cerco un nuovo PP dopo trasformazione rigida
%              attraverso Alfa angolo, P0 vettore.
%               
%
%input: P punto da trasformare
%       Alfa angolo di rotazione
%       P0 vettore di traslazione
%output: P trasformato

%traslazione
PP = [P(1)-P0(1),P(2)-P0(2)];

%rotazione
R = [cos(Alfa) , sin(Alfa) ;-sin(Alfa) , cos(Alfa)];
PP = R*PP';
PP = PP';

end