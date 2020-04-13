%%1curva

clear

P0=[0,0];
P3=[1,0];
T0=[1,2];
T1=[1,-1];
figure()
Main_indPH (P0,P3,T0,T1) 
%%
%1curva

P0=[0,0];
P3=[1,0];
T0=[1,1];
T1=[1,-1];
figure()
Main_indPH (P0,P3,T0,T1)


%% 1curva
P0=[0,0];
P3=[3,-1];
T0=[1/2,-1];
T1=[1,1];
figure()
Main_indPH (P0,P3,T0,T1) 

%% 1curve

P0=[1,3];
P3=[4,2];
T0=[2,-1];
T1=[2,3];
figure()
Main_indPH (P0,P3,T0,T1) 
%% 2curve 

P0=[1,1];
P3=[3,1];
T0=[1/2,-1];
T1=[0,-1];
figure()
Main_indPH (P0,P3,T0,T1) 
%%

P0=[0,0];
P3=[1,0];
T0=[-1,2];
T1=[-1,-1];
figure()
Main_indPH (P0,P3,T0,T1) 
%% OFFSET1

P0=[1,1];
P3=[3,1];
T0=[0.5,-1];
T1=[0,-1];
figure()
Main_indPH (P0,P3,T0,T1) 
%%
%offset2pti
%P0=[1,1];
%P3=[2,0.667];
%T0=[0.5,-1];
%T1=[0.29,0.38];
P0=[2,0.667];
P3=[3,1];
T0=[0.29,0.38];
T1=[0,-1];
figure()
Main_indPH (P0,P3,T0,T1) 
%%
%P=[2,3];
%Alfa=1.2;
%P0=[2,1];
%[PP] = trasf_rigida (P,Alfa,P0)
%PP1= trasf_rigida_inversa (PP,Alfa,P0)

