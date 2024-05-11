clear;
clc;

% Condições Globais do problema
L = 3;
V = 100;
h = 0.15;
Pho = 1.25;
Gamma = 1.4;
K_ar = 0.026;
Cp_ar = 1002;
d = L*1/2;
H = 2*L;
T_dentro = 25;
T_fora = 20;
T_motor = 80;

%Indices de aplicação das diferenças finitas

Lambda = 1.85;
Erro = 0.01;
deltaX = 0.5;
deltaY = deltaX;


%

