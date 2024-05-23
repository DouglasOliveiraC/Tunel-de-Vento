% Propriedades da barra
m = 5;               % Massa 
c = 0.1;             % Coeficiente de amortecimento 
k = 10;             % Rigidez da barra
g = 9.81;            % Aceleração da gravidade

% Condições iniciais causadas pelo peso G
G = m * g;           % Peso da barra
u0 = G/k;            % Deslocamento inicial devido ao peso
a0 = (G - k*u0) / m; % Aceleração inicial após o corte
v0 = 0;              % Velocidade inicial (pois estava em repouso antes do corte)

% Parâmetros de Newmark
gamma = 0.5;
beta = 0.25;
Dt = 0.01;             % Passo de tempo (por exemplo, 0.01 s)

% Inicialização dos vetores
u = u0;
v = v0;
a = a0;
time = 0:Dt:100;        % Simulação por 10 segundos, por exemplo

u_history = zeros(size(time));
v_history = zeros(size(time));
a_history = zeros(size(time));

u_history(1) = u0;
v_history(1) = v0;
a_history(1) = a0;

% Método de Newmark
for i = 1:length(time)-1
    % Fórmulas de Predição
    v_pred = v + (1-gamma)*Dt*a;
    u_pred = u + Dt*v + (0.5-beta)*Dt^2*a;
    
    % Equação de Equilíbrio para o passo i+1
    a_next = (0 - c*v_pred - k*u_pred) / m;
    
    % Correção das acelerações e velocidades
    v = v_pred + gamma*Dt*a_next;
    u = u_pred + beta*Dt^2*a_next;
    a = a_next;
    
    % Armazenar resultados
    u_history(i+1) = u;
    v_history(i+1) = v;
    a_history(i+1) = a;
end

% Plotar os resultados
figure;
subplot(3,1,1);
plot(time, u_history);
title('Deslocamento vs Tempo');
xlabel('Tempo (s)');
ylabel('Deslocamento (m)');

subplot(3,1,2);
plot(time, v_history);
title('Velocidade vs Tempo');
xlabel('Tempo (s)');
ylabel('Velocidade (m/s)');

subplot(3,1,3);
plot(time, a_history);
title('Aceleração vs Tempo');
xlabel('Tempo (s)');
ylabel('Aceleração (m/s^2)');
