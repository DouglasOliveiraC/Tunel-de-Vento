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
d = L*(1/2);
H = 2*L;
T_dentro = 25;
T_fora = 20;
T_motor = 80;

% Lista para armazenar os valores de erromax e iter
erromax_list = [];
iter_list = [];
total_iters = 0;

% Inicialização do objeto de vídeo para a animação
v = VideoWriter('psi_evolution.avi');
open(v);

% Indices de aplicação das diferenças finitas
Lambda = 1.85;
Erro = 0.001;
deltaX = 0.003;
deltaY = deltaX;
iter = 1000;

% **********Elaboração do Domínio do problema************
% Usando até a primeira metade horizontal do domínio devido a simetria
[X, Y] = meshgrid(0:deltaX:(2*d + L)/2, 0:deltaY:H);

% Tamanho da malha
[Ny, Nx] = size(X); 

% Inicialização da matriz psi
psi = zeros(Ny, Nx);

% ======= Aplicação das condições de contorno
% Pontos dentro da semicircunferência no setor esquerdo
R = L/2;  % Raio da semicircunferência
Xc = d+L/2;  % Coordenada X do centro
Yc = h;  % Coordenada Y do centro
carro = ((X - Xc).^2 + (Y - Yc).^2 <= R^2) & (X <= Xc) & (Y >= Yc);

% Definir psi = 0 para esses pontos pertencentes ao carro
psi(carro) = 0;

% Para a regiao de contorno irregular que tem um R = R + deltaX, acima de Yc
cont_Irreg1 = ((X - Xc).^2 + (Y - Yc).^2 <= (R+deltaX)^2) & (X <= Xc) & (Y >= Yc) & (~carro);

% Para a regiao abaixo da semicircunferencia
cont_Irreg2 = ((Y > h-deltaY)&(Y < h)&( X > d - deltaX)&(X <= d+L/2));

% Identificando pontos onde a semicircunferência corta a malha
boundary_points = ((X - Xc).^2 + (Y - Yc).^2 > R^2) & ((X - Xc).^2 + (Y - Yc).^2 < (R + deltaX)^2) & (X <= Xc) & (Y >= Yc);


% Encontrar pontos externos à semicircunferência com y < Yc deltaY e x >
% Xc -R - deltaX
% Encontrar pontos externos ao quarto de circunferência mas dentro do raio R + deltaX
% Encontrar pontos externos ao quarto de circunferência mas dentro do raio R + deltaX
% Encontrar pontos externos ao quarto de circunferência mas dentro do raio R + deltaX
external_mask = ((X - Xc).^2 + (Y - Yc).^2 > R^2) & ((X - Xc).^2 + (Y - Yc).^2 <= (R + sqrt(2)*deltaX/2)^2) & (Y >= Yc) & (X <= Xc);

% Visualização dos pontos externos e o quarto de circunferência
figure;
hold on;
% Plotando os pontos externos que estão dentro do raio R + deltaX
plot(X(external_mask), Y(external_mask), 'ro');

% Plotando o quarto de circunferência
theta = linspace(pi/2, pi, 100);  % ângulos para o quarto de circunferência de 90 a 180 graus
x_circ = Xc + R * cos(theta);
y_circ = Yc + R * sin(theta);
plot(x_circ, y_circ, 'b--', 'LineWidth', 1.5);

% Plotando o centro da circunferência
plot(Xc, Yc, 'gx', 'MarkerSize', 10, 'LineWidth', 2);

xlabel('X');
ylabel('Y');
title('Pontos externos dentro do raio R + \deltaX');
legend('Pontos externos', 'Quarto de circunferência', 'Centro da circunferência');
axis equal;

% Definir os limites do gráfico para evitar distorção
xlim([0, (2*d + L)/2]);
ylim([0, H]);
hold off;


%========================================= Iterações=======================  
convergiu = false;
while ~convergiu
    convergiu = true;
    contador = 0;
    erromax = 0;
    for k = 1:iter
        % Avalia-se primeiro o valor do interior do domínio
        for j = 2:Ny-1
            for i = 2:Nx-1
                psi_old = psi(j,i);
                psi_new = (deltaY^2 * (psi(j,i+1) + psi(j,i-1)) + deltaX^2 * (psi(j+1,i) + psi(j-1,i))) / (2 * (deltaX^2 + deltaY^2));
                psi(j,i) = (1 - Lambda) * psi(j,i) + Lambda * psi_new;
                
                diff = abs(psi(j,i) - psi_old);
                if diff > Erro
                    contador = contador + 1;
                end
                if diff > erromax
                    erromax = diff;
                end
            end
        end
        
        % Avalia-se os valores no teto ( j = Ny)
        j = Ny;
        for i = 2:Nx-1
            psi_old = psi(j,i);
            % Atualizando psi no teto
            psi_new = ((2*psi(j-1,i) + 2*deltaY*V)/deltaY^2 + (psi(j,i-1) + psi(j,i+1))/deltaX^2) / (2/deltaX^2 + 2/deltaY^2);
            psi(j,i) = (1 - Lambda) * psi(j,i) + Lambda * psi_new;
            
            diff = abs(psi(j,i) - psi_old);
            if diff > Erro
                contador = contador + 1;
            end
            if diff > erromax
                erromax = diff;
            end
        end
        
        % Avalia-se os valores na parede esquerda (i = 1)
        i = 1;
        for j = 2:Ny-1
            psi_old = psi(j,i);
            % Atualizando psi na parede esquerda
            psi_new = (deltaX^2 * psi(j-1,i) + deltaX^2 * psi(j+1,i) + 2 * deltaY^2 * psi(j,i+1)) / (2 * deltaX^2 + 2 * deltaY^2);
            psi(j,i) = (1 - Lambda) * psi(j,i) + Lambda * psi_new;
            
            diff = abs(psi(j,i) - psi_old);
            if diff > Erro
                contador = contador + 1;
            end
            if diff > erromax
                erromax = diff;
            end
        end
        
        % Avalia-se o ponto de quina esquerda superior (i = 2, j = Ny-1)
        j = Ny-1;
        i = 2;
        
        psi_old = psi(j,i);
        psi_new = (- V * deltaX^2 * deltaY + psi(j+1,i) * deltaX^2 + psi(j,i-1) * deltaY^2) / (deltaX^2 + deltaY^2);
        psi(j,i) = (1 - Lambda) * psi(j,i) + Lambda * psi_new;
        
        diff = abs(psi(j,i) - psi_old);
        if diff > Erro
            contador = contador + 1;
        end
        if diff > erromax
            erromax = diff;
        end
        
        % Correção das laterais
        psi(1:Ny,Nx) = psi(1:Ny,Nx-1);
        psi(1:Ny,1) = psi(1:Ny,2);
        
        % Reaplicação das condições de contorno
        psi(1,1:Nx) = 0;
        psi(carro) = 0;
    end

    % Tratando pontos vizinhos aos pontos cortados pela semicircunferência usando a equação da imagem
  
    % Tratando pontos vizinhos aos pontos cortados pela semicircunferência usando a equação da imagem
    for j = 2:Ny-1
        for i = 2:Nx-1
            if external_mask(j,i)
                % Encontrar o primeiro ponto anterior que não está na boundary
                % Calcular os ângulos α e β com valores limitados ao intervalo [-1, 1]
                alpha = acos(max(min(Y(j,i) / R, 1), -1));
                beta = asin(max(min((Xc - X(j,i)) / R, 1), -1));

                % Calcular a distância horizontal para frente (direita)
                dist_x_forward = abs(Xc - R * sin(alpha) - X(j, i));
                b = dist_x_forward / deltaX;

                % Calcular a distância vertical para baixo
                dist_y_down = abs(Y(j, i) - R * cos(beta) - Yc);
                a = dist_y_down / deltaY;
                % Condição de ter apenas contorno irregular apenas em x
                if (carro(j,i+1) && ~carro(j-1,i))
                    psi(j,i) = ((psi(j-1,i) + psi(j+1,i)) / deltaY^2 + (b * psi(j,i-1)) / (deltaX^2 * (b^2 + b))) / ...
                        (2 / deltaY^2 + (b * (b + 1)) / (deltaX^2 * (b^2 + b)));
                    % Condição para contorno irregular apenas em y
                elseif (~carro(j,i+1) && carro(j-1,i))
                    psi(j,i) = ((psi(j,i-1) + psi(j,i+1)) / deltaX^2 + (a * psi(j+1,i)) / (deltaY^2 * (a^2 + a))) / ...
                        (2 / deltaX^2 + (a * (a + 1)) / (deltaY^2 * (a^2 + a)));
                else
                    % Para caso de contorno irregular duplo
                    psi(j,i) = (deltaY^2 * psi(j,i-1) + deltaX^2 * psi(j+1,i) + deltaY^2 * a * psi(j,i-1) + deltaX^2 * b * psi(j+1,i)) / ...
                        (deltaX^2 * a + deltaY^2 * a + deltaX^2 * b + deltaY^2 * b + deltaX^2 + deltaY^2 + deltaX^2 * a * b + deltaY^2 * a * b);
                end
            end
        end
    end




    % Verifica se houve alguma atualização significativa
    if contador > 0
        convergiu = false;
    end
    
    % Atualiza o contador total de iterações
    total_iters = total_iters + iter;

    % Armazena o erro máximo da iteração atual e o número da iteração
    erromax_list = [erromax_list, erromax];
    iter_list = [iter_list, total_iters];

    % Exibe o número de atualizações significativas na iteração atual
    disp(['Iteração ', num2str(total_iters), ': Atualizações significativas = ', num2str(contador), ', Erro máximo = ', num2str(erromax)]);

    % Atualiza a animação do gráfico
    figure(1);
    contourf(X, Y, psi, 25); % Atualiza as linhas de contorno (isobaras)
    colorbar;
    xlabel('X');
    ylabel('Y');
    axis image;
    title(['Evolução de Psi - Iteração ', num2str(total_iters)]);
    drawnow;
    
    % Captura o quadro para o vídeo
    frame = getframe(gcf);
    writeVideo(v, frame);
end

% Fecha o objeto de vídeo
close(v);

% Plotando o gráfico de linha do tempo do erro
figure;
plot(iter_list, erromax_list, '-o');
xlabel('Número de Iterações');
ylabel('Erro Máximo');
title('Comportamento do Erro Máximo ao Longo das Iterações');
grid on;

% Adicionando uma linha horizontal em erro = 0.01
hold on;
yline(0.01, '--r', 'Erro = 0.01', 'LabelHorizontalAlignment', 'left', 'LabelVerticalAlignment', 'bottom');
hold off;

%Gráfico de valores de fluxo
figure;
contourf(X, Y, psi, 20); 
hold on;
contour(X, Y, psi, 20, 'k'); 
colorbar;
xlabel('X');
ylabel('Y');
title('Perfil dos Valores de Fluxo');
axis image;
hold off;


% Espelhando a matriz psi para a direita sem repetir a última coluna
psi_completo = [psi, fliplr(psi(:,1:end-1))];

% Plotando a matriz psi completa
figure;
contourf(psi_completo, 25);
colorbar;
xlabel('X');
ylabel('Y');
axis image;
title('Distribuição Completa de Psi');





