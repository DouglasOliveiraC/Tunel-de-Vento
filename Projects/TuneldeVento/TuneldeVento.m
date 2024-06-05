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


% Lista para armazenar os valores de erromax
erromax_list = zeros(1, 1000);

% Inicialização do objeto de vídeo para a animação
v = VideoWriter('psi_evolution.avi');
open(v);

% Indices de aplicação das diferenças finitas

Lambda = 1.85;
Erro = 0.04;
deltaX = 0.007;
deltaY = deltaX;


% **********Elaboração do Domínio do problema************

% Usando até a primeira metade horizontal do domínio devido a simetria

[X, Y] = meshgrid(  0:deltaX:(2*d + L)/2, 0:deltaY:H);

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
%psi(cont_Irreg1) = 10;

% Para a regiao abaixo da semicircunferencia

cont_Irreg2 = ((Y > h-deltaY)&(Y < h)&( X > d - deltaX)&(X <= d+L/2));
%psi(cont_Irreg2) = 10;


%========================================= Iterações  
convergiu = false;
while ~convergiu
    convergiu = true;
    contador = 0;
    erromax =0;
    for iter = 1: 100
           
        
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
        %=======================
        % Correção das laterais
        %=======================
        psi(1:Ny,Nx) = psi(1:Ny,Nx-1);
        psi(1:Ny,1) = psi(1:Ny,2);
           
    
        % ======= Reaplicação das condições de contorno
        
        % Condição inferior em y, psi = 0
        psi(1,1:Nx) = 0;
       
        % Definir psi = 0 para esses pontos pertencentes ao carro
        psi(carro) = 0;
    end

    % Verifica se houve alguma atualização significativa
    if contador > 0
        convergiu = false;
    end
    
     % Armazena o erro máximo da iteração atual
    erromax_list = [erromax_list, erromax];


    % Exibe o número de atualizações significativas na iteração atual
    disp(['Iteração ', num2str(iter), ': Atualizações significativas = ', num2str(contador), ', Erro máximo = ', num2str(erromax)]);

    % Atualiza a animação do gráfico
    figure(1);
    contourf(X, Y, psi, 25); % Atualiza as linhas de contorno (isobaras)
    colorbar;
    xlabel('X');
    ylabel('Y');
    axis image;
    title(['Evolução de Psi - Iteração ', num2str(iter)]);
    drawnow;
    
    % Captura o quadro para o vídeo
    frame = getframe(gcf);
    writeVideo(v, frame);
end

% Fecha o objeto de vídeo
close(v);

% Plotando os dados
figure;
imagesc(0:deltaX:(2*d + L)/2,0:deltaY:H, psi);
colorbar;
xlabel('X');
ylabel('Y');
axis image;
title('Distribuição de Psi com Semicircunferência');

% Ajustar a direção do eixo Y para corresponder ao layout cartesiano usual
set(gca, 'YDir', 'normal');

% Definir os limites dos eixos para corresponder ao domínio real
axis([0 (2*d + L)/2 0 H]);

% Alinhando o zero dos eixos
ax = gca;  % Obtém o handle do eixo atual
ax.XAxisLocation = 'origin';  % Move o X-axis para a origem
ax.YAxisLocation = 'origin';  % Move o Y-axis para a origem


% Plotando as linhas de contorno (isobaras)
figure;
contour(0:deltaX:(2*d + L)/2, 0:deltaY:H, psi, 40); % 20 níveis de contorno
colorbar;
xlabel('X');
ylabel('Y');
axis image;
title('Linhas de Contorno de Psi com Semicircunferência');

% Ajustar a direção do eixo Y para corresponder ao layout cartesiano usual
set(gca, 'YDir', 'normal');

% Definir os limites dos eixos para corresponder ao domínio real
axis([0 (2*d + L)/2 0 H]);

% Alinhando o zero dos eixos
ax = gca;  % Obtém o handle do eixo atual
ax.XAxisLocation = 'origin';  % Move o X-axis para a origem
ax.YAxisLocation = 'origin';  % Move o Y-axis para a origem

% Plotando o gráfico de linha do tempo do erro
figure(4);
plot(1:length(erromax_list), erromax_list, '-o');
xlabel('Número de Iterações');
ylabel('Erro Máximo');
title('Comportamento do Erro Máximo ao Longo das Iterações');
grid on;

% 
% % Plotando a superfície 3D
% figure;
% surf(0:deltaX:(2*d + L)/2, 0:deltaY:H, psi);
% colorbar;
% xlabel('X');
% ylabel('Y');
% zlabel('Psi');
% title('Superfície 3D da Distribuição de Psi com Semicircunferência');
% shading interp; % Interpolação de sombreamento para suavizar a superfície
% 
% % Ajustar a direção do eixo Y para corresponder ao layout cartesiano usual
% set(gca, 'YDir', 'normal');
% 
% % Definir os limites dos eixos para corresponder ao domínio real
% axis([0 (2*d + L)/2 0 H]);
% 

