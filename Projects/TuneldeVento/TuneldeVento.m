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
Erro = 0.00009;
deltaX = 0.006;
deltaY = deltaX;
iter = 1000;

% **********Elaboração do Domínio do problema************
% Usando até a primeira metade horizontal do domínio devido a simetria
[X, Y] = meshgrid(0:deltaX:(2*d + L)/2, 0:deltaY:H);

% Tamanho da malha
[Ny, Nx] = size(X); 

% Inicialização da matriz psi
psi = zeros(Ny, Nx);

% Inicialização=  da matriz psi
errorpsi = zeros(Ny, Nx);


% ======= Aplicação das condições de contorno
% Pontos dentro da semicircunferência no setor esquerdo
R = L/2;  % Raio da semicircunferência
Xc = d+L/2;  % Coordenada X do centro
Yc = h;  % Coordenada Y do centro
carro = ((X - Xc).^2 + (Y - Yc).^2 <= R^2) & (X <= Xc) & (Y >= Yc);

% Definir psi = 0 para esses pontos pertencentes ao carro
psi(carro) = 0;



%=============Determinação dos pontos de contorno irregular =====

% Inicializando a máscara externa
external_mask = false(size(X));

% Percorrendo a malha para identificar os pontos da máscara externa
for j = 2:Ny-1
    for i = 2:Nx-1
        if ~carro(j,i)
            % Verifica se o ponto imediatamente acima em X pertence a `carro`
            if carro(j, i+1)
                external_mask(j, i) = true;
            end
            
            % Verifica se o ponto imediatamente abaixo em Y pertence a `carro`
            if Y(j,i) > Yc
                if carro(j-1, i)
                    external_mask(j, i) = true;
                end
            else
                % Se Y < Yc, verifica se o ponto imediatamente acima é `carro`
                if carro(j+1, i)
                    external_mask(j, i) = true;
                end
            end
        end
    end
end

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
title('Pontos do contorno irregular');
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
    diff = 0;
    
    for k = 1:iter
        
        % Inicialização da variável de erro máximo dentro da iteração
        erromax = 0;
        % Avalia-se primeiro o valor do interior do domínio
        for j = 2:Ny-1
            for i = 2:Nx-1
                
                if carro(j,i)
                    continue; % Pulando atualização de pontos no carro
                end

                psi_old = psi(j,i);
                psi_new = (deltaY^2 * (psi(j,i+1) + psi(j,i-1)) + deltaX^2 * (psi(j+1,i) + psi(j-1,i))) / (2 * (deltaX^2 + deltaY^2));
                psi(j,i) = (1 - Lambda) * psi(j,i) + Lambda * psi_new;
                
                diff = abs((psi(j,i) - psi_old))/psi(j,i);
                errorpsi(j,i) = diff;

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
            
            diff = abs((psi(j,i) - psi_old))/psi(j,i);
            errorpsi(j,i) = diff;
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
            
            diff = abs((psi(j,i) - psi_old))/psi(j,i);
            errorpsi(j,i) = diff;

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
        
        diff = abs((psi(j,i) - psi_old))/psi(j,i);
        errorpsi(j,i) = diff;

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
        


        % =====Tratando pontos de contorno irregular=====
        for j = 2:Ny-1
            for i = 2:Nx
                if external_mask(j,i)
    
                    % Calcular os ângulos α e β com valores limitados ao intervalo [-1, 1]
                    alpha = acos(max(min(Y(j,i) / R, 1), -1));
                    beta = asin(max(min((Xc - X(j,i)) / R, 1), -1));
    
                    % Calcular a distância horizontal para frente (direita)
                    dist_x_forward = abs(Xc - R * sin(alpha) - X(j, i));
                    b = dist_x_forward / deltaX;
    
                    % Calcular a distância vertical para baixo
                    dist_y_down = abs(Y(j, i) - R * cos(beta) - Yc);
                    a = dist_y_down / deltaY;
    
                    if(Y(psi(j,i)>Yc)) %Para o caso de contorno irregular na parte superior
    
                        %Para o caso particular no topo da semicircunferencia
                        if(X(psi(j,i)) == Nx)
                            psi(j,i) = (2*(psi(j,i-1)) / deltaX^2 + (a * psi(j+1,i)) / (deltaY^2 * (a^2 + a))) / ...
                                (2 / deltaX^2 + (a * (a + 1)) / (deltaY^2 * (a^2 + a)));
                        
                        % Condição de ter apenas contorno irregular apenas em x
                        elseif (carro(j,i+1) && ~carro(j-1,i))
                            psi(j,i) = ((psi(j-1,i) + psi(j+1,i)) / deltaY^2 + (b * psi(j,i-1)) / (deltaX^2 * (b^2 + b))) / ...
                                (2 / deltaY^2 + (b * (b + 1)) / (deltaX^2 * (b^2 + b)));
                       
                        % Condição para contorno irregular apenas em y
                        elseif (~carro(j,i+1) && carro(j-1,i))
                            psi(j,i) = ((psi(j,i-1) + psi(j,i+1)) / deltaX^2 + (a * psi(j+1,i)) / (deltaY^2 * (a^2 + a))) / ...
                                (2 / deltaX^2 + (a * (a + 1)) / (deltaY^2 * (a^2 + a)));
                        
                        elseif(carro(j,i+1) && carro(j-1,i))
                        % Para caso de contorno irregular duplo, em x e y
                            psi(j,i) = (deltaY^2 * psi(j,i-1) + deltaX^2 * psi(j+1,i) + deltaY^2 * a * psi(j,i-1) + deltaX^2 * b * psi(j+1,i)) / ...
                                (deltaX^2 * a + deltaY^2 * a + deltaX^2 * b + deltaY^2 * b + deltaX^2 + deltaY^2 + deltaX^2 * a * b + deltaY^2 * a * b);
                        end
                    
                    % %Caso de contorno irregular abaixo da semicircunferência
                    % else
                    %     a = (Yc - Y(j,i)) / deltaY;
                    % 
                    %     psi(j,i) = ((psi(j,i-1) + psi(j,i+1)) / deltaX^2 + (2 * a * psi(j-1,i) * (a^2 + a)) /...
                    %         deltaY^2 )/(2 / deltaX^2 + (2 * a + 2) * (a^2 + a) / deltaY^2);
                    end

                end
            end
        end
        
        % Se for a última iteração, armazena o erro máximo
        if k == iter
            erromax_list = [erromax_list, erromax];
        end
    end

   
    % Verifica se houve alguma atualização significativa
    if contador > 0
        convergiu = false;
    end
    
    % Armazena o erro máximo da iteração atual e o número da iteração
    total_iters = total_iters + iter;
    iter_list = [iter_list, total_iters];

    % Exibe o número de atualizações significativas na iteração atual
    disp(['Iteração ', num2str(total_iters), ': Atualizações significativas = ', num2str(contador), ', Erro máximo = ', num2str(erromax)]);

    % Atualiza a animação do gráfico de psi
    figure(1);
    subplot(1,2,1);
    contourf(X, Y, psi, 25); % Atualiza as linhas de contorno (isobaras)
    colorbar;
    xlabel('X');
    ylabel('Y');
    axis image;
    title(['Evolução de Psi - Iteração ', num2str(total_iters)]);
    
     % Atualiza a animação do gráfico de erro psi
    subplot(1,2,2);
    contourf(X, Y, errorpsi, 50); % Atualiza as linhas de contorno do erro
    colorbar;
    xlabel('X');
    ylabel('Y');
    axis image;
    title(['Erro Psi - Iteração ', num2str(total_iters)]);

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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%===================================
% Calculo do campo de velocidades 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Inicializa matrizes para as componentes de velocidade
u = zeros(Ny, Nx);
v = zeros(Ny, Nx);

% Calcula as derivadas parciais usando diferenças finitas centrais
for j = 2:Ny-1
    for i = 2:Nx-1
        u(j,i) = (psi(j,i+1) - psi(j,i-1)) / (2 * deltaY); % \partial \psi / \partial y
        v(j,i) = -(psi(j+1,i) - psi(j-1,i)) / (2 * deltaX); % -\partial \psi / \partial x
    end
end

% Tratamento das bordas, em i = Nx
for j = 2:Ny-1
    i = Nx;
    u(j,i) = (psi(j,i-1) - psi(j,i-2)) / deltaY; % \partial \psi / \partial y, diferença retroativa
    v(j,i) = -(psi(j+1,i) - psi(j-1,i)) / (2 * deltaX); % -\partial \psi / \partial x, diferença central
end

% Tratamento das bordas, em j = 1 e j = Ny
for i = 2:Nx-1
    % Para j = 1 (borda inferior)
    j = 1;
    u(j,i) = (psi(j,i+1) - psi(j,i-1)) / (2 * deltaY); % \partial \psi / \partial y, diferença central
    v(j,i) = -(psi(j+1,i) - psi(j,i)) / deltaX; % -\partial \psi / \partial x, diferença avançada
    
    % Para j = Ny (borda superior)
    j = Ny;
    u(j,i) = (psi(j,i+1) - psi(j,i-1)) / (2 * deltaY); % \partial \psi / \partial y, diferença central
    v(j,i) = -(psi(j,i) - psi(j-1,i)) / deltaX; % -\partial \psi / \partial x, diferença retroativa
end

% Tratamento dos cantos
% Canto inferior direito (j = 1, i = Nx)
j = 1;
i = Nx;
u(j,i) = (psi(j,i-1) - psi(j,i-2)) / deltaY; % \partial \psi / \partial y, diferença retroativa
v(j,i) = -(psi(j+1,i) - psi(j,i)) / deltaX; % -\partial \psi / \partial x, diferença avançada

% Canto superior direito (j = Ny, i = Nx)
j = Ny;
i = Nx;
u(j,i) = (psi(j,i-1) - psi(j,i-2)) / deltaY; % \partial \psi / \partial y, diferença retroativa
v(j,i) = -(psi(j,i) - psi(j-1,i)) / deltaX; % -\partial \psi / \partial x, diferença retroativa

% Condições de contorno podem precisar ser aplicadas aqui para bordas
% (i.e., i = 1, i = Nx, j = 1, j = Ny)
% Isso pode ser feito usando diferenças finitas avançadas ou retroativas

% Reduzindo a resolução dos vetores para visualização
skip = 20; % Fator de redução da malha (ajuste conforme necessário)
X_reduced = X(1:skip:end, 1:skip:end);
Y_reduced = Y(1:skip:end, 1:skip:end);
u_reduced = u(1:skip:end, 1:skip:end);
v_reduced = v(1:skip:end, 1:skip:end);

% Espelhando as matrizes de velocidade e invertendo a direção das componentes
u_espelhado = [u_reduced, -fliplr(u_reduced)];
v_espelhado = [v_reduced, fliplr(v_reduced)];

% Espelhando as matrizes X e Y para a visualização correta
X_espelhado = [X_reduced, X_reduced + max(max(X_reduced))];
Y_espelhado = [Y_reduced, Y_reduced];  % Apenas duplicando Y, pois o espelhamento é horizontal

% Visualização do campo de velocidades com vetores maiores e dobrados
figure;
quiver(X_espelhado, Y_espelhado, u_espelhado, v_espelhado, 'AutoScaleFactor', 1.5); % Plota os vetores de velocidade
xlabel('X');
ylabel('Y');
title('Campo de Velocidades Completo');
axis equal;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%===================================
% Cálculo das pressoes 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



% Parâmetros conhecidos
p_ref = 101325; % Pressão de referência em Pa (por exemplo, pressão atmosférica)
gamma_ar = 1.4; % Razão de calores específicos para o ar
rho = 1.225; % Densidade do ar em kg/m^3

% Usando o ponto no topo do domínio (j = Ny) e no meio horizontal (i = Nx/2)
j_ref = Ny;
i_ref = Nx;

% Determinando a constante (cte) usando o ponto de referência
velocity_squared_ref = u(j_ref, i_ref)^2 + v(j_ref, i_ref)^2;
cte = (gamma_ar / (gamma_ar - 1)) * (p_ref / rho) + (velocity_squared_ref / 2);

% Inicializando a matriz de pressão
p = zeros(Ny, Nx);

% Calculando a pressão no domínio
for j = 1:Ny
    for i = 1:Nx
        velocity_squared = u(j,i)^2 + v(j,i)^2;
        p(j,i) = rho * (cte - (velocity_squared / 2)) * (gamma_ar - 1) / gamma_ar;
    end
end

% Espelhando a matriz de pressão
p_espelhado = [p, fliplr(p)];

% Criando a nova matriz X para o domínio espelhado
X_espelhado = [X, X + max(max(X))];

% Ajustando Y para a mesma largura de X_espelhado
Y_espelhado = [Y, Y];  % Apenas duplicando Y, pois o espelhamento é horizontal

% Visualização do campo de pressão dobrado
figure;
contourf(X_espelhado, Y_espelhado, p_espelhado, 25);
colorbar;
xlabel('X');
ylabel('Y');
title('Campo de Pressão Completo');
axis equal;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%=================================
% Calculo da pressão na carroceria
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

gamma_ar = 1.4; % Razão de calores específicos para o ar
rho = 1.225; % Densidade do ar em kg/m^3

% Definir o centro da semicircunferência
Xc = L/2; % Coordenada X do centro
Yc = h; % Coordenada Y do centro
R = L/2; % Raio da semicircunferência

% Inicializa a força de sustentação
F_lift = 0;

% Vetores para armazenar as coordenadas e as pressões na carroceria
x_carro = [];
y_carro = [];
p_carro = [];

% Percorre a malha para calcular a pressão nos pontos próximos à carroceria e armazenar os valores
for j = 1:Ny
    for i = 1:Nx
        if carro(j,i)
            % Calcula a normal à superfície em cada ponto (j, i)
            dx = X(j,i) - Xc;
            dy = Y(j,i) - Yc;
            n = [dx; dy] / norm([dx; dy]); % Normal unitária

            % Armazena as coordenadas e a pressão
            x_carro = [x_carro, X(j,i)];
            y_carro = [y_carro, Y(j,i)];
            p_carro = [p_carro, p(j,i)];
        end
    end
end

% Encontrar a menor pressão
[min_p, min_idx] = min(p_carro);
min_x = x_carro(min_idx);
min_y = y_carro(min_idx);

% Exibe a menor pressão
disp(['A menor pressão é: ', num2str(min_p), ' Pa, no ponto (', num2str(min_x), ', ', num2str(min_y), ')']);

% Calcula a pressão média ao longo da carroceria
p_media = mean(p_carro);

% Calcula a força de sustentação assumindo a pressão média ao longo da área total da semicircunferência
area_total = pi * (R^2) / 2; % Área da semicircunferência
F_lift = p_media * area_total;

% Exibe o resultado da força de sustentação
disp(['A força de sustentação é: ', num2str(F_lift), ' N']);

% Plotar a pressão ao longo da carroceria
figure;
scatter(x_carro, y_carro, 50, p_carro, 'filled');
colorbar;
hold on;
% Destacar o ponto de menor pressão
plot(min_x, min_y, 'ro', 'MarkerSize', 10, 'LineWidth', 2);
hold off;
xlabel('X');
ylabel('Y');
title('Pressão ao Longo da Carroceria');
axis equal;
