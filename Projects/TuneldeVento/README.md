![Simulações de Escoamento e Transferência de Calor](https://github.com/DouglasOliveiraC/Tunel-de-Vento/blob/master/Projects/TuneldeVento/psi_completo_example.png)

# Projeto de Simulação de Escoamento e Transferência de Calor em um Veículo

## Descrição do Projeto
Este projeto tem como objetivo a simulação do escoamento de ar e a transferência de calor em um veículo utilizando métodos numéricos. O projeto é dividido em três partes principais:

1. **Simulação do escoamento e cálculo de pressão na superfície do veículo.**
2. **Cálculo da distribuição de temperatura e taxa de calor retirada do veículo.**
3. **Análise da força vertical resultante variando os parâmetros do veículo e a velocidade do escoamento.**

## Estrutura do Projeto

### Parte I: Simulação do Escoamento

- **a)** Plotar a função de corrente ψ do escoamento. [Status: ✅ Concluído]
- **b)** Plotar os vetores de velocidade absoluta do escoamento. [Status: ✅ Concluído]
- **c)** Plotar a pressão no domínio. [Status: ✅ Concluído]
- **d)** Plotar a pressão ao longo da carroceria, explicitando seu valor mínimo. [Status: ✅ Concluído]
- **e)** Calcular a força vertical resultante que atua no veículo. [Status: ✅ Concluído]

### Parte II: Distribuição de Temperatura e Taxa de Calor

- **a)** Calcular a distribuição de temperatura no ar (em °C). [Status: 🔄 Em andamento]
- **b)** Calcular a taxa de calor retirada do carro. [Status: 🔄 Em andamento]

### Parte III: Análise da Força Vertical Resultante

- **a)** Variar a altura do carro (h=0,10; 0,05; 0,20; 0,025 m). [Status: 🔄 Em andamento]
- **b)** Variar a velocidade do escoamento (V=75; 140 km/h). [Status: 🔄 Em andamento]
- **c)** Apresentar os resultados em gráficos com escala adequada (Força x h; Força x V). [Status: 🔄 Em andamento]
- **d)** Discutir brevemente os resultados obtidos: quais as correlações entre as forças e a variação dos parâmetros? [Status: 🔄 Em andamento]

## Condições de Contorno e Parâmetros Utilizados

- **Velocidade do vento**: V=100 km/h
- **Dimensões do carro**: h=0,15 m e L=3 m
- **Propriedades do ar**: 
  - ρ=1,25 kg/m³
  - γ_ar=1,4
  - k_ar=0,026 W/(m·K)
  - c_p_ar=1002 J/(kg·K)
- **Dimensões do domínio**: 
  - d=0,5 L
  - H=2 L
- **Temperaturas**: 
  - T_dentro=25°C
  - T_motor=80°C
  - T_fora=20°C

## Implementação

### Método de Sobrerrelaxação

- **Discretização da malha**: Δx=Δy
- **Parâmetro de sobrerrelaxação**: λ=1,85
- **Tolerância de convergência**: 0,01

## Como Usar

1. Clone este repositório:
   ```sh
   git clone https://github.com/DouglasOliveiraC/Tunel-de-Vento.git
2.  Navegue até o diretório do projeto:
    ```sh
    cd SEU_REPOSITORIO
3. Execute o script principal no MATLAB:


## Resultados
Função de corrente ψ do escoamento: ![Simulações de Escoamento e Transferência de Calor](https://github.com/DouglasOliveiraC/Tunel-de-Vento/blob/master/Projects/TuneldeVento/psi_completo_example.png)
Vetores de velocidade absoluta do escoamento: [[Imagem/screenshot]](https://github.com/DouglasOliveiraC/Tunel-de-Vento/blob/master/Projects/TuneldeVento/campo_vel_completo.png)
Pressão no domínio: [[Imagem/screenshot]](https://github.com/DouglasOliveiraC/Tunel-de-Vento/blob/master/Projects/TuneldeVento/campo_pressoes.png)
Pressão ao longo da carroceria: [[Imagem/screenshot]](https://github.com/DouglasOliveiraC/Tunel-de-Vento/blob/master/Projects/TuneldeVento/pressao_carroceria.png)

Sinta-se à vontade para contribuir com este projeto. Faça um fork do repositório, crie uma nova branch e envie um pull request.

## Licença
Este projeto está licenciado sob a MIT License - veja o arquivo LICENSE para detalhes.


### Conclusão

Este `README.md` está estruturado para fornecer uma visão geral clara do projeto, suas partes componentes e como executar o código. Certifique-se de ajustar os caminhos das imagens e qualquer outro detalhe específico do seu projeto.
