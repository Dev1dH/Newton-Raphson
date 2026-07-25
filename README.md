# Newton-Raphson para Sistemas Elétricos de Potência

O método de Newton-Raphson é a técnica numérica mais utilizada na resolução de equações não lineares no cálculo de fluxo de carga em sistemas elétricos de potência. A principal características que a coloca em destaque em relação ao método de Gauss-Seidel, é a maior velocidade de convergência. Enquanto no Gauss-Seidel a convergência é linear, no Newton-Raphson a convergência é quadrática. Além disso, o Newton-Raphson tem bom desempenho em sistemas maiores.

## Fundamentação matemática

O erro entre as potências especificadas e as potência calculadas nas barras $i$, é dado pelas seguintes equações;

$$\Delta P = P_{esp}-|V_i|\sum_{i=1}^{n} |V_j|(G_{ij}cos\theta_{ij}+B_{ij}sen\theta_{ij})$$

$$\Delta Q = Q_{esp} -|V_i|   \sum_{i=1}^{n} |V_j|(G_{ij}cos\theta_{ij}-B_{ij}sen\theta_{ij})$$

São agrupadas em um  único vetor de erros:

$$
\Delta F =
\begin{bmatrix}
\Delta P \\
\Delta Q
\end{bmatrix}
$$

Utilizando a expansão da série de Taylor, é possível transformar um sistema não linear em uma aproximação linear, mais fácil de ser resolvida. O método Newton-Raphson utiliza apenas a derivada de primeira ordem, 

$$
\begin{bmatrix}
\Delta P \\
\Delta Q
\end{bmatrix}

=
J
\begin{bmatrix}
\Delta \theta \\
\Delta V
\end{bmatrix}
$$

Em que $J$ é a matriz Jacobiana. A Jacobiana contém as derivadas parciais das potências em relação ás incógnitas. Ela é formada por 4 submatrizes H, N, M e L.

$$J=
\begin{bmatrix}
\mathbf{H} & \mathbf{N} \\
\mathbf{J} & \mathbf{L}
\end{bmatrix}
$$

Para os elementos na diagonal da Jacobiana, $j = i$ temos:

$$H_{ij} = |V_i|^2B_{ii}+Q_i$$
$$N_{ii} = -|V_i|^2G_{ii}-P_{i}$$
$$J_{ii} = |V_i|^2G_{ii}-P_{i}$$
$$L_{ii} = |V_i|^2B_{ii}-Q_i$$


Para os elementos fora da diagonal da Jacobiana, $j \neq i$ temos: 

$$H_{ij} = L_{ij} = -|V_i||V_j| (G_{ij}sin\theta_{ij}-B_{ij}cos\theta{ij})$$
$$N_{ij} = J_{ij} = -|V_i||V_j| (G_{ij}sin\theta_{ij}+B_{ij}cos\theta{ij})$$

o processo é repetido até que as correções se tornem tão pequenas que satisfaçam a precisão escolhida.

## Implementação
