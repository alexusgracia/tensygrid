#Migració

## Equacions

Si tenim per exemple les equacions

$$
-\dot{x}_1 -x_1 + u = 0
\\
-\dot{x}_2 -x_2 + x_1x_2 = 0
$$

Això ens donarà les matrius:

$$
S = \begin{bmatrix}
    1 & 0 & 0 & 0 & 0 & 0 \\
    0 & 1 & 0 & 0 & 0 & 0 \\
    0 & 0 & 1 & 0 & 0 & 1 \\
    0 & 0 & 0 & 1 & 0 & 1 \\
    0 & 0 & 0 & 0 & 1 & 0
\end{bmatrix}

P = \begin{bmatrix}
    -1 & 0 & -1 & 0 & 1 & 0 \\
    0 & -1 & 0 & -1 & 0 & 1
\end{bmatrix}
$$

Tal que:

$$
\begin{array}{r|cccccc}
& \dot{x}_1 & \dot{x}_2 & x_1 & x_2 & u & x_1x_2 \\ \hline
\dot{x}_1 & 1 & 0 & 0 & 0 & 0 & 0 \\
\dot{x}_2 & 0 & 1 & 0 & 0 & 0 & 0 \\
x_1 & 0 & 0 & 1 & 0 & 0 & 1 \\
x_2 & 0 & 0 & 0 & 1 & 0 & 1 \\
u & 0 & 0 & 0 & 0 & 1 & 0 \\ \hline
P_1 & -1 & 0 & -1 & 0 & 1 & 0 \\
P_2 & 0 & -1 & 0 & -1 & 0 & 1
\end{array}
$$

Perquè la $x_1x_2$ és combinació de les altres : $x_1$ i $x_2$.

Si per exemple tinguessim:

$$
\dot{x}_1 = -x_1 + u + \textcolor{red}{2x_2}
\\
\dot{x}_2 = -x_2 + x_1x_2
$$

Aleshores fem:

$$
-\dot{x}_1 -x_1 + 2x_2 + u = 0
\\
-\dot{x}_2 - x_2 + x_2 x_2 = 0
$$

Això ens donarà les matrius:

$$
S = \begin{bmatrix}
    1 & 0 & 0 & 0 & 0 & 0 \\
    0 & 1 & 0 & 0 & 0 & 0 \\
    0 & 0 & 1 & 0 & 0 & 1 \\
    0 & 0 & 0 & 1 & 0 & 1 \\
    0 & 0 & 0 & 0 & 1 & 0
\end{bmatrix}

P = \begin{bmatrix}
    -1 & 0 & -1 & \textcolor{red}{2} & 1 & 0 \\
    0 & -1 & 0 & -1 & 0 & 1
\end{bmatrix}
$$

Tal que:

$$
\begin{array}{r|cccccc}
& \dot{x}_1 & \dot{x}_2 & x_1 & x_2 & u & x_1x_2 \\ \hline
\dot{x}_1 & 1 & 0 & 0 & 0 & 0 & 0 \\
\dot{x}_2 & 0 & 1 & 0 & 0 & 0 & 0 \\
x_1 & 0 & 0 & 1 & 0 & 0 & 1 \\
x_2 & 0 & 0 & 0 & 1 & 0 & 1 \\
u & 0 & 0 & 0 & 0 & 1 & 0 \\ \hline
P_1 & -1 & 0 & -1 & \textcolor{red}{2} & 1 & 0 \\
P_2 & 0 & -1 & 0 & -1 & 0 & 1
\end{array}
$$

Per tant, la idea és aconseguir passar de les equacions a la matriu.

## Todo

- [x] Lector d'equacions d'entrada
- [x] Lògica linearitzar
- [x] Generar classe tests
- [] Repassar tests
- [] Depurar codi
- [] Classe per a linearitzar equacions d'entrada
- [] Simplificar procés

Alguns exemples d'equacions

```python

'-dx1 - x1 + u',
'-dx2 - x2 + x1*x2',

'-dx1 - x1 + u+ 2x2'
'-dx2 - x2 + x1*x2'

'-3dx1 + 2x1 - 5u',


"2*x + y - 5",
"x - 3*y + z - 2",
"4*z - 8",

```

# Coses a afegir
# Explicació mètode
Un iMTI es defineix com un sistema d'equacions on igualem funcions a zero. D'altra manera, és una funció vectorial de vàries variables on totes s'anul·len pel conjunt solució:
$$
\langle \mathcal{F}| \mathcal{M}(\mathrm{x})\rangle=\mathrm{f}(\mathrm{x})=0.
$$

Podem expandir $\mathrm{f}(\mathrm{x})$ al voltant d'un punt d'operaxió $\mathrm{x_0}$ tal que
$$
\mathrm{f}(\mathrm{x})\approx\mathrm{f}(\mathrm{x}-\mathrm{x}_0)+\underbrace{\left(\nabla\mathrm{f}(\mathrm{x})\right)_{\mathrm{x_0}}}_{\mathrm{J}_\mathrm{f}(\mathrm{x})}(\mathrm{x}-\mathrm{x_0}),
$$

on el jacobià $\mathrm{J}_\mathrm{f}(\mathrm{x})$ queda definit com la matriu $\underbrace{r}_{\#\text{ eqs}}\times \underbrace{(n+m)}_{\#\textrm{ variables}}$

$$
\mathrm{J}_\mathrm{f}(\mathrm{x}):=\left(\dfrac{\partial f_i}{\partial x_j}\right),\ i=\{1,...,r\},\ j=\{1,...,n+m\}.
$$

Per definició, $\mathrm{f}(\mathrm{x})=0$. Respecte $\mathrm{f}(\mathrm{x}-\mathrm{x}_0)$, com estem expandint al voltant d'un punt d'operació, qualsevol $\textrm{x}':=\textrm{x}-\textrm{x}_0$ ha de formar part del conjunt solució, i.e.: $\mathrm{x}'\in\{\mathrm{x}'\text{ such that }\mathrm{f}(\mathrm{x'})=0\}$, pel que també s'anul·la. Llavors, ens queda
$$
\mathrm{J}_\mathrm{f}(\mathrm{x})|_{\mathrm{x}_0}(\mathrm{x}-\mathrm{x}_0)\approx0.
$$

Al que ens respecta a nosaltres, tindrem (abusant de notació usaré l'igual) redefinint $\mathrm{x}$ com el vector $(\dot{\mathrm{x}},\mathrm{x},\mathrm{u})$ amb dimensions $\mathbb{R}^{2n+m}$

$$
0=E(\mathrm{\dot{x}}-\underbrace{\mathrm{\dot{x}}_0}_{=0 \text{ (punt extrem)}})+A(\mathrm{x}-\mathrm{x}_0)+B(\mathrm{u}-\mathrm{u}_0)=E\Delta\mathrm{\dot{x}}+A\Delta\mathrm{x}+B\Delta\mathrm{u},
$$

que és

$$
-E\mathrm{\dot{x}}=A\Delta\mathrm{x}+B\Delta\mathrm{u}.
$$