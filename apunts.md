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