%% Linearization and Perturbation Test for iMTI models
% Basado en la representación CPN para comparación con Python
clear; clc;

%% 1. Configuración del Modelo (Singular Case)
% Dimensiones
n = 3; % x1, y1, z1 (total de variables del sistema)
m = 1; % u1 (entrada)
r = 5; % Número de términos/productos en las ecuaciones
N = 2*n + m; % Vector v: [dx1; dx2; dx3; x1; y1; z1; u1]

% --- Definición de Matrices P y S ---
% Eq 1: 3*dx1*y1 + 6*(0.5 + 0.5*z1)*(0.66 - 0.33*u1)*x1
% Eq 2: y1 - 2*x1
% Eq 3: z1 + 0.5*y1 - 1

% S: Estructura de factores (Filas: dx1..3, x1..3, u1)
% Col 1: dx1 * y1
% Col 2: x1 * z1 * u1 (Término no lineal aproximado para la estructura)
% Col 3: y1
% Col 4: x1
% Col 5: z1
S = [1 0 0 0 0;  % dx1
     0 0 0 0 0;  % dx2
     0 0 0 0 0;  % dx3
     0 1 0 1 0;  % x1
     1 0 1 0 1;  % y1
     0 1 0 0 1;  % z1
     0 1 0 0 0]; % u1
S = sparse(S);

% P: Coeficientes globales
P = sparse([3  6  0  0  0;   % Eq 1
            0  0  1 -2  0;   % Eq 2
            0  0  0 0.5 1]);  % Eq 3

% --- Punto de Operación (v) ---
dx = [1; 0; 0]; 
x  = [2; 4; 5]; % x1, y1, z1
u  = [3]; 
v = [dx; x; u];

%% 2. Linealización Analítica (Algoritmo Scalable)
X = S.*v - abs(S);
izero = sparse(S.*v == abs(S) & S);
imone = sparse(X == -1);
ctwo = sum(imone) > 1;
X = spfun(@(x) x+1, X);
X(izero | imone) = 1;
[~, coli, val] = find(X);
Y = accumarray(coli, val, [r 1], @prod)'; 
Y(ctwo) = 0;
X = spfun(@(x) 1./x, X);
F = S.*Y.*X;

% Matriz LTI Combinada
EABC_analitico = full(P * F');

%% 3. Prueba de Perturbación (Linealización Numérica)
% Calculamos f(v) y f(v + delta) para aproximar la derivada
h = 1e-7; % Perturbación pequeña
EABC_num = zeros(size(EABC_analitico));

for i = 1:N
    v_plus = v;
    v_plus(i) = v_plus(i) + h;
    
    % Evaluar f(v_plus)
    Xp = 1 - abs(S) + S.*v_plus;
    Yp = prod(Xp, 1);
    f_plus = P * Yp';
    
    % Evaluar f(v)
    X0 = 1 - abs(S) + S.*v;
    Y0 = prod(X0, 1);
    f0 = P * Y0';
    
    % Diferencia finita (Jacobiano)
    EABC_num(:, i) = (f_plus - f0) / h;
end

%% 4. Comparación de Resultados
disp('--- MATRIZ EABC ANALÍTICA (Python/CPN) ---')
disp(EABC_analitico)

disp('--- MATRIZ EABC NUMÉRICA (Perturbación) ---')
disp(EABC_num)

diff_max = max(max(abs(EABC_analitico - EABC_num)));
fprintf('Error máximo entre métodos: %e\n', diff_max);

%% 5. Extracción y Estabilidad
E = -EABC_analitico(:, 1:n);
A = EABC_analitico(:, n+1:2*n);

disp('--- Análisis de Estabilidad ---')
lambda = eig(A, E);
lambda_finito = lambda(isfinite(lambda));
disp('Autovalores del sistema (finitos):')
disp(lambda_finito)