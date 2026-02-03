%% mss
% Multilinear state-space model (explicit)

%% Description
% Use |mss| to create an explicit real-valued multilinear time-invariant 
% state-space model.
%
% A multilinear time-invariant (MTI) state-space model is a mathematical
% representation of a physical system as a set of input, output, and state
% variables related by first-order differential equations. The |mss| model 
% object can represent SISO or MIMO MTI state-space 
% models in continuous time or discrete time.
%
% A continuous-time explict MTI state-space model can be
% represented by contracted tensor products:
%
% $$
% \dot{\mathbf{x}} = \left\langle \mathcal{F},
% \mathcal{M}\left(\mathbf{x},\mathbf{u}\right)\right\rangle ,
% $$
% 
% $$
% \mathbf{y} = \left\langle \mathcal{G},
% \mathcal{M}\left(\mathbf{x},\mathbf{u}\right)\right\rangle
% $$
%
% Here, $\mathbf{x}\in\mathrm{R}^n$, $\mathbf{u}\in\mathrm{R}^m$, and 
% $\mathbf{y}\in\mathbf{R}^p$ represent the states, inputs, and outputs of 
% sizes $n\in\mathbf{Z}_{>=0}$, $m\in\mathbf{Z}_{>=0}$, and $p\in\mathbf{Z}_{>=0}$
% respectively, while $\mathrm{F}\in\mathbf{R}^{\times^{(n+m)}2\times n}$ and
% $\mathrm{G}\in\mathbf{R}^{\times^{(n+m)}2\times p}$ are the state-space 
% parameter tensors. While $\mathrm{F}$ and $\mathrm{G}$ can have an
% arbitrary rank, monomial tensors are always rank 1, e.g., 
% $\mathrm{M}\left(\mathbf{x},\mathbf{u}\right)$ canonical
% polyadic decomposition is:
%
% $$
% \mathrm{M}\left(\mathbf{x},\mathbf{u}\right) = \left[
% \left(\begin{array}{c}1\\x_1\end{array}\right),
% \cdots,
% \left(\begin{array}{c}1\\x_n\end{array}\right),
% \left(\begin{array}{c}1\\u_1\end{array}\right),
% \cdots,
% \left(\begin{array}{c}1\\u_m\end{array}\right)
% \right],
% $$
%
% where each factor has a single column, i.e., rank 1.
% 
% Alternatively, |mss| also supports MTI state-space
% systems in matrix format:
% 
% $$
% \dot{\mathbf{x}} = \mathbf{F}
% \mathbf{m}\left(\mathbf{x},\mathbf{u}\right)
% $$
%
% $$
% \mathbf{y} = \mathbf{G}
% \mathbf{m}\left(\mathbf{x},\mathbf{u}\right)
% $$
%
% Here, the parameter matrices are given as $\mathrm{F}\in\mathbf{R}^{n\times 2^{(n+m)}}$ and
% $\mathrm{G}\in\mathbf{R}^{p\times 2^{(n+m)}}$, and the monomial vector is
% 
% $$
% \mathbf{m}\left(\mathbf{x},\mathbf{u}\right) =
% \left(\begin{array}{c}1\\u_n\end{array}\right)\otimes
% \cdots\otimes
% \left(\begin{array}{c}1\\u_1\end{array}\right)\otimes
% \left(\begin{array}{c}1\\x_n\end{array}\right)\otimes
% \cdots\otimes
% \left(\begin{array}{c}1\\x_1\end{array}\right),
% $$
% 
% where $\otimes$ is the Kronecker product.
%
% For efficieny reasons, as seen from the size scaling of each format, 
% |mss| handles all systems internally in the tensor format.
%
% You can create a MTI state-space model object by 
% either specifying the system and output parameters or by converting a 
% model of another type (such as a linear time-invariant state-space model 
% |ss|) to  MTI state-space form. You can use an
% |mss| model object to:
%
% * Simulate an MTI system
% * Combine with other MTI models to represent a more complex system
% * Linearize an MTI system
% 
% For further reference see [1].
% 
%% Creation
% <html><h3>Syntax</h3></html> 
% 
%   |sys = mss(Ftens,Gtens)|
%   |sys = mss(Fmat,Gmat)|
%   |sys = mss(Ftens,Gtens,ts)|
%   |sys = mss(Fmat,Gmat,ts)|
% 
% <html><h2>Description</h2></html>
% 
% |sys = mss(Ftens,Gtens)| creates a continuous-time mutlilinear
% time-invariant state-space model object of the following form:
% 
% $$
% \dot{\mathbf{x}} = \left\langle \mathrm{F}_\mathrm{tens},
% \mathrm{M}\left(\mathbf{x},\mathbf{u}\right)\right\rangle ,
% $$
%
% $$
% \mathbf{y} = \left\langle \mathrm{G}_\mathrm{tens},
% \mathrm{M}\left(\mathbf{x},\mathbf{u}\right)\right\rangle .
% $$
% 
% For instance, consider a plant with $n$ states, $m$ inputs, and $p$
% outputs. The state-space tensors are:
% 
% * $\mathrm{F}_\mathrm{tens}$ is a 2 by $(n+m-1)$ times 2 by $n$ tensor
% * $\mathrm{G}_\mathrm{tens}$ is a 2 by $(n+m-1)$ times 2 by $p$ tensor
%
% Note that |Gtens| is an optional input.
% 
% |sys = mss(Fmat,Gmat)| creates a continuous-time mutlilinear
% time-invariant state-space model object of the following form:
% 
% $$
% \dot{\mathbf{x}} = \mathbf{F}_\mathrm{mat}
% \mathbf{m}\left(\mathbf{x},\mathbf{u}\right) ,
% $$
%
% $$
% \mathbf{y} = \mathbf{G}_\mathrm{mat}
% \mathbf{m}\left(\mathbf{x},\mathbf{u}\right) .
% $$
% 
% For instance, consider a plant with $n$ states, $m$ inputs, and $p$
% outputs. The state-space tensors are:
% 
% * $\mathbf{F}_\mathrm{mat}$ is a $n \times 2^{(n+m)}$ matrix
% * $\mathbf{G}_\mathrm{mat}$ is a $p \times 2^{(n+m)}$ matrix
% 
% Note that |Gmat| is an optional input with a 
% default output $\mathbf{y}=\mathbf{x}$
% 
% |sys = mss(Ftens,Gtens,ts)| creates a discrete-time mutlilinear
% time-invariant state-space model object of the following form with the 
% sample time |ts| (in seconds):
% 
% $$
% \mathbf{x}\left[k+1\right] = \left\langle \mathrm{F}_\mathrm{tens},
% \mathrm{M}\left(\mathbf{x}\left[k\right],\mathbf{u}\left[k\right]\right)\right\rangle
% $$
%
% $$
% \mathbf{y}\left[k\right] = \left\langle \mathrm{G}_\mathrm{tens},
% \mathrm{M}\left(\mathbf{x}\left[k\right],\mathbf{u}\left[k\right]\right)\right\rangle
% $$
% 
% with discrete time variable $k$. For instance, consider a plant with $n$ 
% states, $m$ inputs, and $p$ outputs. The state-space tensors are:
% 
% * $\mathrm{F}_\mathrm{tens}$ is a 2 by $(n+m-1)$ times 2 by $n$ tensor
% * $\mathrm{G}_\mathrm{tens}$ is a 2 by $(n+m-1)$ times 2 by $p$ tensor
%
% Note that |Gtens| is an optional input.
%
% |sys = mss(Fmat,Gmat,ts)| creates a discrete-time mutlilinear
% time-invariant state-space model object of the following form with the 
% sample time |ts| (in seconds):
% 
% $$
% \mathbf{x}\left[k+1\right] = \mathbf{F}_\mathrm{mat}
% \mathbf{m}\left(\mathbf{x}\left[k\right],\mathbf{u}\left[k\right]\right) ,
% $$
%
% $$
% \mathbf{y}\left[k\right] = \mathbf{G}_\mathrm{mat}
% \mathbf{m}\left(\mathbf{x}\left[k\right],\mathbf{u}\left[k\right]\right)
% ,
% $$
% 
% with discrete time variable $k$. For instance, consider a plant with $n$ 
% states, $m$ inputs, and $p$ outputs. The state-space tensors are:
% 
% * $\mathbf{F}_\mathrm{mat}$ is a $n$ by 2 to the power of $(n+m)$ matrix
% * $\mathbf{G}_\mathrm{mat}$ is a $p$ by 2 to the power of $(n+m)$ matrix
% 
% Note that |Gmat| is an optional input.
%
% <html><h2>Input Arguments</h2></html>
%
% <html><h3>Ftens - System tensor</h3></html> 
% 
% System tensor, specified as a 2 by $(n+m-1)$ times 2 by $n$ tensor,
% where, $n$ is the number of states and $m$ is the number of inputs. This
% input sets the value of property |Ftens|.
%
% |Ftens| accepts the following tensor formats:
%
% * multidimensional array 
% * <CPN1Doc.html |CPN1|>
% * <otvlDoc.html |otvl|>
% * |ktensor| from <https://www.tensortoolbox.org/ Tensor Toolbox>
%
% <html><h3>Gtens - Output tensor</h3></html> 
% 
% Output tensor, specified as a 2 by $(n+m-1)$ times 2 by $p$ tensor,
% where, $n$ is the number of states, $m$ is the number of inputs, and $p$ 
% is the number of outputs. This input sets the value of property |Gtens|.
%
% |Gtens| accepts the following tensor formats:
%
% * multidimensional array 
% * <CPN1Doc.html |CPN1|>
% * |ktensor| from <https://www.tensortoolbox.org/ Tensor Toolbox>
%
% <html><h3>Fmat - System matrix</h3></html> 
% 
% System matrix, specified as a $n$ by 2 to the power of $(n+m)$ matrix,
% where, $n$ is the number of states and $m$ is the number of inputs. This
% input sets the value of property |Fmat|.
%
% <html><h3>Gmat - Output matrix</h3></html> 
% 
% Output matrix, specified as a $p$ by 2 to the power of $(n+m)$ matrix,
% where, $n$ is the number of states, $m$ is the number of inputs, and $p$ 
% is the number of outputs. This input sets the value of property |Gmat|.
%
% <html><h3>ts - Sample time</h3></html> 
% 
% Sample time in seconds, specified as a scalar.
%
% <html><h2>Output Arguments</h2></html>
%
% <html><h3>sys - Output system model</h3></html> 
%
% Output system model returned as:
%
% * A multilinear time-invariant state-space (|mss|) model object, when the
% inputs $\mathrm{F}$ and $\mathrm{G}$ are numeric tensors (or matrices), 
% or when converting from another model object type.
%
%% Properties
% 
% <html><h3>F - System tensor</h3></html> 
% 
% System tensor, specified as a 2 by $(n+m-1)$ times 2 by $n$ tensor,
% where, $n$ is the number of states and $m$ is the number of inputs. This 
% property is internally stored in <CPN1Doc.html |CPN1|>
% format.
%
% <html><h3>G - Output tensor</h3></html> 
% 
% Output tensor, specified as a 2 by $(n+m-1)$ times 2 by $p$ tensor,
% where, $n$ is the number of states, $m$ is the number of inputs, and $p$ 
% is the number of outputs. This property is internally stored in <CPN1Doc.html |CPN1|>
% format.
%
% <html><h3>n - Number of states</h3></html> 
% 
% Number of states, specified as an integer $n\geq 0$.
%
% <html><h3>m - Number of inputs</h3></html> 
% 
% Number of inputs, specified as an integer $m\geq 0$.
%
% <html><h3>p - Number of outputs</h3></html> 
% 
% Number of outputs, specified as an integer $p\geq 0$.
%
% <html><h3>ts - Sample time</h3></html> 
% 
% Sample time in seconds, specified as a double $ts\geq 0$.
%
%% Object Functions
% The following lists contain a representative subset of the functions you
% can use with or in context of |mss| model objects.
%
% <html><h3>Multilinear Simulation</h3></html> 
% 
% <html>
% <table border=1>
% <tr><td><tt><a href="msimDoc.html">msim</a></tt></td><td>Plot simulated time response of dynamic 
% system to arbitrary inputs; simulated response data</td></tr>
% <tr><td><tt><a href="linearizeDoc.html">linearize</a></tt></td><td>Linearize an explicit multilinear model around an operating point</td></tr>
% <tr><td><tt><a href="mlinearizeDoc.html">mlinearize</a></tt></td><td>Obtains a continuous or discrete multilinear time-invariant state-space
% model, within an operating domain, of a Simulink model.</td></tr>
% <tr><td><tt><a href="timeseries2CpnAlsDoc.html">mlgreyest</a></tt></td><td>GreyBox parameter identification of MTI model from input-output data in time domain</td></tr>
% <tr><td><tt><a href="rmssDoc.html">rmss</a></tt></td><td>Generate a random continuous multilinear states-space model</td></tr>
% <tr><td><tt><a href="c2dDoc.html">c2d</a></tt></td><td>Transform continous time to discrete time mti system</td></tr>
% <tr><td><tt><a href="MEXaccelerateDoc.html">c2d</a></tt></td><td>Enable MEX acceleration (experimental)</td></tr>
% </table>
% </html>
%
%% Examples
%
% <html><h3>SISO Multilinear Time-Invariant State-Space Model</h3></html>
%
% Create a SISO MTI state-space model defined by the following system
% equations:
%
% $$
% \dot{x_1} = x_1 x_2 + 0.5 u
% $$
%
% $$
% \dot{x_2} = 2 x_1 u + 7
% $$
%
% $$
% y = 3 x_1 x_2 u +5
% $$
%
% The model has the following parameter matrices
%
% $$
% \mathbf{F} = \left(\begin{array}{cccccccc}
% 0 & 0 & 0 & 1 & 0.5 & 0 & 0 & 0 \\ 
% 7 & 0 & 0 & 0 & 0 & 2 & 0 & 0
% \end{array}\right),
% $$
%
% $$
% \mathbf{G} = \left(\begin{array}{cccccccc}
% 5 & 0 & 0 & 0 & 0 & 0 & 0 & 3
% \end{array}\right),
% $$
%
% and monomial vector
% 
% $$
% \mathbf{m}\left(x_1,x_2,u\right) = \left(\begin{array}{c}
% 1 \\ x_1 \\ x_2 \\ x_1 x_2 \\ u \\ x_1 u \\ x_2 u \\ x_1 x_2 u
% \end{array}\right).
% $$
%
% Specify the $\mathbf{F}$ and $\mathbf{G}$ matrices, and create the MTI
% state-space model.
%
F = [ 0  0  0  1  0.5 0  0  0 ;
      7  0  0  0  0   2  0  0];
G = [ 5  0  0  0  0   0  0  3 ;];
sys = mss(F, G)
%%
% The same system can be represented in
% <CPN1Doc.html |CPN1|> format:
%
% $$
% \mathbf{F}_\mathrm{U} = \left(\begin{array}{cccc}
% 1 & 0 & 1 & 0 \\ 
% 1 & 0 & 0 & 0 \\
% 0 & 1 & 1 & 0
% \end{array}\right),
% $$
%
% $$
% \mathbf{F}_\mathrm{phi} = \left(\begin{array}{cccc}
% 1 & 0.5 & 0 & 0 \\ 
% 0 & 0 & 2 & 7
% \end{array}\right),
% $$
%
% $$
% \mathbf{G}_\mathrm{U} = \left(\begin{array}{cc}
% 1 & 0 \\ 
% 1 & 0 \\
% 1 & 0
% \end{array}\right),
% $$
%
% $$
% \mathbf{G}_\mathrm{phi} = \left(\begin{array}{cc}
% 3 & 5
% \end{array}\right).
% $$
%
% Initialize the $\mathrm{F}$ and $\mathrm{G}$ tensors in <CPN1Doc.html |CPN1|>
% format , and create the MTI state-space model.
%
F_U = [1  0  1  0;
       1  0  0  0;
       0  1  1  0];
F_phi = [1  0.5 0  0;
         0  0   2  7];
F=CPN1(F_U,F_phi);
G_U = [1  0;
       1  0;
       1  0];
G_phi = [3  5];
G=CPN1(G_U,G_phi);
sys = mss(F, G)
%% References
%
% [1] Lichtenberg, Gerwald; Pangalos, Georg; Cateriano Yáñez, Carlos; Luxa,
% Aline; Jöres, Niklas; Schnelle, Leona; Kaufmann, Christoph (2022): 
% Implicit multilinear modeling. In at - Automatisierungstechnik 70 (1), 
% pp. 13–30. DOI: 10.1515/auto-2021-0133.
%
%% See Also
%
% msim