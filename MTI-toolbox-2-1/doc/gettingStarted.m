%% Getting Started with the MTI Toolbox
% Tools for Multilinear Modeling, Simulation, and Analysis  
%  
%% Theory
% 
% Multilinear functions are nonlinear, a subclass of polynomials and a
% superclass of linear as well as Boolean functions. 
% Using these as right hand sides of nonlinear state space models, 
% they can represent deterministic multi-valued or binary dynamics  
% in the same way that linear time-invariant (LTI) models can be 
% represented as matrices of their parameters 
% and linear algebra is the mathematical theory behind them, 
% the parameters of multilinear time-invariant (MTI) models 
% can be represented as tensors 
% and multilinear algebra holds the mathematical background, which was first introduced in [1].
% 
% Early results showed that the dynamics of hybrid systems consisting of 
% continuous-valued as well as a discrete-valued parts can be represented 
% in an adquate way by MTI models and reduction methods can be applied [2]. 
%
% This allows the use of modern decomposition methods for the parameter tensors, 
% like canonical polyadic (CP), which are sums of outer products.  
% Because of the special size $n \times 2 \times 2 \times \cdots \times 2$ of the tensors,
% normalization of the factors leads to further reduced representations (CPN), 
% first introduced in [3] and used throughout this toolbox.
%
% Moreover, new methods could be beneficially used in engineering applications like 
% multilinearization of nonlinear models, first introduced in [4]. 
% This technique allows the approximation of nonlinear dynanmic behaviour
% in within a operating region and not only a point as linearization does.
% 
% The class of explict MTI models is not closed under composition, i.e. 
% a connection of two MTI models might be no longer multilinear, but polynomial. 
% In contrast, implicit MTI models, also called descriptor models are a closed class, [5]. 
% Many engineering tasks already profit from implicit modeling, 
% e.g. causality is not necessary to predefine, or conservation laws can be enforced.
% This comes with an increased complexity of the methods but still, 
% a MTI descriptor model can be described by a single parameter tensor.  
%   
% Tensor decompositions show - starting in the decade 2010-2020 - their ability 
% to break the curse of dimensionality for many relevant applications. 
% MTI modeling enables high performance nonlinear modeling for large scale hybrid systems, 
% for simulation, analysis and design.   
% Don't confuse: Multilinear are not multi-linear models, the latter are multiple linear
% models, e.g. for every operating point an individual LTI model.  
% These could also be represented as tensors and within the MTI toolbox, 
% but suffer from the curse of dimensionality. 
%% Syntax  
% Most of the MTI-Toolbox commands have the same or a similar syntax 
% as the Control System Toolbox, some only differ by the prefix |m|, 
% many are overloaded methods and share the exact syntax. 
% Many basic commands are available, e.g. 
%
% * defining a multilinear state space model |msys = <output/mssDoc.html  mss>(F,G,Ts)| by  
% transition and output tensor |F| and |G|, 
% which works like |lsys = ss(A,B,C,D,Ts)|.
% * simulating the model by |[y,tOut,x] = <output/msimDoc.html  msim>(msys,u,t,x0)|
% similar to |lsim|.
% * defining an implicit descriptor state space model |dmsys = <output/dmssDoc.html  dmss>(H)| by the tensor |H|, 
% comparable to |dlsys = dss(A,B,C,D,E,Ts)|.
% * simulating this model by |[y,tOut,x] = <output/dmsimDoc.html dmsim>(dmsys,u,t,x0)|
% similar to |lsim|.
% * state transformation |<output/mss2mssDoc.html  mss2mss>| similar to |ss2ss|
% * discretization w.r.t time by |<output/dmssc2dDoc.html  c2d>| or the inverse by |<output/d2cDoc.html d2c>|
% * linearization by |<output/linearizeDoc.html linearize>|, which have fast and sparse implementations 
% * grey box parameter estimation by |<output/mlgreyestDoc.html mlgreyest>|
% 
% Moreover, there are features exclusive for MTI models, e.g.
%
% * representation of CPN factors of parameter tensors as |<output/CPN1Doc.html CPN1>| for a 1-norm
% * multilinearization by |<output/mlinearizeDoc.html mlinearize>|, i.e. a multilinear approximation of 
% a SIMULINK block 
% * converting an LTI to a MTI model by |<output/ss2mssDoc.html ss2mss>| which allows 
% factorized representations.  
%
%% Applications
% Application of MTI models started with HVAC systems [2].
% Within projects, even real-time applicabilty has been shown. 
% New results show the applicability to power networks. 
% This opens the door to large scale multi-energy systems. 
% More information can be found on <http://mti.systems mti.systems>
%
%
%% References
%
% [1] G. Lichtenberg (2012): Hybrid Tensor Systems, Habilitation, TU
% Hamburg.
%
% [2] G. Pangalos, A. Eichler, G. Lichtenberg (2015):  
% Hybrid Multilinear Modeling and Applications,
% https://doi.org/10.1007/978-3-319-11457-6_5 
%
% [3] K. Kruppa, G. Pangalos,  G. Lichtenberg (2014):  
% Multilinear approximation of nonlinear state space models,
% https://doi.org/10.3182/20140824-6-za-1003.00455
% 
% [4] L. Schnelle, G. Lichtenberg, C. Warnecke (2022): 
% Using Low-rank Multilinear Parameter Identification for Anomaly Detection of Building Systems,
% https://doi.org/10.1016/j.ifacol.2022.07.173
%
% [5] G. Lichtenberg, G. Pangalos, C. Cateriano Yáñez, A. Luxa, N. Jöres, L. Schnelle, C. Kaufmann (2022):
% Implicit multilinear modeling: An introduction with application to energy systems
% https://doi.org/10.1515/auto-2021-0133
