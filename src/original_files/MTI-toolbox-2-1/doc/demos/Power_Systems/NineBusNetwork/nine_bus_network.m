%% Write the bus connections.

% There are nv source nodes. 
% The node nv+1 is the ground node.

% Convention followed
%   List all the conections from ground to generating node first.
%   List all the conections from ground to load node second.
%   List all the conections from node to node in the desired order, third.

% Numer of blocks
nb = 3;
% Number of nodes
nv = 3*nb;
% Number of loads
nl = 2*nb;
% Number of lines
Lines = 2*nb + nb; % nb-1 if it is an open loop.

% Initialize bus matrix.
Bus = zeros(nv+nl+Lines,2);

% Node connections can be loaded manually, however in this network structure
% the busses follow the sequence [node-load, node, node-load, ...].

% Sources, indicate (ground node) -> (to generation bus)
for i=1:nv
    Bus(i,:) = [nv+1 i];
end

% Loads, indicate (ground node) -> (to laod)
% Use aux variables
x = 1; % In each block, loads are at nodes 1 and,
y = 3; % at node 3.

% Initialize an aux vector with the nl zero components.
aux = zeros(nl,1);
aux(1:2) = [1;3]; % the first two loads are in busses one and three.

% Exploit the sequence [node-load, node, node-load,.....] to alocate all
% other loads.
for i=2:nl/2 
    x = x + 3;
    y = x + 2;
    aux(2*i-1:2*i) = [x;y];
end

% Now, indicate (ground node) -> (to laod)
Bus(nv+1:nv+nl,:) = [(nv+1).*ones(nl,1) aux];

% Lines are conected sequentially, indicate (node j) -> (node k)
for i=1:Lines
    if (i+1)==nv+1 % when the loop is closed node 9 just goes to node 1
        Bus(i+nv+nl,:) = [i 1];
    else
        Bus(i+nv+nl,:) = [i i+1];
    end
end

% Create a cell for the voltage labels
node_labels = cell(1,nv+1);
for i=1:nv
    node_labels{i} = sprintf('%s%d','V',i);
end
node_labels{end} = 'V0';

% The weights are used to sort the connections in the chosen order.
weights = 1:1:length(Bus);

% Use matlab's graph functions to make life easier
G = digraph(Bus(:,1)',Bus(:,2)',weights,node_labels);

% Usally the order is no the same as the one we picked.
order = G.Edges.Weight;
% Sort from min to max, as in weights, save the permutation.
[~,indx]=sort(order);

% Create the incidence matrix, it's sparse, so make it full. 
Ah=incidence(G);
A = full(Ah);

% Delete the last row for ground conections, it is not needed anymore.
A = A(1:end-1,indx);

% A = [Agen Aloads Abar]
% Agen source to node incidence
% Aload load to node incidence
% Abar node to node incidence

plot(G,'Layout','force','EdgeLabel',G.Edges.Weight)
title('Network graph')
disp('Bus V0 is ground.')
 
%% Incidence Matrices
% Use Kronecker product to augment the incidence matrix,
% in this case use the I_2 matrix, since working with dq-coordinates.
Agen = kron(A(:,1:nv),eye(2));
Aload = kron(A(:,nv+1:nv+nl),eye(2));
Abar = kron(A(:,nv+nl+1:end),eye(2));
B = -Abar; % Matlab's sign convetion is flipped.



%% Symbolic variables
% Busses or nodes
Nodes = nv; %  for this example, all busses/nodes have sources.
I = symDQvars('I',Nodes);
v = symDQvars('v',Nodes);
% Transmission Lines
f = symDQvars('f',Lines);
% Voltage Sources
Igen=symDQvars('Igen',nv);
% Loads
Il=symDQvars('Il',nl);
%xload=symDQvars('xload',nl);

%% Transmission Lines Parameters
w0 = 2*pi*50;
% Blocks
R1 = 0.23;
X1 = 0.1;
L1 = X1/w0;
R2 = 0.35;
X2 = 0.58;
L2 = X2/w0;

% Block connections
R12 = 0.4;
L12 = 2e-3;

% Overall Inductance Matrix
Ldq = zeros(2,2,Lines);
Ldq(:,:,1) = L1.*eye(2);
Ldq(:,:,2) = L2.*eye(2);
Ldq(:,:,3) = L12.*eye(2);

% Same for the other blocks
for i=4:Lines
    Ldq(:,:,i) = Ldq(:,:,mod(i-1,3)+1);
end

% Build the matrix
L = [];
for i=1:Lines
    L = blkdiag(L,Ldq(:,:,i));
end

% Overall Resistence Matrix
Rdq = zeros(2,2,Lines);
Rdq(:,:,1) = R1.*eye(2);
Rdq(:,:,2) = R2.*eye(2);
Rdq(:,:,3) = R12.*eye(2);

% Same for the other blocks
for i=4:Lines
    Rdq(:,:,i) = Rdq(:,:,mod(i-1,3)+1);
end

% Build the matrix
R = [];
for i=1:Lines
    R = blkdiag(R,Rdq(:,:,i));
end

%% Loads

% Loads are found in nodes [1 3 4 6 7 9]

% Before disturbance

% Define resistive loads for all blocks, find the inverse.
rl = [25 20 25 20 25 20]'; % rl = [rl11 rl12 rl21 rl22 rl31 rl32]
iRl = 1./rl;

% Augment for two component dq vectors 
iRl = kron(diag(iRl),eye(2));

% Set iRl between Aload to get a conductance matrix where 
% the values are only in the nodes where the load shoud be.
G1 = Aload*iRl*Aload';

% After disturbance, same logic.
rl2 = [30 20 15 20 25 25]';
iRl2 = 1./rl2;
iRl2 = kron(diag(iRl2),eye(2));
G2 = Aload*iRl2*Aload';

% Load jump
 
% Mix of Conductance matrices
u1 = sym('u1','real');
G = u1*G1 + (1-u1)*G2;

% I belive just outputing these values in the command window is much more
% intuitive, take a look.

%% Virtual resistor to ground
% Look into section III-C in Mahmoud Kabalan' paper (...Large-Signal Nonlinear Analysis).
Rn = 1e3*eye(2*Nodes);

%% Equations in dq

% Rotation matrix to account for the dq-frame
Wlines = [];
for i=1:Lines
    Wlines = blkdiag(Wlines,wcom*[0 -1;1 0]);
end

% Network equations
eq1bl = L*diff(f,t)== B'*v - (R + Wlines*L)*f;
eq2bl = v == Rn*(Agen*Igen - G*v - B*f);

disp('Network equations in workspace')