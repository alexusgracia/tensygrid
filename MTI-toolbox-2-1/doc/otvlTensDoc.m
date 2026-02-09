%% otvlTens
%  OTVL tensor
%
%% Syntax
% |sys = otvlTens(FPhi, Fa, Fc, varargin)|
%

%% Description
% Use |sys = otvlTens(FPhi, Fa, Fc, varargin)| to create an OTVL tensor. 
%
% An OTVL tensor contains multiple OTVLs, which represent the structure of
% state equations of a multilinear model.
% 

%% Input arguments
% |FPhi| parameter of dimension total number of OTVL rows x 1
%
% |Fa| parameter of dimension 2 x number of variables
%
% |Fc| parameter of dimension number of Variables x 1
%
% |varargin| input of variable length containing an
% <otvlDoc.html |otvl|>-object for each
% state equation of the multilinear model

%% Outut arguments
% |sys| OTVL tensor containing state equations of MTI model

%% Multilinear models and OTVL
% OTVL impose structural restrictions on the contained multilinear
% equations, since they only allow certain multilinear combinations, see 
% also [1].
% The state equations of a multilinear model whose equation structure is 
% contained in the OTVL:
%
% <html>
% <table>
% <tr><th><var>x<sub>1</sub></var></th><th><var>x<sub>2</sub></var></th><th><var>x<sub>3</sub></var></th></tr>
% <tr><td>1</td><td>1</td><td>0</td></tr>
% <tr><td>1</td><td>-</td><td>1</td></tr>
% <tr><td>0</td><td>0</td><td>0</td></tr>
% </table>
% </html>
%
% <html>
% <table>
% <tr><th><var>x<sub>1</sub></var></th><th><var>x<sub>2</sub></var></th><th><var>x<sub>3</sub></var></th></tr>
% <tr><td>-</td><td>1</td><td>1</td></tr>
% <tr><td>-</td><td>0</td><td>0</td></tr>
% </table>
% </html>
%
% <html>
% <table>
% <tr><th><var>x<sub>1</sub></var></th><th><var>x<sub>2</sub></var></th><th><var>x<sub>3</sub></var></th></tr>
% <tr><td>-</td><td>-</td><td>0</td></tr>
% </table>
% </html>
%
% are represented by
%
% $$x_1' =  \phi_1 (x_1 + a_{11})(x_2+a_{12})(1-x_3-a_{23}) 
% + \phi_2(x_1 + a_{11})(x_3+a_{13}) 
% + \phi_3(1-x_1-a_{21})(1-x_2-a_{22})(1-x_3-a_{23})+c_1$$
%
% $$x_2' = \phi_4(x_2 + a_12)(x_3 + a_13) +
% \phi_5(1-x_2-a_{22})(1-x_3-a_{23}) + c_2$$
%
% $$x_3' = \phi_6(1-x_3-a_{23})$$ 
% 
% The dimensions of the parameters are
%
% Fphi: total number of TVL rows x 1
%
% Fa: 2 x number of Variables
%
% Fc: number of Variables x 1

%% Example create OTVL tensor
% An OTVL tensor contains equation structure and parameters of a
% multilinear model. It is created by using the constructor method of the
% class otvlTens.

%create TVL structures for state equations
myTvlStruct1= false(3,3,2);
myTvlStruct1(:,:,1) = [0 0 0; 0 1 0; 0 0 0]; %position of dont cares
myTvlStruct1(:,:,2) = [1 1 0; 1 0 1; 0 0 0]; %Boolean values
myTvlStruct2= false(2,3,2);
myTvlStruct2(:,:,1) = [1 0 0; 1 0 0]; %position of dont cares
myTvlStruct2(:,:,2) = [0 1 1;0 0 0]; %Boolean values
myTvlStruct3= false(1,3,2);
myTvlStruct3(:,:,1) = [1 1 0]; %position of dont cares
myTvlStruct3(:,:,2) = [0 0 0]; %Boolean values

%create  OTVLs
myOTvl1 = otvl(myTvlStruct1);
myOTvl2 = otvl(myTvlStruct2);
myOTvl3 = otvl(myTvlStruct3);

%create OTVL tensor
myC = zeros(3,1); 
myPhi = ones(6,1); 
myA = zeros(2,3);
oTens = otvlTens(myPhi, myA, myC, myOTvl1, myOTvl2, myOTvl3);

%% References
% [1] M. Engels, G. Lichtenberg, and S. Knorn. "An approach to 
% structured multilinear modeling with relaxed Boolean output functions", 
% in 22nd IFAC World Congress, Yokohama, Japan, 2023, pp.7920-7925

%% See Also
% <otvlDoc.html otvl>