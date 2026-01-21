classdef cpnTens < mtiTens
    % Enrico Uhlenberg, Christoph Kaufmann, Torben Warnecke, Carlos Cateriano Yáñez, Leona Schnelle - 12/06/2024
    properties (Access = public)
        U {mustBeNumericOrLogical}
        phi {mustBeNumericOrLogical}
    end
    
    properties (Access = protected )
        % Sparse helper variables - this will be implemented differently as
        % soon as we get to sparse/non-sparse differentiation.
        RowInd 
        ColInd
        DataVec
        SignVec
        rowsU
        colsU
        processXUfun % property that holds the logic calculating the next simulation step. See constructor of CPN1 for example
    end
    
    
    methods (Access = protected)
        [rows, cols] = getSizeU(sys)

        [nZRows, nZColumns, nZSignVector] = getSparseIndices(sys)
    end

    methods (Access = public)
         updateSparseIndicesAndSize(sys)
        
         F = cpn2Fulltens(sys)
        
         lsys = mLinSparse(msys,x0,u0)

         cpTens = cpn2Cp(obj)

         lin = cpn2Lin(obj)
    end

    methods (Static, Access = public) 
            [U,phi,rank] = mat2Cpn(M,normtype)
            [U,phi,rank,maxrank] = fulltens2Cpn(F,normtype,decompose, reduce)
            [F_U,F_phi,G_U,G_phi] = randCpn(n,p,m,rank,sparsity,boolean,normtype)
            [f0, f1, f2, f3, fn, Neq] = ktensor2Cpn(T)
            [F_U,F_phi,G_U,G_phi] = createRandomCpnTens(n,p,m,rank,sparsity,boolean,normtype)
        	cost=OptiMTI(W,x,u,tv,n,m,r,x0,ts,os,lsq,ntype)
            [FU,FPhi,cost]=timeseries2CpnAls(data,r,options)
            [FU,FPhi,cost]=timeseries2CpnAls_nonneg(data,r,options)
            [FU,FPhi,cost]=timeseries2CpnNonlin(data,r,options)

    end


end