function [F_U,F_phi,G_U,G_phi] = randCpn(n,p,m,r,boolean,normtype)
            %DRMSS  Generate random discrete-time multilinear state-space models.
            %
            % input parameter:
            %   n - number of states, defaults to random number between 1 and
            %        a positive random integer drawn from a standard normal distribution
            %       with standard deviation 10
            %   p - number of outputs, defaults to 1
            %   m - number of inputs, , defaults to 1
            %   r - rank of the resulting decomposed tensor, defaults to
            %          (n+m)/2
            %   sparsity - sparsity of the resulting decomposed  tensor,
            %              current defaults is 0.2 but it will be removed
            %              in the future 
            %   boolean - if set to true the structure matrix (U) will be logical
            %            instead of double
            %   normtype - normalization type of structure matrix (U),
            %           defaults to '1', the only normalization currently supported
            %    
            %
            %
            % Example: [F_U,F_phi,G_U,G_phi] = randCpn(500,1,2,0,'1')
            if nargin < 1
                n=max([1,round(abs(10*randn(1,1)))]);
            end
            if nargin < 2
                p=1;
            end
            if nargin < 3
                m=1;
            end
            if nargin < 4
                r=round(n+m/2);
            end
            
            if nargin < 5
                boolean = false;
            end
            if nargin < 6
                normtype = '1';
            end
            if normtype  ~= '1'
                error('Requested normtype not yet implemented')
            end
            
            nm=n+m;

           
            
            if nm>r
                F_U = sparse(randperm(nm),[randperm(r),randi(r,1,nm-r)],rand(nm,1));
            else %nm<r
                F_U=sparse([randperm(nm),randi(nm,1,r-nm)],randperm(r),rand(r,1));
            end
            
            if boolean
                F_U = F_U>0;
            end

            if n>r   % minimal sparsity implementation, as many non-zero enrties as states 
                F_phi = sparse(randperm(n),[randperm(r),randi(r,1,n-r)],rand(n,1));
            else % n=<r
                F_phi = sparse([randperm(n),randi(n,1,r-n)],randperm(r),rand(r,1));
            end

            G_U = sprand(n+m,r,0.2);
            if boolean
                G_U = G_U>0;
            end
            G_phi = sprand(p,r,0.7);


        end