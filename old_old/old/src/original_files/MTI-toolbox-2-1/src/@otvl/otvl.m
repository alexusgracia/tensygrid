classdef otvl
    %OTVL Construct an orthogonal ternary vector list (otvl) containing 
    % structural information of a multilinear equation. The OTVL is given 
    % as a 3D array.
    %
    % obj = otvl(otvlArray) constructs an otvl object.
    %   Input:                                                                 
    %       otvlArray       3D otvl array       
    %
    %   Properties: 
    %       ttable          3D otvl array
    %                             
    % For detailed documentation, click <a href="matlab:open((which('otvlDoc.html')))">here</a>
    %
    % See also otvlTens
    
    % Marah Engels - 11/06/2024
    
    properties (Access = public)
        ttable {mustBeNumericOrLogical}
    end

    properties (Access = protected )
        otvlPosDontCare %position of don't care operators
        bvlCols %columns of the BVL corresponding to the TVL
        bvlRows %rows of the BVL corresponding to the TVL
        otvlRows %rows of the OTVL
        nodeshift {mustBeNumeric} = 0
    end


    methods 
        function nDontCare = get.otvlPosDontCare(obj)
            nDontCare = sum(obj.ttable(:,:,1) == 1, 2);
        end
        function cols = get.bvlCols(obj)
            cols = size(obj.ttable,2);
        end
        function bRows = get.bvlRows(obj)
            bRows = sum(2.^obj.otvlPosDontCare);
        end
        function tRows = get.otvlRows(obj)
            tRows = size(obj.ttable,1);
        end

    end

    methods (Access = public)
        function obj = otvl(otvlArray) %construct OTVL
            if nargin == 0
                obj.ttable = false(0,0,0);
            else
                obj.ttable = logical(otvlArray);
            end
        end
        %bvlObj = otvl2bvl(obj) %method transform otvl2bvl
        %tvlCpObj = otvl2cp(obj)  %method transform otvl2cp (not implemented)
    end
end