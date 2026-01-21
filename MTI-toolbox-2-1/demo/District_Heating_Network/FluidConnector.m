classdef FluidConnector < handle & matlab.mixin.SetGet
    properties
        from_object
        from_type
        to_object
        to_type
        connectorModel;
    end
    methods 
        function obj = FluidConnector(props)
            arguments
                props.?FluidConnector
            end
                set(obj,props)
                obj.connectorModel = obj.TotalconnectorModel();
        end

        function connectorModel = TotalconnectorModel(obj)
            H1 = hyCPN1();
            H1.F.algebraic = sparse([1:6],[1:6],1,6,6);
            H1.phi.equality = [1,0,0,1,0,0; 0,1,0,0,-1,0; 0,0,1,0,0,-1];
            connectorModel = dmss(H1,0);
            connectorModel.algebraicName = [obj.from_object + ["_m_", "_t_", "_p_"] + obj.from_type,...
                obj.to_object + ["_m_", "_t_", "_p_"] + obj.to_type] ;
        end
    end
end     