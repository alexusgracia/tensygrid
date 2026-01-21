classdef FluidDirectionBool < handle & matlab.mixin.SetGet
    properties
        inletObject;
        booleanName;
        fluidDirectionModel;
    end

    methods
        function obj = FluidDirectionBool(props)
            arguments
                props.?FluidDirectionBool
            end
                set(obj,props)
                obj.fluidDirectionModel = obj.TotalfluidDirectionModel();
        end
        function fluidDirectionModel=TotalfluidDirectionModel(obj)

            H1 = hyCPN1();
            H1.F.input = [1,1];
            H1.F.boolean = [1,0];
            H1.phi.inequality = [2, -1];
            fluidDirectionModel = dmss(H1,0);
            fluidDirectionModel.inputName = obj.inletObject+"_m_inlet";
            fluidDirectionModel.booleanName = obj.booleanName;
        end
    end
end