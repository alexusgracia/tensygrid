classdef MassFlowAbsolute < dynamicprops

    properties
        booleanName = "mass_flow_bool";
        pipeName = "1";
        fluidAbsModel;
    end
    
    methods
        function obj = MassFlowAbsolute(opts)
            % Enable the possibility of using name-value arguments to set
            % the class property
            arguments
                opts.?MassFlowAbsolute
            end
            % Now replace all the values of the property depending on the
            % given values. If not given, default values are used
            for prop = string(fieldnames(opts))'
                obj.(prop) = opts.(prop);
            end

            obj.fluidAbsM()
        end


        function fluidAbsM(obj)
            H1 = hyCPN1();
            H1.F.input = [1,1];
            H1.F.boolean = [1,0];
            H1.phi.inequality = [2, -1];
            fluidAbsM = dmss(H1,0);
            fluidAbsM.inputName = [obj.pipeName+"_m_inlet"];
            fluidAbsM.booleanName = [obj.booleanName];
            obj.fluidAbsModel = fluidAbsM;
        end
    end

end