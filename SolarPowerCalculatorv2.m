classdef SolarPowerCalculatorv2
    properties
        % PV Module Specifications loaded from library
        module_specs
    end

    methods
        function obj = SolarPowerCalculatorv2(module_name)
            % Constructor to initialize with specific module model
            % Args:
            %   module_name: Name of module model from library

            PV_LIBRARY = getPVModuleLibrary();
            if ~isfield(PV_LIBRARY, module_name)
                error('Selected PV module not found in library');
            end
            obj.module_specs = PV_LIBRARY.(module_name);
        end

        function power = calculatePowerOutput(obj, data, hasOperationalData, N_PV)
            % Calculate power output based on available data
            % Args:
            %   data: Structure containing either operational data or environmental data
            %   hasOperationalData: Boolean indicating if operational data exists
            %   N_PV: Surface area of PV modules

            try
                if hasOperationalData
                    power = SolarPowerCalculatorv2.processOperationalData(data);
                else
                    power = obj.calculateFromEnvironmental(data, N_PV);
                end
            catch e
                error('Error in power calculation: %s', e.message);
            end
        end
    end

    methods (Static)
        function power = processOperationalData(data)
            % Process data (operational or environmental)
            % Args:
            %   data: Table/array containing power generation data or environmental data

            if ~isnumeric(data)
                data = table2array(data);
            end
            power = data;
        end
    end

    methods
        function power = calculateFromEnvironmental(obj, data, N_PV)
            % Calculate power from environmental data
            % Args:
            %   data: Structure with fields 'temperature' and 'irradiance'
            %   N_PV: Surface area of PV modules

            validateattributes(data(:,1), {'numeric'}, {'vector'});
            validateattributes(data(:,2), {'numeric'}, {'vector'});

            % Calculate cell temperature
            Tc = data(:,1) + (data(:,2) .* ((obj.module_specs.ncot - 20) / 800));

            % Calculate temperature losses
            temp_loss = abs(1 + (obj.module_specs.temp_coeff .* (Tc - obj.module_specs.ref_temp)));

            % Calculate Total PV Plant power output
            power = data(:,2) .* N_PV.*obj.module_specs.ratedpowerkw.* obj.module_specs.SurfaceAream2.*temp_loss .* ...
                obj.module_specs.panel_efficiency .* obj.module_specs.Degrading_Factor;
        end

        function costs = calculateLifecycleCosts(obj, N_PV, Pr_life, int, inf)
            % Calculate lifecycle costs
            % Args:
            %   N_PV: Surface area of PV modules
            %   Pr_life: Project lifetime in years
            %   int: Interest rate
            %   inf: Inflation rate

            N_PVTR = Pr_life / obj.module_specs.life;

            costs.capital = N_PV* obj.module_specs.IC * obj.module_specs.ratedpowerkw;

            % Replacement Costs
            if obj.module_specs.life == Pr_life
                costs.replacement = 0;
            else
                j = [1:1:(round(N_PVTR) + 1 - 1)];
                PWF = (((1 + int) / (1 + inf)) .^ ((Pr_life * j) / ((round(N_PVTR) - 1) + 1)));
                costs.replacement = (N_PV* obj.module_specs.ratedpowerkw * obj.module_specs.RC) * PWF;
            end

            % Operation and Maintenance Costs
            costs.maintenance = (N_PV* obj.module_specs.ratedpowerkw * obj.module_specs.OMC) * (((1 + inf) / (int - inf)) * (1 - (((1 + inf) / (1 + int)) ^ Pr_life)));
        end
    end
end
