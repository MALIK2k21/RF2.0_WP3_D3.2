classdef BiomassGeneratorPowerCalculator
    properties
        % Biomass Generator Specifications
        generator_specs
    end

    methods
        function obj = BiomassGeneratorPowerCalculator(generator_specs)
            % Constructor to initialize with specific generator model
            % Args:
            %   generator_specs: Structure containing generator specifications

            obj.generator_specs = generator_specs;
        end

        function power = calculatePowerOutput(obj, data, hasOperationalData)
            % Calculate power output based on available data
            % Args:
            %   data: Structure containing either operational data or bio waste data
            %   hasOperationalData: Boolean indicating if operational data exists

            try
                if hasOperationalData
                    power = obj.processOperationalData(data);
                else
                    power = obj.calculateFromBioWaste(data);
                end
            catch e
                error('Error in power calculation: %s', e.message);
            end
        end

        function power = processOperationalData(data)
            % Process existing operational data
            % Args:
            %   data: Table/array containing power generation data

            if ~isnumeric(data)
                data = table2array(data);
            end
            power = data;
        end

        function power = calculateFromBioWaste(obj, data)
            % Calculate power from bio waste data
            % Args:
            %   data: Structure with field 'bioWaste'

            validateattributes(data, {'numeric'}, {'scalar'});
            bioWaste = data;

            % Get generator parameters from specs
            calorific_value = obj.generator_specs.calorific_value;
            total_conversion_efficiency = obj.generator_specs.total_conversion_efficiency;
            operation_hours_per_day = obj.generator_specs.operation_hours_per_day;

            % Calculate power output
            power = (bioWaste * calorific_value * total_conversion_efficiency * operation_hours_per_day) / 876;
        end

        function costs = calculateLifecycleCosts(obj, N_gen, Pr_life, int, inf)
            % Calculate lifecycle costs for biomass generators
            % Args:
            %   N_gen: Number of biomass generators
            %   Pr_life: Project lifetime in years
            %   int: Interest rate
            %   inf: Inflation rate

            % Calculate lifecycle costs using similar logic as in WindTurbinePowerCalculator


            N_GEN = Pr_life / obj.generator_specs.life;

            % Capital Costs
            costs.capital = N_gen * obj.generator_specs.rated_power * obj.generator_specs.capital_cost;

            % Replacement Costs
            if obj.generator_specs.life == Pr_life
                costs.replacement = 0;
            else
                j = [1:1:(round(N_GEN) + 1 - 1)];
                PWF_GEN = (((1 + inf) / (1 + int)) .^ ((Pr_life * j) / ((round(N_GEN) - 1) + 1)));
                costs.replacement = (N_gen * obj.generator_specs.rated_power * obj.generator_specs.replacement_cost) * PWF_GEN;
            end

            % Operation and Maintenance Costs
            costs.maintenance = (N_gen * obj.generator_specs.rated_power * obj.generator_specs.om_cost) * (((1 + inf) / (int - inf)) * (1 - (((1 + inf) / (1 + int)) ^ Pr_life)));

        end

    end
end
