classdef WindTurbinePowerCalculator
    properties
        % Wind Turbine Specifications loaded from library
        turbine_specs
    end

    methods
        function obj = WindTurbinePowerCalculator(turbine_name)
            % Constructor to initialize with specific turbine model
            % Args:
            %   turbine_name: Name of turbine model from library

            WT_LIBRARY = getWindTurbineLibrary();
            if ~isfield(WT_LIBRARY, turbine_name)
                error('Selected wind turbine not found in library');
            end
            obj.turbine_specs = WT_LIBRARY.(turbine_name);
        end

        function power = calculatePowerOutput(obj, data, hasOperationalData)
            % Calculate power output based on available data
            % Args:
            %   data: Structure containing either operational data or wind speed data
            %   hasOperationalData: Boolean indicating if operational data exists

            try
                if hasOperationalData
                    power = obj.processOperationalData(data);
                else
                    power = obj.calculateFromEnvironmental(data);
                end
            catch e
                error('Error in wind power calculation: %s', e.message);
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

        function power = calculateFromEnvironmental(obj, data)
            % Calculate power from wind speed data
            % Args:
            %   data: Structure with field 'windspeed'

            validateattributes(data, {'numeric'}, {'vector'});
            WS = data;

            % Initialize power output array
            power = zeros(size(WS));

            % Get turbine parameters from specs
            WTS_CI = obj.turbine_specs.cut_in_speed;
            WTS_CO = obj.turbine_specs.cut_out_speed;
            WTS_RATED = obj.turbine_specs.rated_speed;
            WT_RATED = obj.turbine_specs.rated_power;

            % Calculate power based on wind speed regions
            for j = 1:length(WS)
                if WS(j) < WTS_CI || WS(j) >= WTS_CO
                    % Below cut-in or above cut-out speed
                    power(j) = 0;
                elseif WS(j) > WTS_RATED && WS(j) <= WTS_CO
                    % At rated power
                    power(j) = WT_RATED;
                elseif WS(j) >= WTS_CI && WS(j) < WTS_RATED
                    % Power curve between cut-in and rated speed
                    power(j) = (WT_RATED) * ...
                        ((WS(j)^3 - WTS_CI^3) / ...
                        (WTS_RATED^3 - WTS_CI^3));
                end
            end
        end

        function [totalPower, dumpedPower] = calculateSystemOutput(obj, power, N_wt)
            % Calculate total system output and dumped power
            % Args:
            %   power: Power output per turbine
            %   N_wt: Number of wind turbines

            ratedSystemPower = (obj.turbine_specs.rated_power) * N_wt;
            totalPower = power * N_wt;
            dumpedPower = max(0, totalPower - ratedSystemPower);
            totalPower = min(totalPower, ratedSystemPower);
        end

        function costs = calculateLifecycleCosts(obj, N_wt, Pr_life, int, inf)
            % Calculate lifecycle costs for wind turbines
            % Args:
            %   N_wt: Number of wind turbines
            %   Pr_life: Project lifetime in years
            %   int: Interest rate
            %   inf: Inflation rate

            N_WTR = Pr_life / obj.turbine_specs.life;

            % Capital Costs
            costs.capital = N_wt * obj.turbine_specs.rated_power * obj.turbine_specs.IC;

            % Replacement Costs
            if obj.turbine_specs.life == Pr_life
                costs.replacement = 0;
            else
                j = [1:1:(round(N_WTR) + 1 - 1)];
                PWF_WT = (((1 + inf) / (1 + int)) .^ ((Pr_life * j) / ((round(N_WTR) - 1) + 1)));
                costs.replacement = (N_wt * obj.turbine_specs.rated_power * obj.turbine_specs.RC) * PWF_WT;
            end

            % Operation and Maintenance Costs
            costs.maintenance = (N_wt * obj.turbine_specs.rated_power * obj.turbine_specs.OMC) * (((1 + inf) / (int - inf)) * (1 - (((1 + inf) / (1 + int)) ^ Pr_life)));

        end

    end
end