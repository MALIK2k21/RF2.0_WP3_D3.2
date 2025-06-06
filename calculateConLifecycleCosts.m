function costs = calculateConLifecycleCosts(rated_power,Pr_life, int, inf)
            % Calculate lifecycle costs for wind turbines
            % Args:
            %   N_wt: Number of wind turbines
            %   Pr_life: Project lifetime in years
            %   int: Interest rate
            %   inf: Inflation rate
            life_unit=20;
            % rated_power=25; %kw

            IC=300;          % IN capital Cost($/kW)
            RC=0.8*IC;          % IN replacement Cost($/kW)
            OMC=0.02*IC;       % Operation & Maintenance Cost($/kW/year)---2% of Initial Cost

            N_TR = Pr_life / life_unit;

            % Capital Costs
            costs.capital = rated_power * IC;

            % Replacement Costs
            if life_unit == Pr_life
                costs.replacement = 0;
            else
                j = [1:1:(round(N_TR) + 1 - 1)];
                PWF_WT = (((1 + inf) / (1 + int)) .^ ((Pr_life * j) / ((round(N_TR) - 1) + 1)));
                costs.replacement = ( rated_power * RC) * PWF_WT;
            end

            % Operation and Maintenance Costs
            costs.maintenance = ( rated_power * OMC) * (((1 + inf) / (int - inf)) * (1 - (((1 + inf) / (1 + int)) ^ Pr_life)));

        end