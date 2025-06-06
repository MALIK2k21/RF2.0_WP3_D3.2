function costs = calculateESSLifecycleCosts(rated_cap_BAT,N_BT,ESS_life,IC,Pr_life, int, inf)
            % Calculate lifecycle costs for geothermal plant
            % Args:
            
            %   Pr_life: Project lifetime in years
            %   int: Interest rate
            %   inf: Inflation rate

            
        


            % Replacement cost is equal to capital costs
            OMC=0.02*IC;       % Operation & Maintenance Cost($/W/year)---2% of Initial Cost
            RC=0.8*IC;  
            N_TR = Pr_life / ESS_life;

            % Capital Costs
            costs.capital = N_BT*(rated_cap_BAT*IC) ;

            % Replacement Costs
            if ESS_life == Pr_life
                costs.replacement = 0;
            else
                j = [1:1:(round(N_TR) + 1 - 1)];
                PWF_WT = (((1 + inf) / (1 + int)) .^ ((Pr_life * j) / ((round(N_TR) - 1) + 1)));
                costs.replacement = (N_BT*(rated_cap_BAT* RC)) * PWF_WT;
            end

            % Operation and Maintenance Costs
            costs.maintenance = (N_BT*(rated_cap_BAT* OMC)) * (((1 + inf) / (int - inf)) * (1 - (((1 + inf) / (1 + int)) ^ Pr_life)));

        end