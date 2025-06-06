function costs = calculateGTLifecycleCostsv1(N_well,well_depth,GTcap,Pr_life, int, inf)
            % Calculate lifecycle costs for geothermal plant
            % Args:
            % GTcap: geothrmal plant capc
            %   Pr_life: Project lifetime in years
            %   int: Interest rate
            %   inf: Inflation rate

            life_unit=30;
        
            % IC_Well=2000;          % capital Cost for drilling and "NO" stimulating cost due to hydrothermal system
            % IC_Well=2000+2800;          % capital Cost for drilling and stimulating well($/m)
            % IC_Gen=2250;          % IN capital Cost($/kW)
            IC_Gen=8000;          % IN capital Cost($/kW)
            % Replacement cost is equal to capital costs
            OMC=0.02*IC_Gen;       % Operation & Maintenance Cost($/kW/year)---2% of Initial Cost

            N_TR = Pr_life / life_unit;
            % Here N_well essentially emans a new 680kW plant which comes
            % with 10 more tehrmalwells

            % Capital Costs
            % costs.capital = (N_well *10* well_depth * (IC_Well)) + (GTcap*IC_Gen);
            costs.capital = (N_well*(GTcap*IC_Gen));
            % Replacement Costs
            if life_unit == Pr_life
                costs.replacement = 0;
            else
                j = [1:1:(round(N_TR) + 1 - 1)];
                PWF_WT = (((1 + inf) / (1 + int)) .^ ((Pr_life * j) / ((round(N_TR) - 1) + 1)));
                costs.replacement = ((N_well*(GTcap*IC_Gen*0.8))) * PWF_WT;
            end

            % Operation and Maintenance Costs
            costs.maintenance = (N_well*GTcap* OMC) * (((1 + inf) / (int - inf)) * (1 - (((1 + inf) / (1 + int)) ^ Pr_life)));

        end