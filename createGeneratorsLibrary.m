%%
% Feasibility study for size optimization of a geothermal/PV/wind/diesel hybrid power plant using the harmony search algorithm
% Majid Reza Naseh
%----------------------- BIOMASS---------------------
 % Optimal sizing of hybrid PV–diesel–biomass gasification plants for electrification of off-grid communities:
        % An efficient approach based on Benders’ decomposition
        % Feasibility study of an islanded microgrid in rural area consisting of PV, wind, biomass, and battery energy storage system
  % https://iea.blob.core.windows.net/assets/b5b73936-ee21-4e38-843b-8ba7430fbe92/TheFutureofGeothermal.pdf
% Geothermal carbon emission: Greenhouse Gas Emissions from Geothermal Power Production
  %%

function generators = createGeneratorsLibrary()
    generators = struct(...
        'biomassGen', struct(...
            'type', 'Gasification Plant', ...
            'rated_power', 50, ... % kW
            'number_of_units', 1, ...
            'biomass_consumption', 1.3, ... % kg/kWh
            'biochar_production', 0.15, ... % kg/kg of biomass
            'min_power', 0.3 * 50, ... % kW
            'min_startup_time', 3, ... % hr
            'min_shutdown_time', 72, ... % hr
            'ramping_rate', 0.3 * 50, ... % kW/hr
            'life', 15, ... % Years
            'efficiency', 0.21, ... % Efficiency
            'IC', 3000, ... % $/kW
            'OMC', 0.15 * 3000, ... % $/kWh/year
            'RC', 2000, ... % $/kW
            'SUC', 2000 ... % Startup Cost $
        ), ...
        'geothermal_flash_steam', struct(...
            'type', 'Flash steam power plant', ...
            'rated_power', 1000, ... % kW
            'number_of_units', 1, ...
            'input_temperature', 300, ... % °C
            'output_temperature', 48, ... % °C
            'life', 20, ... % Years
            'efficiency', 0.16, ... % Efficiency
            'IC', 500000, ... % $
            'OMC', 0.0005, ... % $/kWh
            'RC', 500000 * 0.05 ... % $
        ), ...
        'geothermal_binary_cycle', struct(...
            'type', 'Binary cycle power plant', ...
            'rated_power', 1000, ... % kW
            'number_of_units', 1, ...
            'input_temperature', 180, ... % °C
            'output_temperature', 48, ... % °C
            'life', 20, ... % Years
            'efficiency', 0.09, ... % Efficiency
            'IC', 500000, ... % $
            'OMC', 0.0005, ... % $/kWh
            'RC', 500000 * 0.05 ... % $
        ), ...
        'geothermal_dry_steam', struct(...
            'type', 'Dry steam power plant', ...
            'rated_power',  1000, ... % kW
            'number_of_units', 1, ...
            'input_temperature', 300, ... % °C
            'output_temperature', 48, ... % °C
            'life', 20, ... % Years
            'efficiency', 0.18, ... % Efficiency
            'IC', 500000, ... % $
            'OMC', 0.0005, ... % $/kWh
            'RC', 500000 * 0.05 ... % $
        ), ...
        'PHS', struct(...
            'type', 'pumped hydroelectric system', ...
            'rated_power',  1000, ... % kW
            'pumping_load', 1, ...% meter
            'input_temperature', 300, ... % °C
            'output_temperature', 48, ... % °C
            'life', 20, ... % Years
            'turbine_efficiency', 0.85, ... % Efficiency
            'pump_efficiency', 0.85, ... % Efficiency            
            'IC', 500000, ... % $
            'OMC', 0.0005, ... % $/kWh
            'RC', 500000 * 0.05 ... % $
        ), ...
            'hydrothermal_geothermal', struct(...
            'type', 'standard power plant', ...
            'N_wells', 10, ... % number of wells
            'embedded_CO2', 10+122, ...%gCO2/kWh (Plant cycle + Fuel Cyclce emissions)
            'well_depth', 3, ... % kilometer
            'Thermal_gradient_min', 10, ... % K/km 
            'Thermal_gradient_max', 30, ... % K/km 
            'in_flowrate', 80, ... % kg/s for 10 wells
            'in_flowrate_unit', 8, ... % kg/s for 1 wells
            'thermal_recovery_factor', 0.25, ... % https://doi.org/10.1007/s11053-022-10138-4
            'capacity_factor', 80, ... % electricity generation capacity factor   
            'input_temperature', 180, ... % °C
            'output_temperature', 48, ... % °C
            'life', 20, ... % Years  (for electricity generation), 25 for heat gen.
            'efficiency_BC', 0.12, ... % Binary Cycle Efficiency
            'efficiency_DS', 0.18, ... % Dry Steam Cycle Efficiency
            'efficiency_FC', 0.16, ... % Flash Steam Cycle Efficiency
            'Drilling_cost', 2000, ... % $/m
            'Stimualtion_cost', 2800, ... % $/m
            'Power_generation_cost', 2250, ... % $/kW
            'OMC', 2250*0.02, ... % $/kWh/year
            'RC', 2250*0.8 ... % $
        )...
    );
end
