function f = objective_function(x, PVSF, PLC, ESSAg, BSD, SSD, NLPS, battery_type, facility, ESS_conf, L1, power_profiles,factor, varargin)



% Parse and validate inputs
config = parseAndValidateInputs(x, PVSF, PLC, ESSAg, BSD, SSD, NLPS, battery_type, facility, ESS_conf,L1, power_profiles,factor);

% Initialize system parameters and specifications
[systemParams, batterySpecs, supercapSpecs] = initializeSystemParameters(config);

% Load external data and generate power profiles
powerProfiles = generatePowerProfiles(systemParams, config);


%%
% Design power allocation filter for hybrid ESS
filterParams = designPowerAllocationFilter(powerProfiles, config);

% Run main simulation loop
simulationResults = runMainSimulation(systemParams, batterySpecs, supercapSpecs, powerProfiles, filterParams, config);

% Calculate lifecycle costs and emissions
[Total_costs, costs,grid_costs,penalties,PSC,Total_emissions] = calculateFinalCostsAndEmissions(simulationResults, systemParams, config,powerProfiles);
Total_emissions=Total_emissions/1000; %gCO-eq/kWh-->kgCO2-eq/kWh
Total_costs1=Total_costs+PSC;
PVR=((max(simulationResults.PLF_seq)-min(simulationResults.PLF_seq))/(max(simulationResults.PLF_seq)))*100;  %Peak valley ratio of net load
f=[Total_costs1,Total_emissions];
end


%% Configuration and Initialization Functions

function config = parseAndValidateInputs(x,PVSF, PLC, ESSAg, BSD, SSD, NLPS, battery_type, facility, ESS_conf, L1, power_profiles, factor,varargin)
% Existing parameters
if strcmp(ESS_conf, 'BAT')
    config.N_BT = ceil(x(1));
    config.freq = x(2);
    config.N_SC=0;
    
    % if strcmp(power_profiles.wind.type, 'production')
    %     config.N_WT = calculateWTUnits(config.WPSF);
    % elseif strcmp(power_profiles.wind.type, 'meteorological')
    %     config.N_WT = x(3);
    % else
    %     error('Unknown wind type: %s. Use "production" or "meteorological"', power_profiles.wind.type);
    % end
    % -------------------------------------------------------
elseif strcmp(ESS_conf, 'SC')
    config.N_SC = ceil(x(1));
    config.freq = x(2);
    config.N_BT=0;
    
    % if strcmp(power_profiles.wind.type, 'production')
    %     config.N_WT = calculateWTUnits(config.WPSF);
    % elseif strcmp(power_profiles.wind.type, 'meteorological')
    %     config.N_WT = x(3);
    % else
    %     error('Unknown wind type: %s. Use "production" or "meteorological"', power_profiles.wind.type);
    % end
    % -------------------------------------------------------
elseif strcmp(ESS_conf, 'HYB')

    config.N_BT = ceil(x(1));
    config.N_SC = ceil(x(2));
    config.freq = x(3);
    % WPSF=0;
    % if strcmp(power_profiles.wind.type, 'production')
    %     config.N_WT = calculateWTUnits(config.WPSF);
    % elseif strcmp(power_profiles.wind.type, 'meteorological')
    %     config.N_WT = x(4);
    % else
    %     error('Unknown wind type: %s. Use "production" or "meteorological"', power_profiles.wind.type);
    % end

end

WPSF=0;
% System configuration
% config.PVSF = x(4);
config.PVSF = PVSF;
config.WPSF = WPSF;
config.do_PLC = PLC;
config.do_ESSAg = ESSAg;
config.C_BSD = BSD;
config.C_SSD = SSD;
config.do_NLPS = NLPS;
config.battery_type = battery_type;
config.facility = facility;
config.ESS_conf = ESS_conf;
config.factor = factor;
% PV
% Solar
if strcmp(power_profiles.solar.type, 'production')
    config.N_PV = calculatePVUnits(config.PVSF);
elseif strcmp(power_profiles.solar.type, 'meteorological')
    config.N_PV = x(5);
else
    error('Unknown solar type: %s. Use "production" or "meteorological"', power_profiles.solar.type);
end

% Wind
if strcmp(power_profiles.wind.type, 'production')
    config.N_WT = calculateWTUnits(config.WPSF);
elseif strcmp(power_profiles.wind.type, 'meteorological')
    config.N_WT = x(4);
else
    error('Unknown wind type: %s. Use "production" or "meteorological"', power_profiles.wind.type);
end


%Geothermal
config.N_wells=0;


% Load and power profiles
config.L1 = L1;
config.power_profiles = power_profiles;

% Validate power profiles structure
if ~isfield(power_profiles, 'type')
    error('power_profiles must have a "type" field: "production" or "meteorological"');
end

% Handle varargin for optional parameters
if ~isempty(varargin)
    for idx = 1:2:length(varargin)
        config.(varargin{idx}) = varargin{idx+1};
    end
end
end


function N_PV = calculatePVUnits(PVSF)
% Calculate number of PV units based on scaling factor
PV_PLANT_CAPACITY_KW = 250;
pvLib = getPVModuleLibrary();
module_rating = pvLib.Mono_SI_SOLARIA.ratedpowerkw;
N_PV = PV_PLANT_CAPACITY_KW * PVSF / module_rating;
end

function N_WT = calculateWTUnits(WPSF)
% Calculate number of PV units based on scaling factor
WP_PLANT_CAPACITY_KW = 250;
WTLib = getWindTurbineLibrary();
module_rating = WTLib.wt27.rated_power;
N_WT = WP_PLANT_CAPACITY_KW * WPSF / module_rating;
end



function [systemParams, batterySpecs, supercapSpecs] = initializeSystemParameters(config)
% Initialize all system parameters and component specifications

% Time and financial parameters
systemParams = struct();
systemParams.HOURS_PER_YEAR = 8760;
systemParams.PROJECT_LIFE_YEARS = 30;
systemParams.INTEREST_RATE = 0.034;        % 3.4%
systemParams.INFLATION_RATE = 0.024;       % 2.4%
systemParams.TIME_STEP_HOURS = 1;
systemParams.INVERTER_EFFICIENCY = 0.95+(0.95*config.factor);
systemParams.BDCF = 168;                   % Battery degradation calculation frequency (weekly)

% Battery system parameters
[batterySpecs, supercapSpecs] = initializeESSSpecifications(config, systemParams);

% Store configuration in system parameters
systemParams.config = config;
systemParams.batteff=batterySpecs.efficiency;
systemParams.supercapeff=supercapSpecs.efficiency;
end

function [batterySpecs, supercapSpecs] = initializeESSSpecifications(config, systemParams)
% Initialize battery and supercapacitor specifications
(getESSLibrary().(config.battery_type).SOC_max-getESSLibrary().(config.battery_type).SOC_min)/100;
ESS_lib = getESSLibrary();

% Battery specifications
batterySpecs = struct();
batterySpecs.rated_cap_unit = 100;          % kWh per unit
batterySpecs.C_rate = 1;                    % 1C discharge rate
batterySpecs.rated_power_unit = batterySpecs.rated_cap_unit * batterySpecs.C_rate;

batData = ESS_lib.(config.battery_type);
batterySpecs.SOC_max = config.N_BT * batterySpecs.rated_cap_unit * (batData.SOC_max/100);
batterySpecs.SOC_min = config.N_BT * batterySpecs.rated_cap_unit * (batData.SOC_min/100);
batterySpecs.efficiency = batData.round_trip_efficiency/100+(config.factor*(batData.round_trip_efficiency/100));
batterySpecs.max_charge = config.N_BT * batterySpecs.rated_power_unit;
batterySpecs.max_discharge = config.N_BT * batterySpecs.rated_power_unit;
batterySpecs.initial_capacity = config.N_BT * batterySpecs.rated_cap_unit;
batterySpecs.self_discharge_rate = (batData.self_discharge_daily / 24) * config.C_BSD;
batterySpecs.aging_coeffs = getBatteryAgingCoefficients(config.battery_type);
batterySpecs.specs = batData;

% Supercapacitor specifications
supercapSpecs = struct();
supercapSpecs.rated_cap_unit = 50;          % kWh per unit
supercapSpecs.C_rate = 1;
supercapSpecs.rated_power_unit = supercapSpecs.rated_cap_unit * supercapSpecs.C_rate;

scData = ESS_lib.SC;
supercapSpecs.SOC_max = config.N_SC * supercapSpecs.rated_cap_unit * (scData.SOC_max/100);
supercapSpecs.SOC_min = config.N_SC * supercapSpecs.rated_cap_unit * (scData.SOC_min/100);
supercapSpecs.efficiency = scData.round_trip_efficiency/100+(config.factor*(scData.round_trip_efficiency/100));
supercapSpecs.max_charge = config.N_SC * supercapSpecs.rated_power_unit;
supercapSpecs.max_discharge = config.N_SC * supercapSpecs.rated_power_unit;
supercapSpecs.initial_capacity = config.N_SC * supercapSpecs.rated_cap_unit;
supercapSpecs.self_discharge_rate = scData.self_discharge * config.C_SSD;
supercapSpecs.specs = scData;
end

%% Power Profile Generation
function powerProfiles = generatePowerProfiles(systemParams, config)
% Generate or load power generation and demand profiles
% Supports hybrid input: mix of production profiles and meteorological data

powerProfiles = struct();

% Load demand profile
powerProfiles.load = config.L1 / systemParams.INVERTER_EFFICIENCY;

% Check if using new hybrid structure or legacy structure
if isfield(config.power_profiles, 'solar') && isstruct(config.power_profiles.solar)
    % New hybrid structure: each source has its own type
    powerProfiles = handleHybridPowerProfiles(powerProfiles, config, systemParams);
else
    % Legacy structure: global type for all sources
    switch lower(config.power_profiles.type)
        case 'production'
            powerProfiles = handleDirectPowerProfiles(powerProfiles, config, systemParams);
        case 'meteorological'
            powerProfiles = handleMeteorologicalData(powerProfiles, config, systemParams);
        otherwise
            error('Unknown power_profiles.type: %s. Use "production" or "meteorological"', config.power_profiles.type);
    end
end

% Calculate total renewable generation and net load
powerProfiles.total_renewable = powerProfiles.solar + powerProfiles.wind + powerProfiles.geothermal;
powerProfiles.net_load = powerProfiles.load - powerProfiles.total_renewable;
end

function powerProfiles = handleHybridPowerProfiles(powerProfiles, config, systemParams)
% Handle hybrid power profiles where each source can be production or meteorological

profiles = config.power_profiles;

% --- Solar Power Processing ---
if isfield(profiles, 'solar') && isstruct(profiles.solar)
    switch lower(profiles.solar.type)
        case 'production'
            % Direct production profile
            if isfield(profiles.solar, 'data') && ~isempty(profiles.solar.data)
                powerProfiles.solar = profiles.solar.data * config.PVSF * systemParams.INVERTER_EFFICIENCY;
                validateProfileLength(powerProfiles.solar, systemParams.HOURS_PER_YEAR, 'Solar production');
            else
                powerProfiles.solar = zeros(systemParams.HOURS_PER_YEAR, 1);
            end

        case 'meteorological'
            % Generate from meteorological data
            if isfield(profiles.solar, 'data') && size(profiles.solar.data, 2) >= 2
                T = profiles.solar.data(:,1);
                IR = profiles.solar.data(:,2);
                powerProfiles.solar = generateSolarFromMeteo(T, IR, config, systemParams);
            elseif isfield(profiles.solar, 'T') && isfield(profiles.solar, 'IR')
                powerProfiles.solar = generateSolarFromMeteo(profiles.solar.T, profiles.solar.IR, config, systemParams);
            else
                error('Solar meteorological data requires T and IR fields or data matrix with 2 columns');
            end

        otherwise
            error('Unknown solar type: %s. Use "production" or "meteorological"', profiles.solar.type);
    end
else
    powerProfiles.solar = zeros(systemParams.HOURS_PER_YEAR, 1);
end

% --- Wind Power Processing ---
if isfield(profiles, 'wind') && isstruct(profiles.wind)
    switch lower(profiles.wind.type)
        case 'production'
            % Direct production profile
            if isfield(profiles.wind, 'data') && ~isempty(profiles.wind.data)
                powerProfiles.wind = profiles.wind.data * systemParams.INVERTER_EFFICIENCY;
                validateProfileLength(powerProfiles.wind, systemParams.HOURS_PER_YEAR, 'Wind production');
            else
                powerProfiles.wind = zeros(systemParams.HOURS_PER_YEAR, 1);
            end

        case 'meteorological'
            % Generate from wind speed data
            if isfield(profiles.wind, 'data') && ~isempty(profiles.wind.data)
                powerProfiles.wind = generateWindFromMeteo(profiles.wind.data, config,systemParams);
            elseif isfield(profiles.wind, 'WS_ref') && ~isempty(profiles.wind.WS_ref)
                powerProfiles.wind = generateWindFromMeteo(profiles.wind.WS_ref, config,systemParams);
            else
                error('Wind meteorological data requires WS_ref field or data vector');
            end

        otherwise
            error('Unknown wind type: %s. Use "production" or "meteorological"', profiles.wind.type);
    end
else
    powerProfiles.wind = zeros(systemParams.HOURS_PER_YEAR, 1);
end

% --- Geothermal Power Processing ---
if isfield(profiles, 'geothermal') && isstruct(profiles.geothermal)
    switch lower(profiles.geothermal.type)
        case 'production'
            % Direct production profile
            if isfield(profiles.geothermal, 'data') && ~isempty(profiles.geothermal.data)
                powerProfiles.geothermal = profiles.geothermal.data * systemParams.INVERTER_EFFICIENCY;
                validateProfileLength(powerProfiles.geothermal, systemParams.HOURS_PER_YEAR, 'Geothermal production');
            else
                powerProfiles.geothermal = zeros(systemParams.HOURS_PER_YEAR, 1);
            end

        case 'meteorological'
            % Generate from meteorological data (if applicable)
            powerProfiles.geothermal = generateGeothermalProfile(config, systemParams.HOURS_PER_YEAR);

        otherwise
            error('Unknown geothermal type: %s. Use "production" or "meteorological"', profiles.geothermal.type);
    end
else
    powerProfiles.geothermal = generateGeothermalProfile(config, systemParams.HOURS_PER_YEAR);
end

% Apply bounds checking
powerProfiles.solar = max(powerProfiles.solar, 0);
powerProfiles.wind = max(powerProfiles.wind, 0);
powerProfiles.geothermal = max(powerProfiles.geothermal, 0);
end

function validateProfileLength(profile, expected_length, profile_name)
% Validate that profile has correct length
if length(profile) ~= expected_length
    error('%s profile length (%d) must equal HOURS_PER_YEAR (%d)', ...
        profile_name, length(profile), expected_length);
end
end

% Keep your existing helper functions unchanged:
function powerProfiles = handleDirectPowerProfiles(powerProfiles, config, systemParams)
% Handle direct power production profiles (legacy support)

profiles = config.power_profiles;

% Solar power profile
if isfield(profiles, 'solar') && ~isempty(profiles.solar)
    powerProfiles.solar = profiles.solar * config.PVSF * systemParams.INVERTER_EFFICIENCY;
    validateProfileLength(powerProfiles.solar, systemParams.HOURS_PER_YEAR, 'Solar production');
else
    powerProfiles.solar = zeros(systemParams.HOURS_PER_YEAR, 1);
end

% Wind power profile
if isfield(profiles, 'wind') && ~isempty(profiles.wind)
    powerProfiles.wind = profiles.wind * systemParams.INVERTER_EFFICIENCY;
    validateProfileLength(powerProfiles.wind, systemParams.HOURS_PER_YEAR, 'Wind production');
else
    powerProfiles.wind = zeros(systemParams.HOURS_PER_YEAR, 1);
end

% Geothermal power profile
if isfield(profiles, 'geothermal') && ~isempty(profiles.geothermal)
    powerProfiles.geothermal = profiles.geothermal * systemParams.INVERTER_EFFICIENCY;
    validateProfileLength(powerProfiles.geothermal, systemParams.HOURS_PER_YEAR, 'Geothermal production');
else
    powerProfiles.geothermal = zeros(systemParams.HOURS_PER_YEAR, 1);
end

% Apply bounds checking
powerProfiles.solar = max(powerProfiles.solar, 0);
powerProfiles.wind = max(powerProfiles.wind, 0);
powerProfiles.geothermal = max(powerProfiles.geothermal, 0);
end

function powerProfiles = handleMeteorologicalData(powerProfiles, config, systemParams)
% Handle meteorological data to generate power profiles (legacy support)

profiles = config.power_profiles;

% Solar power from temperature and irradiance
if isfield(profiles, 'T') && isfield(profiles, 'IR')
    powerProfiles.solar = generateSolarFromMeteo(profiles.T, profiles.IR, config, systemParams);
else
    powerProfiles.solar = zeros(systemParams.HOURS_PER_YEAR, 1);
end

% Wind power from wind speed
if isfield(profiles, 'WS_ref') && ~isempty(profiles.WS_ref)
    powerProfiles.wind = generateWindFromMeteo(profiles.WS_ref, config,systemParams);
else
    powerProfiles.wind = zeros(systemParams.HOURS_PER_YEAR, 1);
end

% Geothermal (typically constant)
powerProfiles.geothermal = generateGeothermalProfile(config, systemParams.HOURS_PER_YEAR);
end

% Keep your existing helper functions unchanged:
function solarProfile = generateSolarFromMeteo(T, IR, config, systemParams)
% Generate solar power from meteorological data (original method)

% Validate input lengths
if length(T) ~= systemParams.HOURS_PER_YEAR || length(IR) ~= systemParams.HOURS_PER_YEAR
    error('Temperature and irradiance profiles must be %d hours long', systemParams.HOURS_PER_YEAR);
end

spc = SolarPowerCalculatorv2('Mono_SI_SOLARIA');
N_PV = calculatePVUnits(config.PVSF);
rawPV = spc.calculatePowerOutput([T IR], false, N_PV);
solarProfile = systemParams.INVERTER_EFFICIENCY * rawPV;
end

function windProfile = generateWindFromMeteo(WS_ref, config, systemParams)
% Generate wind power from wind speed (original method)

if config.N_WT == 0
    windProfile = zeros(length(WS_ref), 1);
    return;
end

windProfile = computeWindPower(WS_ref, config.N_WT, systemParams);
end

function geothermalProfile = generateGeothermalProfile(config, hours_per_year)
% Generate geothermal power profile
geothermalProfile = computeGeothermalPower(config.N_wells, hours_per_year);
end
function [windGen] = computeWindPower(WS_ref,N_WT, systemParams)
libWT = getWindTurbineLibrary();
hub_h = libWT.wt27.hub_height;
ref_h = 16;  % original reference height
WS = WS_ref .* ((hub_h/ref_h)^(1/7));
wtc = WindTurbinePowerCalculator('wt27');
windSeries = wtc.calculatePowerOutput(WS, false);
windGen= wtc.calculateSystemOutput(windSeries, N_WT)*systemParams.INVERTER_EFFICIENCY;
% Note: original code used N_wt=0, so net wind=0 unless user changes.
end

%% Power Allocation Filter Design

function filterParams = designPowerAllocationFilter(powerProfiles, config)
% Design power allocation filter for hybrid ESS operation

filterParams = struct();

switch config.ESS_conf
    case 'HYB'
        filterParams = designHybridFilter(powerProfiles.net_load, config.freq);
    case 'BAT'
        filterParams.P_slow = powerProfiles.net_load;
        filterParams.P_fast = zeros(size(powerProfiles.net_load));
    case 'SC'
        filterParams.P_fast = powerProfiles.net_load;
        filterParams.P_slow = zeros(size(powerProfiles.net_load));
end
end

function filterParams = designHybridFilter(net_load, freq)
% Design low-pass filter for hybrid ESS power allocation

% Filter design parameters
Fs = 1;                              % Sampling frequency (samples per day)
Fn = Fs / 2;                         % Nyquist frequency
cutoff_norm = freq / Fn;             % Normalized cutoff frequency
% cutoff_norm =round((cutoff_norm ),4);
% if cutoff_norm > 1
%     cutoff_norm_incorrect=cutoff_norm;
%     cutoff_norm=1;
% end
filter_order = 2;
% cutoff_norm =round((cutoff_norm ),4);

% Design Butterworth filter
[b, a] = butter(filter_order, cutoff_norm, 'low');

% Apply zero-phase filtering
P_slow_raw = filtfilt(b, a, net_load);
filterParams.P_slow_raw=P_slow_raw ;
filterParams.P_fast_raw=net_load -P_slow_raw ;
alpha_raw = P_slow_raw ./ (net_load + eps);
alpha = min(max(alpha_raw, 0), 1);

% Calculate power allocation
filterParams.P_slow = alpha .* net_load;
filterParams.P_fast = net_load - filterParams.P_slow;
filterParams.alpha = alpha;
% filterParams.FcIn=cutoff_norm_incorrect;
filterParams.Fc=cutoff_norm;
end

%% Main Simulation Loop

function results = runMainSimulation(systemParams, batterySpecs, supercapSpecs, powerProfiles, filterParams, config)
% Run main hourly simulation loop

% Initialize simulation state
simState = initializeSimulationState(systemParams, batterySpecs, supercapSpecs);

% Pre-calculate weekly indices for efficiency
weeklyIndices = calculateWeeklyIndices(systemParams.HOURS_PER_YEAR, systemParams.BDCF);

% Main hourly loop
for hour = 2:systemParams.HOURS_PER_YEAR
    % Update self-discharge losses
    simState = updateSelfDischargeLosses(simState, hour, batterySpecs, supercapSpecs);

    % Determine charging or discharging mode and calculate power flows
    simState = calculatePowerFlows(simState, hour, filterParams, batterySpecs, supercapSpecs, systemParams);

    % Update SOC for both ESS types
    simState = updateSOCStates(simState, hour, batterySpecs, supercapSpecs);

    % Calculate system-level metrics
    simState = updateSystemMetrics(simState, hour, powerProfiles, systemParams);

    % Perform weekly aging calculations if enabled
    if config.do_ESSAg && (mod(hour, systemParams.BDCF) == 0 || hour == systemParams.HOURS_PER_YEAR)
        simState = performWeeklyAging(simState, hour, weeklyIndices, batterySpecs, systemParams);
    end
end

results = simState;
end

function simState = initializeSimulationState(systemParams, batterySpecs, supercapSpecs)
% Initialize all simulation state variables

len = systemParams.HOURS_PER_YEAR;

% SOC tracking
simState.SOC_bat = zeros(len, 1);
simState.SOC_sc = zeros(len, 1);
simState.SOC_bat(1) = batterySpecs.SOC_min;
simState.SOC_sc(1) = supercapSpecs.SOC_min;

% Power flow tracking
simState.P_ch_bat = zeros(len, 1);      % Battery charging power
simState.P_dh_bat = zeros(len, 1);      % Battery discharging power
simState.P_ch_sc = zeros(len, 1);       % Supercap charging power
simState.P_dh_sc = zeros(len, 1);       % Supercap discharging power

% System metrics
simState.LSBG = zeros(len, 1);          % Load served by grid
simState.dumped_power = zeros(len, 1);   % Dumped excess power
simState.P_net = zeros(len, 1);         % Net grid power
simState.P_net_noGrid = zeros(len, 1);  % Net power without grid
simState.PLF_seq = zeros(len, 1);       % Peak load following sequence

% Loss tracking
simState.self_discharge_loss_bat = zeros(len, 1);
simState.self_discharge_loss_sc = zeros(len, 1);

% State counters
simState.idle_cnt_bat = zeros(len, 1);
simState.idle_cnt_sc = zeros(len, 1);
simState.u_state = zeros(len, 1);       % 0=charging, 1=discharging

% Violation counters
simState.SOC_min_bat_violations = 0;
simState.SOC_min_sc_violations = 0;
simState.EOL_violations = 0;

% Aging tracking (if needed)
num_weeks = ceil(len / systemParams.BDCF);
simState.deg_weekly = zeros(num_weeks, 1);
simState.cap_hist = zeros(num_weeks, 1);
simState.cum_deg_hist = zeros(num_weeks, 1);
simState.degradation_cost = 0;
end

function weeklyIndices = calculateWeeklyIndices(hours_per_year, bdcf)
% Pre-calculate weekly indices for efficient aging calculations

num_weeks = ceil(hours_per_year / bdcf);
weeklyIndices = cell(num_weeks, 1);

for week = 1:num_weeks
    start_idx = (week - 1) * bdcf + 1;
    end_idx = min(week * bdcf, hours_per_year);
    weeklyIndices{week} = start_idx:end_idx;
end
end

%% Power Flow Calculation Functions

function simState = calculatePowerFlows(simState, hour, filterParams, batterySpecs, supercapSpecs, systemParams)
% Calculate power flows for current hour based on charging/discharging mode

% Determine power allocation for hybrid ESS
[bat_power_target, sc_power_target] = calculatePowerTargets(hour, filterParams, systemParams);

if bat_power_target <= 0 && sc_power_target <= 0
    % Charging mode (negative power means charging)
    simState = handleChargingMode(simState, hour, -bat_power_target, -sc_power_target, batterySpecs, supercapSpecs, systemParams);
else
    % Discharging mode (positive power means discharging)
    simState = handleDischargingMode(simState, hour, bat_power_target, sc_power_target, batterySpecs, supercapSpecs, systemParams);
end

% Update idle counters
simState = updateIdleCounters(simState, hour);

end

function [bat_target, sc_target] = calculatePowerTargets(hour, filterParams, systemParams)
% Calculate power targets for battery and supercapacitor

bat_target = filterParams.P_slow(hour);
sc_target = filterParams.P_fast(hour);
end

function simState = handleChargingMode(simState, hour, bat_charge_target, sc_charge_target, batterySpecs, supercapSpecs, systemParams)
% Handle charging mode for both ESS types

% Apply inverter efficiency to charging targets
bat_charge_available = bat_charge_target * systemParams.INVERTER_EFFICIENCY;
sc_charge_available = sc_charge_target * systemParams.INVERTER_EFFICIENCY;

% Calculate actual charging power considering constraints
simState.P_ch_bat(hour) = min([
    bat_charge_available * batterySpecs.efficiency,
    batterySpecs.SOC_max - simState.SOC_bat(hour-1),
    batterySpecs.max_charge * systemParams.TIME_STEP_HOURS
    ]);

simState.P_ch_sc(hour) = min([
    sc_charge_available * supercapSpecs.efficiency,
    supercapSpecs.SOC_max - simState.SOC_sc(hour-1),
    supercapSpecs.max_charge * systemParams.TIME_STEP_HOURS
    ]);

% Calculate dumped power
total_available = bat_charge_available + sc_charge_available;
total_used = (simState.P_ch_bat(hour) / batterySpecs.efficiency) + (simState.P_ch_sc(hour) / supercapSpecs.efficiency);
simState.dumped_power(hour) = max(0, total_available - total_used);

% Reset discharging variables
simState.P_dh_bat(hour) = 0;
simState.P_dh_sc(hour) = 0;
simState.LSBG(hour) = 0;
simState.u_state(hour) = 0;  % Charging state
end

function simState = handleDischargingMode(simState, hour, bat_discharge_target, sc_discharge_target, batterySpecs, supercapSpecs, systemParams)
% Handle discharging mode for both ESS types

% Check for SOC violations
if simState.SOC_bat(hour-1) < batterySpecs.SOC_min
    simState.SOC_min_bat_violations = simState.SOC_min_bat_violations + 1;
end
if simState.SOC_sc(hour-1) < supercapSpecs.SOC_min
    simState.SOC_min_sc_violations = simState.SOC_min_sc_violations + 1;
end

% Calculate available energy from both ESS
bat_available = max(0, simState.SOC_bat(hour-1) - batterySpecs.SOC_min)* systemParams.INVERTER_EFFICIENCY * batterySpecs.efficiency;
sc_available = max(0, simState.SOC_sc(hour-1) - supercapSpecs.SOC_min) * systemParams.INVERTER_EFFICIENCY * supercapSpecs.efficiency;
% bat_discharge_target=bat_discharge_target/systemParams.INVERTER_EFFICIENCY/batterySpecs.efficiency;
% sc_discharge_target=sc_discharge_target/systemParams.INVERTER_EFFICIENCY/supercapSpecs.efficiency;
% Calculate actual discharging power
simState.P_dh_bat(hour) = min([
    bat_discharge_target,
    bat_available,
    batterySpecs.max_discharge * systemParams.TIME_STEP_HOURS
    ]);

simState.P_dh_sc(hour) = min([
    sc_discharge_target,
    sc_available,
    supercapSpecs.max_discharge * systemParams.TIME_STEP_HOURS
    ]);

% Calculate remaining deficit (load served by grid)
total_deficit = bat_discharge_target + sc_discharge_target;
total_supplied = simState.P_dh_bat(hour) + simState.P_dh_sc(hour);
simState.LSBG(hour) = max(0, total_deficit - total_supplied);

% Reset charging variables
simState.P_ch_bat(hour) = 0;
simState.P_ch_sc(hour) = 0;
simState.dumped_power(hour) = 0;
simState.u_state(hour) = 1;  % Discharging state
end

%% State Update Functions

function simState = updateSelfDischargeLosses(simState, hour, batterySpecs, supercapSpecs)
% Calculate self-discharge losses for both ESS types

% Battery self-discharge
if simState.idle_cnt_bat(hour-1) > 0
    discharge_factor = 1 - (1 - batterySpecs.self_discharge_rate)^simState.idle_cnt_bat(hour-1);
else
    discharge_factor = batterySpecs.self_discharge_rate;
end
simState.self_discharge_loss_bat(hour) = simState.SOC_bat(hour-1) * discharge_factor;

% Supercapacitor self-discharge
if simState.idle_cnt_sc(hour-1) > 0
    discharge_factor = 1 - (1 - supercapSpecs.self_discharge_rate)^simState.idle_cnt_sc(hour-1);
else
    discharge_factor = supercapSpecs.self_discharge_rate;
end
simState.self_discharge_loss_sc(hour) = simState.SOC_sc(hour-1) * discharge_factor;

% Apply self-discharge losses to previous SOC
simState.SOC_bat(hour-1) = simState.SOC_bat(hour-1) - simState.self_discharge_loss_bat(hour);
simState.SOC_sc(hour-1) = simState.SOC_sc(hour-1) - simState.self_discharge_loss_sc(hour);
end

function simState = updateSOCStates(simState, hour, batterySpecs, supercapSpecs)
% Update SOC states and enforce limits

% Update SOC based on power flows
simState.SOC_bat(hour) = simState.SOC_bat(hour-1) + (simState.P_ch_bat(hour) - simState.P_dh_bat(hour));
simState.SOC_sc(hour) = simState.SOC_sc(hour-1) + (simState.P_ch_sc(hour) - simState.P_dh_sc(hour));

% Enforce SOC limits for battery
if simState.SOC_bat(hour) < batterySpecs.SOC_min
    simState.LSBG(hour) = simState.LSBG(hour) + (batterySpecs.SOC_min - simState.SOC_bat(hour));
    simState.SOC_min_bat_violations = simState.SOC_min_bat_violations + 1;
    simState.SOC_bat(hour) = batterySpecs.SOC_min;
elseif simState.SOC_bat(hour) > batterySpecs.SOC_max
    simState.dumped_power(hour) = simState.dumped_power(hour) + (simState.SOC_bat(hour) - batterySpecs.SOC_max);
    simState.SOC_bat(hour) = batterySpecs.SOC_max;
end

% Enforce SOC limits for supercapacitor
if simState.SOC_sc(hour) < supercapSpecs.SOC_min
    simState.LSBG(hour) = simState.LSBG(hour) + (supercapSpecs.SOC_min - simState.SOC_sc(hour));
    simState.SOC_min_sc_violations = simState.SOC_min_sc_violations + 1;
    simState.SOC_sc(hour) = supercapSpecs.SOC_min;
elseif simState.SOC_sc(hour) > supercapSpecs.SOC_max
    simState.dumped_power(hour) = simState.dumped_power(hour) + (simState.SOC_sc(hour) - supercapSpecs.SOC_max);
    simState.SOC_sc(hour) = supercapSpecs.SOC_max;
end
end

function simState = updateIdleCounters(simState, hour)
% Update idle counters for both ESS types

% Battery idle counter
if simState.P_ch_bat(hour) == 0 && simState.P_dh_bat(hour) == 0
    simState.idle_cnt_bat(hour) = simState.idle_cnt_bat(hour-1) + 1;
else
    simState.idle_cnt_bat(hour) = 0;
end

% Supercapacitor idle counter
if simState.P_ch_sc(hour) == 0 && simState.P_dh_sc(hour) == 0
    simState.idle_cnt_sc(hour) = simState.idle_cnt_sc(hour-1) + 1;
else
    simState.idle_cnt_sc(hour) = 0;
end
end

function simState = updateSystemMetrics(simState, hour, powerProfiles, systemParams)
% Update system-level performance metrics

% Calculate net power flows
simState.P_net(hour) = (powerProfiles.load(hour) + simState.P_ch_bat(hour) + simState.P_ch_sc(hour)) - ...
    (simState.P_dh_bat(hour) + simState.P_dh_sc(hour) + powerProfiles.total_renewable(hour));
%
% simState.P_net_noGrid(hour) = powerProfiles.load(hour) - ...
%     (simState.P_dh_bat(hour) + simState.P_dh_sc(hour) + powerProfiles.total_renewable(hour));

simState.P_net_noGrid(hour) = powerProfiles.load(hour) -(powerProfiles.total_renewable(hour)-((simState.P_ch_bat(hour)/ systemParams.INVERTER_EFFICIENCY / systemParams.batteff) + (simState.P_ch_sc(hour)/ systemParams.INVERTER_EFFICIENCY / systemParams.supercapeff))) ...
    -(simState.P_dh_bat(hour) + simState.P_dh_sc(hour));

% Calculate peak load following metric
simState.PLF_seq(hour) = simState.P_net_noGrid(hour) - simState.P_net_noGrid(hour-1);


% --- Insert energy balance calculation ---
simState.energy_balance(hour) = ...
    powerProfiles.load(hour) ...
    - powerProfiles.total_renewable(hour) ...
    - powerProfiles.geothermal(hour) ...
    - (simState.P_dh_bat(hour) + simState.P_dh_sc(hour)) ...
    + (simState.P_ch_bat(hour) / systemParams.INVERTER_EFFICIENCY / systemParams.batteff) ...
    + (simState.P_ch_sc(hour) / systemParams.INVERTER_EFFICIENCY / systemParams.supercapeff) ...
    - simState.LSBG(hour) ...
    - simState.dumped_power(hour) ...
    + simState.self_discharge_loss_bat(hour) ...
    + simState.self_discharge_loss_sc(hour);
end

%% Aging and Degradation Functions

function simState = performWeeklyAging(simState, hour, weeklyIndices, batterySpecs, systemParams)
% Perform weekly battery aging calculations

week_idx = ceil(hour / systemParams.BDCF);
week_hours = weeklyIndices{week_idx};

% Calculate normalized SOC for the week
SOC_normalized = (simState.SOC_bat(week_hours) - batterySpecs.SOC_min) / ...
    (batterySpecs.SOC_max - batterySpecs.SOC_min);
SOC_normalized = min(max(SOC_normalized, 0), 1);

% Calculate cycling degradation using rainflow counting
cycling_degradation = calculateCyclingDegradation(SOC_normalized, batterySpecs.aging_coeffs);

% Calculate calendar degradation
avg_SOC = mean(SOC_normalized);
calendar_degradation = calculateCalendarDegradation(avg_SOC, batterySpecs.aging_coeffs, length(weeklyIndices));

% Update degradation history
weekly_degradation = calendar_degradation + cycling_degradation;
simState.deg_weekly(week_idx) = weekly_degradation;

if week_idx == 1
    simState.cum_deg_hist(week_idx) = weekly_degradation;
else
    simState.cum_deg_hist(week_idx) = simState.cum_deg_hist(week_idx-1) + weekly_degradation;
end
simState.cum_deg_hist(week_idx) = min(simState.cum_deg_hist(week_idx), 1);

% Update capacity and SOC limits
remaining_capacity = (1 - simState.cum_deg_hist(week_idx)) * batterySpecs.initial_capacity;
simState.cap_hist(week_idx) = remaining_capacity;

% Check end-of-life criterion (75% of initial capacity)
EOL_threshold = 0.75 * batterySpecs.initial_capacity;
if remaining_capacity < EOL_threshold
    remaining_capacity = EOL_threshold;
    simState.EOL_violations = simState.EOL_violations + 1;
end

% Update SOC limits based on degraded capacity
batterySpecs.SOC_max = remaining_capacity * (batterySpecs.specs.SOC_max/100);
batterySpecs.SOC_min = remaining_capacity * (batterySpecs.specs.SOC_min/100);

% Calculate degradation cost
capacity_loss = simState.cum_deg_hist(week_idx) * batterySpecs.initial_capacity;
simState.degradation_cost = capacity_loss * batterySpecs.specs.installation_cost;
end

function cycling_deg = calculateCyclingDegradation(SOC_normalized, aging_coeffs)
% Calculate cycling degradation using rainflow counting

cycles = rainflow(SOC_normalized * 100);
if ~isempty(cycles)
    cycle_counts = cycles(:,1);
    depth_of_discharge = cycles(:,2)/100;
    cycling_deg = sum((aging_coeffs.A_cyc * depth_of_discharge.^2 + aging_coeffs.B_cyc * depth_of_discharge) .* cycle_counts);
else
    cycling_deg = 0;
end
end

function calendar_deg = calculateCalendarDegradation(avg_SOC, aging_coeffs, num_weeks)
% Calculate calendar (idle) degradation

calendar_deg = (aging_coeffs.A_idl * avg_SOC^2 + aging_coeffs.B_idl * avg_SOC + aging_coeffs.C_idl) / num_weeks;
end

%% Cost and Emissions Calculation

function [Total_costs, costs,grid_costs,penalties,PSC,Total_emissions] = calculateFinalCostsAndEmissions(results, systemParams, config,powerProfiles)
% Calculate final lifecycle costs and emissions

% Calculate component lifecycle costs
costs = calculateComponentCosts(results, systemParams, config);

% Calculate grid costs
grid_costs = calculateGridCosts(results, systemParams, config);

% Calculate penalties
penalties = calculatePenalties(results, systemParams, config);
PSC       = calculatepeakshavingfactor(results, systemParams, config);

% Calculate emissions
emissions = calculateEmissions(results, systemParams, config,powerProfiles);

% Sum total costs
Total_costs = costs.capital + costs.replacement + costs.maintenance + ...
    grid_costs + results.degradation_cost + penalties;

Total_emissions = emissions.total;
% simState.costs=costs;
% simState.gridcosts=grid_costs;
% simState.penaltiescosts=penalties;
% simState.PeakShavcosts=PSC;

% Putting multiple objectives in array
f=[Total_costs, PSC];
end

function costs = calculateComponentCosts(results, systemParams, config)
% Calculate lifecycle costs for all system components

% Component costs
solar_costs = computePVLifecycleCosts(config.N_PV, systemParams.PROJECT_LIFE_YEARS, ...
    systemParams.INTEREST_RATE, systemParams.INFLATION_RATE);
wind_costs = computeWindLifecycleCosts(config.N_WT, systemParams.PROJECT_LIFE_YEARS, ...
    systemParams.INTEREST_RATE, systemParams.INFLATION_RATE);
battery_costs = computeESSLifecycleCosts(config.N_BT, 100, config.battery_type, ...
    systemParams.PROJECT_LIFE_YEARS, systemParams.INTEREST_RATE, systemParams.INFLATION_RATE);
supercap_costs = computeESSLifecycleCosts(config.N_SC, 10, 'SC', ...
    systemParams.PROJECT_LIFE_YEARS, systemParams.INTEREST_RATE, systemParams.INFLATION_RATE);
solar_power=config.power_profiles.solar.data;
wind_power=config.power_profiles.wind.data;
GT_power=config.power_profiles.geothermal.data;
P_BESS=results.P_ch_bat+results.P_dh_bat;
P_SCESS=results.P_ch_sc+results.P_dh_sc;
LSBG=results.LSBG;
Pr_life=systemParams.PROJECT_LIFE_YEARS;
int_rate=systemParams.INTEREST_RATE;        % 3.4%
inf_rate=systemParams.INFLATION_RATE;       % 2.4%;
[conv_costs, trans_costs] = computeConverterTransformerCosts( ...
    solar_power,wind_power , GT_power, P_BESS,P_SCESS, LSBG, Pr_life, int_rate, inf_rate);

% Sum costs
costs.capital = solar_costs.capital + wind_costs.capital + battery_costs.capital + supercap_costs.capital+conv_costs.capital+trans_costs.capital;
costs.replacement = sum(solar_costs.replacement) + sum(wind_costs.replacement) + ...
    sum(battery_costs.replacement) + sum(supercap_costs.replacement)+ sum(conv_costs.replacement)   + sum(trans_costs.replacement);;
costs.maintenance = solar_costs.maintenance + wind_costs.maintenance + ...
    battery_costs.maintenance + supercap_costs.maintenance+ conv_costs.maintenance     + trans_costs.maintenance;
costs.solar_costs=solar_costs.capital+sum(solar_costs.replacement)+solar_costs.maintenance;
costs.batt_costs=battery_costs.capital+sum(battery_costs.replacement)+battery_costs.maintenance;
costs.sc_costs=supercap_costs.capital+sum(supercap_costs.replacement)+supercap_costs.maintenance;
costs.batt_ag_cost=results.degradation_cost;
costs.convrtr_costs=conv_costs;

costs.transformer_costs=trans_costs;
end

function grid_costs = calculateGridCosts(results, systemParams, config)
% Calculate grid energy costs based on facility

EURO_TO_USD = 1.09;
annual_net_load = sum(results.P_net_noGrid);
peak_load = max(results.P_net_noGrid);
total_grid_energy = sum(results.LSBG);

switch upper(config.facility)
    case 'case1'
        grid_price = case1_GRID_PRICE(annual_net_load, peak_load, config.do_PLC);
        grid_costs = EURO_TO_USD * grid_price * total_grid_energy;
    case 'case3'
        grid_price = case3_GRID_PRICE(total_grid_energy);
        grid_costs = EURO_TO_USD * grid_price ;
    case 'case2'
        grid_price = case2_GRID_PRICE(total_grid_energy);
        grid_costs = EURO_TO_USD * grid_price ;
    otherwise
        error('Unknown facility: %s', config.facility);
end


end
function PSC=calculatepeakshavingfactor(results, systemParams, config)
% Peak load following penalty
if config.do_NLPS
    PSC = sum(0.5 * (results.PLF_seq).^2);
else
    PSC=0;
end
end

function penalties = calculatePenalties(results, systemParams, config)
% Calculate system penalties

penalties = 0;

% Energy balance penalty (should be near zero for feasible solutions)
energy_balance_error = sum(round(abs(results.energy_balance),4));
if energy_balance_error > 1e-3
    penalties = penalties + 1e12;  % Large penalty for infeasible solutions
end
end

function emissions = calculateEmissions(results, systemParams, config,powerProfiles)
% Calculate total CO2 emissions

% Emission factors (gCO2/kWh)
GRID_CO2 = 329; %gCO2/kWh %https://www.eea.europa.eu/en/analysis/indicators/greenhouse-gas-emission-intensity-of-1/greenhouse-gas-emission-intensity-of-electricity-generation-country-level
PV_CO2 = getPVModuleLibrary().Mono_SI.Embedded_CO2;
WIND_CO2 = getWindTurbineLibrary().wt27.embedded_CO2;
BATTERY_CO2 = getESSLibrary().(config.battery_type).Embedded_CO2 * 1000;  % kg to g
SC_CO2 = getESSLibrary().SC.Embedded_CO2 * 1000;  % kg to g
TR_CO2=6955.40;% kg CO2 eq/MVA Assessing the Life Cycle CO2 Emissions of Various Transformers Under Different Scenarios with Multiple Functional Units
TR_CO2=6955.40/1000*1000; %g CO2 eq/kVA


wind_power=config.power_profiles.wind.data;
TR_rating=((max(wind_power)/0.8) * 1.20)+20; %kVA, wind + grid transformer
% Here GWP is provided per kVA of Transformer total capcaity
emissions.transformers=TR_rating*TR_CO2; %g CO2 eq


% Calculate emissions by source
emissions.grid = sum(results.LSBG) * GRID_CO2;
% Here GWP is provided per kWh of enrgy geenration by wind plant
emissions.pv = sum(powerProfiles.solar) * PV_CO2;  % Simplified
% Here GWP is provided per kWh of enrgy geenration by wind plant
emissions.wind = sum(powerProfiles.wind) * WIND_CO2;
% Find toatl batteries useable cacpaity(kWh) and then multiply GWP(g CO-eq/kWh of useable capcaity)
% Battery system parameters
[batterySpecs, supercapSpecs] = initializeESSSpecifications(config, systemParams);
emissions.battery = (config.N_BT*batterySpecs.rated_cap_unit*(getESSLibrary().(config.battery_type).SOC_max-getESSLibrary().(config.battery_type).SOC_min)/100)* BATTERY_CO2;
% Here GWP is provided per kWh of total enrgy passed through supercap
emissions.supercap = sum(results.P_ch_sc + results.P_dh_sc) * SC_CO2;

emissions.total = emissions.grid + emissions.pv + emissions.wind + emissions.battery+emissions.supercap;
end
function coeff = getBatteryAgingCoefficients(battery_type)
coefficients = struct( ...
    'LFP', struct('A_idl', 6.02e-06, 'B_idl', 1.35e-05, 'C_idl', 1.85e-05, ...
    'A_cyc', -4.72e-05, 'B_cyc', 9.62e-05), ...
    'LMO', struct('A_idl', 6.81e-05, 'B_idl', 4.02e-05, 'C_idl', 1.63e-05, ...
    'A_cyc', -1.21e-04, 'B_cyc', 4.01e-04), ...
    'NMC', struct('A_idl', 8.07e-06, 'B_idl', 3.41e-06, 'C_idl', 2.83e-05, ...
    'A_cyc', -4.05e-05, 'B_cyc', 1.01e-04), ...
    'LTO', struct('A_idl', 3.03e-06, 'B_idl', 2.81e-05, 'C_idl', 5.02e-06, ...
    'A_cyc', -1.57e-05, 'B_cyc', 4.40e-05) ...
    );
coeff = coefficients.(battery_type);
end

%% -----------------------------------------------------------------------------
%  Helper: PV lifecycle costs (calls SolarPowerCalculatorv2.calculateLifecycleCosts)
% -----------------------------------------------------------------------------

function costs = computePVLifecycleCosts(N_PV, Pr_life, i_rate, f_rate)
spc = SolarPowerCalculatorv2('Mono_SI_SOLARIA');
costs = spc.calculateLifecycleCosts(N_PV, Pr_life, i_rate, f_rate);
end


%% -----------------------------------------------------------------------------
%  Helper: Wind lifecycle costs
% -----------------------------------------------------------------------------
function costs = computeWindLifecycleCosts(N_wt, Pr_life, i_rate, f_rate)
wtc = WindTurbinePowerCalculator('wt27');
costs = wtc.calculateLifecycleCosts(N_wt, Pr_life, i_rate, f_rate);
end


%% -----------------------------------------------------------------------------
%  Helper: Geothermal lifecycle costs
% -----------------------------------------------------------------------------
function costs = computeGeothermalLifecycleCosts(N_wells, Pr_life, i_rate, f_rate)
if N_wells <= 0
    costs.capital     = 0;
    costs.replacement = 0;
    costs.maintenance = 0;
    return;
end
GTGen = createGeneratorsLibrary().hydrothermal_geothermal;
well_depth_m = GTGen.well_depth * 1000;
GTcap = N_wells * 100;
costs = calculateGTLifecycleCostsv1(N_wells, well_depth_m, GTcap, Pr_life, i_rate, f_rate);
end


%% -----------------------------------------------------------------------------
%  Helper: ESS lifecycle costs
% -----------------------------------------------------------------------------
function costs = computeESSLifecycleCosts(N_ESS, rated_cap_unit, ESS_type, Pr_life, i_rate, f_rate)
ESSSpecs = getESSLibrary().(ESS_type);
% rated_ESS_energy = rated_cap_unit * N_BT;
costs = calculateESSLifecycleCosts(rated_cap_unit, N_ESS, ...
    ESSSpecs.calendar_life, ESSSpecs.installation_cost, Pr_life, i_rate, f_rate);
end


%% -----------------------------------------------------------------------------
%  Helper: converter + transformer lifecycle costs
% -----------------------------------------------------------------------------
function [conv_costs, trans_costs] = computeConverterTransformerCosts( ...
    solar_series, wind_series, GT_power, P_BESS,P_SCESS, LSBG, Pr_life, i_rate, f_rate)

% 1) Sizing (kW)
PC_PV   = max(solar_series) * 1.20;
PC_WT   = 0;
PC_GT   = 0;
PC_BT   = max(P_BESS)   * 1.20;
PC_SC   = max(P_SCESS)   * 1.20;
PC_KARA = 0;
PC_Grid = 0;

% 2) Converter lifecycle
conv_costs.PV   = calculateConLifecycleCosts(PC_PV,   Pr_life, i_rate, f_rate);
conv_costs.WT   = calculateConLifecycleCosts(PC_WT,   Pr_life, i_rate, f_rate);
conv_costs.GT   = calculateConLifecycleCosts(PC_GT,   Pr_life, i_rate, f_rate);
conv_costs.BT   = calculateConLifecycleCosts(PC_BT,   Pr_life, i_rate, f_rate);
conv_costs.SC   = calculateConLifecycleCosts(PC_SC,   Pr_life, i_rate, f_rate);
conv_costs.KARA = calculateConLifecycleCosts(PC_KARA, Pr_life, i_rate, f_rate);
conv_costs.Grid = calculateConLifecycleCosts(PC_Grid, Pr_life, i_rate, f_rate);
conv_costs.capital=conv_costs.PV.capital+conv_costs.WT.capital+conv_costs.GT.capital+conv_costs.BT.capital+conv_costs.SC.capital+conv_costs.KARA.capital+conv_costs.Grid.capital;
conv_costs.replacement=conv_costs.PV.replacement + conv_costs.WT.replacement + conv_costs.GT.replacement + conv_costs.BT.replacement+ conv_costs.SC.replacement + conv_costs.KARA.replacement + conv_costs.Grid.replacement;
conv_costs.maintenance=conv_costs.PV.maintenance + conv_costs.WT.maintenance + conv_costs.GT.maintenance + conv_costs.BT.maintenance+conv_costs.SC.maintenance + conv_costs.KARA.maintenance + conv_costs.Grid.maintenance;

% 3) Transformer lifecycle (80% power factor)
WTTr   = (max(wind_series)/0.8) * 1.20;
GTTr   = (max(GT_power)/0.8)   * 1.20;
GridTr = 20;  % kVA
trans_costs.WT   = calculateTrLifecycleCosts(WTTr,   Pr_life, i_rate, f_rate);
trans_costs.GT   = calculateTrLifecycleCosts(GTTr,   Pr_life, i_rate, f_rate);
trans_costs.Grid = calculateTrLifecycleCosts(GridTr,  Pr_life, i_rate, f_rate);
trans_costs.capital = trans_costs.WT.capital+trans_costs.GT.capital+trans_costs.Grid.capital;
trans_costs.replacement = trans_costs.WT.replacement+trans_costs.GT.replacement+trans_costs.Grid.replacement;
trans_costs.maintenance = trans_costs.WT.maintenance+trans_costs.GT.maintenance+trans_costs.Grid.maintenance;

% zero out unused
trans_costs.PV   = struct('capital',0,'replacement',0,'maintenance',0);
trans_costs.BT   = struct('capital',0,'replacement',0,'maintenance',0);
trans_costs.KARA = struct('capital',0,'replacement',0,'maintenance',0);
end
