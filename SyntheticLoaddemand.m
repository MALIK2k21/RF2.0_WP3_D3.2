function L1=SyntheticLoaddemand(hours)
% MATLAB code to generate hourly profiles for 8760 hours (one year)
% Number of hours in a year


% Time vectors
hourOfDay = mod((0:hours-1), 24)';
dayOfYear = floor((0:hours-1)/24)';
% Load demand profile (L1) for industrial facility with 3-shift operation
% Average demand: 1500 kW with shift-based variations
L1 = zeros(hours,1);
for i = 1:hours
    hour = hourOfDay(i);
    baseLoad = 1500; % Average load in kW
    if hour >= 6 && hour < 14
        % Shift 1 (6 AM - 2 PM): High demand period
        L1(i) = baseLoad * 1.2 + 50*randn();
    elseif hour >= 14 && hour < 22
        % Shift 2 (2 PM - 10 PM): Medium demand period
        L1(i) = baseLoad + 50*randn();
    else
        % Shift 3 (10 PM - 6 AM): Low demand period
        L1(i) = baseLoad * 0.8 + 50*randn();
    end
end
L1(L1 < 0) =end 0; % Ensure no negative load demand