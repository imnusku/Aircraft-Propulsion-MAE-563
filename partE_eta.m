clc; clear; clear all;

% Initialize altitude and Mach number ranges
z_range = 2000:500:20000; % Altitudes in meters
M1_range = 0.8:0.1:5.0;   % Flight Mach numbers

% Preallocate arrays for results
optimal_M1_efficiency = zeros(size(z_range));
optimal_M1_TSFC = zeros(size(z_range));
max_efficiency_results = zeros(size(z_range));
min_TSFC_results = zeros(size(z_range));

for z_idx = 1:length(z_range)
    z = z_range(z_idx); % Set altitude for this iteration

    % Initialize temporary storage for results
    efficiency_for_M1 = zeros(size(M1_range));
    TSFC_for_M1 = zeros(size(M1_range));
    
    for M1_idx = 1:length(M1_range)
        M1 = M1_range(M1_idx); % Set Mach number for this iteration
        
        %% module 1
        % Constants
        gamma = 1.4; % Specific heat ratio
        R = 286.9; % Specific gas constant in J/(kg*K)
        z_star = 8404; % Scale height in meters
        T_s = 288.0; % Standard temperature at sea level in Kelvin
        p_s = 101.3; % Standard pressure at sea level in kPa

        % Calculate T and P for each altitude
        if z < 7958 % Within the troposphere
            T1 = T_s * (1 - (((gamma - 1) / gamma) * (z / z_star)));
            p1 = p_s * ((1 - (((gamma - 1) / gamma) * (z / z_star)))^(gamma / (gamma - 1)));
        else % In the tropopause
            T1 = 210.0; % Constant temperature in tropopause
            p1 = 33.6 * exp(-(z - 7958) / 6605);
        end

        % Total-to-static relations
        Tt1 = T1 * (1 + (gamma - 1) / 2 * M1^2);
        pt1 = p1 * (1 + (gamma - 1) / 2 * M1^2)^(gamma / (gamma - 1));

        % Sound speed and velocity
        a1 = sqrt(gamma * R * T1);
        V1 = M1 * a1;

        %% module 2
        gamma = 1.4; % Specific heat ratio
        eta_d = 0.92; % Inlet/diffuser efficiency
        M2 = 0.15;

        Tt2 = Tt1;
        T2 = Tt2 / (1 + (gamma - 1) / 2 * M2^2);
        pt2 = p1 * (1 + (eta_d * (gamma - 1) / 2) * M1^2)^(gamma / (gamma - 1));
        p2 = pt2 / (1 + (gamma - 1) / 2 * M2^2)^(gamma / (gamma - 1));

        %% module 3
        gamma = 1.3; % Changed after the combustor
        Tt3_max = 2400;

        % Choking check
        Tt3_choked = Tt2 * (1 / (2 * (gamma + 1))) * ...
            ((1 / (M2^2) * ((1 + gamma * M2^2)^2))) * ...
            ((1 + (((gamma - 1) / 2) * M2^2)))^(-1);

        if Tt3_choked < Tt3_max
            Tt3 = Tt3_choked;
            M3 = 1;
        else
            Tt3 = Tt3_max;
            C = (Tt3 / Tt2) * ((1 + (gamma - 1) / 2 * M2^2) / ((1 + gamma * M2^2)^2)) * M2^2;
            a = C * (gamma^2) - ((gamma - 1) / 2);
            b = 2 * C * gamma - 1;
            c = C;
            M3_roots = roots([a, b, c]);
            M3_squared = M3_roots(M3_roots > 0 & M3_roots < 1);
            M3 = sqrt(M3_squared);
        end

        q23 = 986 * (Tt3 - Tt2) + 0.5 * 0.179 * (Tt3^2 - Tt2^2);
        T3 = Tt3 ./ (1 + (gamma - 1) ./ 2 * M3.^2);
        p3 = p2;
        pt3 = p3 * ((1 + (gamma - 1) / 2 * M3.^2).^(gamma / (gamma - 1)));

        %% module 4
        Ae = 0.015;
        eta_n = 0.94;
        Tte = Tt3;

        test_M = sqrt(2 / (gamma - 1)) * sqrt((eta_n * ...
            ((1 - (p1 / pt3)^((gamma - 1) / gamma))) / ...
            (1 - eta_n * (1 - (p1 / pt3)^((gamma - 1) / gamma)))));

        if test_M < 1
            Me = test_M;
            pe = p1;
        else
            Me = 1;
            pe = pt3 * (1 - (1 / eta_n) * (gamma - 1) / (gamma + 1))^(gamma / (gamma - 1));
        end

        Te = Tte / (1 + (gamma - 1) / 2 * Me^2);
        ae = sqrt(gamma * R * Te);
        Ve = Me * ae;
        rho_e = pe * 1000 / (R * Te);
        m_dot_e = rho_e * Ve * Ae;

        %% module 6
        g0 = 9.81;
        m_dot_i = m_dot_e / (1 + q23 / 43.2e6);
        m_dot_f = m_dot_e - m_dot_i;
        f = m_dot_f / m_dot_i;
        jet_thrust = m_dot_i .* (1 + f) .* Ve - m_dot_i .* V1;
        pressure_thrust = (pe - p1) * 1000 * Ae;
        total_thrust = jet_thrust + pressure_thrust;
        Veq = Ve + ((pe - p1) * 1000 * Ae / m_dot_e);
        TSFC = (m_dot_f ./ total_thrust) * 3600;
        thermal_efficiency = (m_dot_e .* Veq^2 ./ 2 - m_dot_i .* V1^2 ./ 2) ./ (m_dot_i .* q23);
        propulsive_efficiency = 2 ./ (1 + (Veq ./ V1));
        overall_efficiency = thermal_efficiency * propulsive_efficiency;

        %% Store temporary results
        efficiency_for_M1(M1_idx) = overall_efficiency;
        TSFC_for_M1(M1_idx) = TSFC;
    end

    %% Find optimal M1 for this altitude
    [max_efficiency, max_eff_idx] = max(efficiency_for_M1);
    [min_TSFC, min_TSFC_idx] = min(TSFC_for_M1);

    % Store optimal results
    optimal_M1_efficiency(z_idx) = M1_range(max_eff_idx);
    optimal_M1_TSFC(z_idx) = M1_range(min_TSFC_idx);
    max_efficiency_results(z_idx) = max_efficiency;
    min_TSFC_results(z_idx) = min_TSFC;
end

%% Plot the results
figure;
plot(z_range, optimal_M1_efficiency, 'LineWidth', 2);
xlabel('Altitude (m)');
ylabel('Optimal Mach Number for Max Efficiency');
title('Optimal Mach Number vs. Altitude for Max Efficiency');
ylim([0 7]); % Set y-axis limits
grid on;

figure;
plot(z_range, optimal_M1_TSFC, 'LineWidth', 2);
xlabel('Altitude (m)');
ylabel('Optimal Mach Number for Min TSFC');
title('Optimal Mach Number vs. Altitude for Min TSFC');
ylim([0 7]); % Set y-axis limits
grid on;

