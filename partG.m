clc; clear; clear all;

% Initialize nozzle efficiency range
eta_n_range = 0.5:0.05:1.0;

% Preallocate arrays for results
overall_efficiency_results = zeros(size(eta_n_range));
thrust_results = zeros(size(eta_n_range));
TSFC_results = zeros(size(eta_n_range));

% Fixed parameters
z = 4300; % Altitude in meters
M1 = 2.4; % Mach number
eta_d = 0.92; % Fixed inlet/diffuser efficiency

for idx = 1:length(eta_n_range)
    eta_n = eta_n_range(idx); % Current nozzle efficiency

    %% module 1
    % Constants
    gamma = 1.4; % Specific heat ratio
    R = 286.9; % Specific gas constant in J/(kg*K)
    z_star = 8404; % Scale height in meters
    T_s = 288.0; % Standard temperature at sea level in Kelvin
    p_s = 101.3; % Standard pressure at sea level in kPa

    % Calculate T and P for the given altitude
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

    %% Store results
    overall_efficiency_results(idx) = overall_efficiency;
    thrust_results(idx) = total_thrust;
    TSFC_results(idx) = TSFC;
end

%% Plot the results
figure;
plot(eta_n_range, overall_efficiency_results, 'LineWidth', 2);
xlabel('Nozzle Efficiency (\eta_n)');
ylabel('Overall Efficiency');
title('Overall Efficiency vs. Nozzle Efficiency');
grid on;

figure;
plot(eta_n_range, thrust_results, 'LineWidth', 2);
xlabel('Nozzle Efficiency (\eta_n)');
ylabel('Thrust (N)');
title('Thrust vs. Nozzle Efficiency');
grid on;

figure;
plot(eta_n_range, TSFC_results, 'LineWidth', 2);
xlabel('Nozzle Efficiency (\eta_n)');
ylabel('TSFC (kg/hr/N)');
title('TSFC vs. Nozzle Efficiency');
grid on;
