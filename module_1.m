clc;clear;clear all;
%% module 1

% Constants
gamma = 1.4;        % Specific heat ratio
R = 286.9;          % Specific gas constant in J/(kg*K)
z_star = 8404;      % Scale height in meters
T_s = 288.0;        % Standard temperature at sea level in Kelvin
p_s = 101.3;        % Standard pressure at sea level in kPa

% Input Parameters
z = 10000;
M1=0.82;

% Calculate T and P for each altitude
    %if z < 7958  % Within the troposphere
        T1 = T_s * (1 - (((gamma - 1) / gamma) * (z / z_star)));
        p1 = p_s * ((1 - (((gamma - 1) / gamma) * (z / z_star)))^(gamma / (gamma - 1)));
    %else             % In the tropopause
     %   T1 = 210.0;  % Constant temperature in tropopause
      %  p1 = 33.6 * exp(-(z - 7958) / 6605);
    %end

% Total-to-static relations
Tt1 = T1 * (1 + (gamma - 1) / 2 * M1^2);
pt1 = p1 * (1 + (gamma - 1) / 2 * M1^2)^(gamma / (gamma - 1));

% Sound speed and velocity
a1 = sqrt(gamma * R * T1);
V1 = M1 * a1;
%% module 2

% Constants
gamma = 1.4;        % Specific heat ratio
R = 286.9;          % Specific gas constant in J/(kg*K)

% Inputs from Module 1 or given data
%M2 = 0.15;          % Mach number at State 2
M2=0.4;              %vary for case 2
eta_d = 0.97;       % inlet/diffuser efficiency

% Module 2 Calculations
% Total temperature remains constant (no work or heat transfer)
Tt2 = Tt1;

% Compute static temperature at State 2
T2 = Tt2 / (1 + (gamma - 1) / 2 * M2^2);

% Compute total and static pressures
pt2 = p1 * (1 + (eta_d *(gamma - 1) / 2) * M1^2)^(gamma / (gamma - 1));
p2 = pt2 / (1 + (gamma - 1) / 2 * M2^2)^(gamma / (gamma - 1));

% Compute entropy change across the diffuser
cp2 = 1004; % Specific heat at constant pressure (J/(kg*K)) for air
Delta_s_12 = cp2* log(Tt2 / Tt1) - R * log(pt2 / pt1); % Entropy change (J/(kg*K))

% Compute speed of sound and velocity at State 2
a2 = sqrt(gamma * R * T2); % Speed of sound at State 2
V2 = M2 * a2;              % Velocity at State 2

% %% module 3
% 
% % Constants
% gamma = 1.3;                 % Specific heat ratio changed to 1.3 from 1.4                    
% Tt3_max = 2400;                % Maximum allowable total temperature at combustor exit (K)
% 
% 
% % Module 3 Calculations
% % Step 1: Check if the combustor is thermally choked
% Tt3_choked = Tt2 * (1 / (2 * (gamma + 1))) * ...
%                    ((1 / (M2^2) * ((1 + gamma* M2^2)^2))) * ...
%                    ((1 + (((gamma - 1) / 2) * M2^2)))^(-1);
% 
% % If thermally choked, adjust the maximum temperature
% if Tt3_choked < Tt3_max
%     Tt3 = Tt3_choked;
%     M3=1;
% else
%     Tt3 = Tt3_max;
%     % Solve for M3 in the non-choked case
%      % Use quadratic equation to solve for M3
%      C = (Tt3 / Tt2) * ((1 + (gamma - 1) / 2 * M2^2) / ((1 + gamma * M2^2)^2)) * M2^2;
% 
% 
%     % Quadratic coefficients
%     a = C * (gamma^2) - ((gamma - 1)/2);
%     b = 2 * C * gamma - 1;
%     c = C;
% 
%    % Solve the quadratic equation
%     M3_roots = roots([a, b, c]);
% 
%         if M2<1 
%         M3_squared = M3_roots(M3_roots > 0 & M3_roots <= 1);
%         else
%         M3_squared = M3_roots( M3_roots >=1);
%         end
%      
% 
%     % Take the square root to find M3
%     M3 = sqrt(M3_squared);
% end
% 
% % Step 2: Compute heat added (q23)
% q23 = 986 * (Tt3 - Tt2)+0.5*0.179*(Tt3^2 - Tt2^2);
% 
% % Step 3: Compute static properties at combustor exit
% T3 = Tt3 / (1 + (gamma - 1) / 2 * M3.^2);
% p3 = p2;  % Static pressure remains the same in constant-pressure combustion
% pt3 = p3 * ((1 + (gamma - 1) / 2 * M3.^2).^(gamma/(gamma-1)));
% 
% 
% % Step 4: Compute speed of sound and velocity at combustor exit
% a3 = sqrt(gamma * R * T3); % Speed of sound at State 3
% V3 = M3 * a3;                             % Velocity at State 3
% 
% %need cp3 from T3
% cp3=986+0.179*T3;
% 
% % Step 5: Compute entropy increase across the combustor
% Delta_s_23 = cp3* log(Tt3 / Tt2) - R * log(pt3 / pt2);
% 
% %entropy chnage 1-3
% Delta_s_31 = cp3* log(Tt3 / Tt1) - R * log(pt3 / pt1);
% 
% %% module 4
% 
% % Inputs 
% Ae=0.015;
% eta_n=0.94;
% Tte = Tt3;                    % Total temperature at nozzle entrance
% 
% 
% % Step 1: Compute test Mach number
% test_M = sqrt(2 / (gamma - 1)) * sqrt((eta_n * ...
%           ((1 - (p1 / pt3)^((gamma - 1) / gamma))) / ...
%           (1 - eta_n * (1- (p1 / pt3)^((gamma - 1) / gamma)))));
% 
% 
% % Step 2: Check choking condition
% if test_M < 1
%     % Nozzle is not choked
%     Me = test_M;
%     pe = p1;  % Exit pressure equals ambient pressure
% else
%     % Nozzle is choked
%     Me = 1;
%     pe = pt3 * (1 - (1 / eta_n) * (gamma - 1) / (gamma + 1))^(gamma / (gamma - 1));
% end
% 
% % Step 3: Compute static properties at nozzle exit
% Te = Tte / (1 + (gamma - 1) / 2 * Me^2); % Static temperature
% pte= pe*((1 + (gamma - 1) / 2 * Me^2)^(gamma / (gamma - 1)));
% 
% %Compute speed of sound and velocity at nozzle exit
% ae = sqrt(gamma * R * Te);              % Speed of sound
% Ve = Me * ae;                           % Velocity
% 
% % Step 4: Compute entropy increase
% cpe = 986 + 0.179 * Te; % Specific heat at exit (J/(kg*K))
% Delta_s_3e = cpe * log(Tte / Tt3) - R * log(pte / pt3);
% 
% 
% % Mass flux
% rho_e = pe *1000/ (R * Te); % Convert Pe to Pascals if needed
% m_dot_e = rho_e * Ve * Ae;
% 
% %% module 5
% Tt4=Tte;
% 
% % Step 2: Apply exit criterion for eta_n_ext
% if test_M < 1
%     eta_n_ext = 1;
% else
%     eta_n_ext = test_M.^ (-0.3);
% end
% 
% % Step 3: Compute T4 (Static Temperature at State 4)
% T4 = Tte * (1 - eta_n_ext * (1 - (p1 / pte)^((gamma - 1) / gamma)));
% 
% % Step 4: Compute M4 (Mach Number at State 4)
% M4 = sqrt((2 / (gamma - 1)) * ((Tt4 / T4) - 1));
% 
% 
% % Step 6: Compute p4 (Static Pressure at State 4)
% p4 = p1; % Assume static pressure matches ambient pressure
% 
% % Step 7: Compute pt4 (Total Pressure at State 4)
% pt4 = p4 * (1 + (gamma - 1) / 2 * M4^2)^(gamma / (gamma - 1));
% 
% % Step 8: Compute velocity at State 4
% a4 = sqrt(gamma * R * T4);    % Speed of sound at State 4
% V4 = M4 * a4;                % Velocity at State 4
% 
% % Step 9: Compute entropy increase across the nozzle
% cp4 = 986 + 0.179 * T4;
% Delta_s_4e = cp4 * log(Tt4 / Tte) - R * log(pt4 / pte);
% Delta_s_41 = Delta_s_12+Delta_s_23+Delta_s_3e+Delta_s_4e;
% 
% %% module 6
% 
% % Constants
% g0 = 9.81;                 % Gravitational acceleration (m/s^2)
% q_f = 43.2e6;              % Heating value of fuel (J/kg)
% 
% 
% % Step 1: Compute air mass flow rate (mi)
% m_dot_i = m_dot_e / (1 + q23 / q_f);
% 
% % Step 2: Compute fuel mass flow rate (mf)
% m_dot_f = m_dot_e - m_dot_i;
% 
% % Step 3: Compute fuel-to-air ratio (f)
% f = m_dot_f / m_dot_i;
% 
% % Step 4: Compute thrust
% jet_thrust = m_dot_i * (1 + f) * Ve - m_dot_i * V1; % Jet thrust
% pressure_thrust = (pe - p1)*1000 * Ae;                   % Pressure thrust
% total_thrust = jet_thrust + pressure_thrust;                % Total thrust
% 
% % Step 5: Compute equivalent velocity (Veq)
% Veq = Ve + ((pe - p1) * 1000*Ae / m_dot_e);            % Equivalent velocity
% 
% % Step 6: Compute TSFC
% TSFC = (m_dot_f / total_thrust) * 3600;                    % TSFC in (kg/hr)/N
% 
% % Step 7: Compute specific impulse
% Isp = total_thrust / (m_dot_f * g0);                       % Specific impulse (s)
% 
% % Step 8: Compute efficiencies
% thermal_efficiency = (m_dot_e * Veq^2 / 2 - m_dot_i * V1^2 / 2) / (m_dot_i * q23);
% propulsive_efficiency = 2 / (1+(Veq /V1));
% overall_efficiency = thermal_efficiency * propulsive_efficiency;
% 
% %Propulsive power
% Prop_power=total_thrust*V1;
% 
% 
