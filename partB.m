clc;clear;clear all;
%% module 1

% Constants
gamma = 1.4;        % Specific heat ratio
R = 286.9;          % Specific gas constant in J/(kg*K)
z_star = 8404;      % Scale height in meters
T_s = 288.0;        % Standard temperature at sea level in Kelvin
p_s = 101.3;        % Standard pressure at sea level in kPa

% Input Parameters
z = linspace(0, 22000, 500);  % Altitude range from 0 to 25,000 meters

% Initialize arrays for temperature and pressure
T1 = zeros(size(z));
p1 = zeros(size(z));

% Calculate T and P for each altitude
for i = 1:length(z)
    if z(i) < 7958  % Within the troposphere
        T1(i) = T_s * (1 - (((gamma - 1) / gamma) * (z(i) / z_star)));
        p1(i) = p_s * ((1 - (((gamma - 1) / gamma) * (z(i) / z_star)))^(gamma / (gamma - 1)));
    else             % In the tropopause
        T1(i) = 210.0;  % Constant temperature in tropopause
        p1(i) = 33.6 * exp(-(z(i) - 7958) / 6605);
    end
end

% International Standard Atmosphere (ISA) data
z_isa = [0, 500, 1000, 1500, 2000, 2500, 3000, 3500, 4000, 4500, 5000, ...
         5500, 6000, 6500, 7000, 7500, 8000, 8500, 9000, 9500, 10000, ...
         10500, 11000, 11500, 12000, 12500, 13000, 13500, 14000, ...
         14500, 15000, 15500, 16000, 16500, 17000, 17500, 18000, ...
         18500, 19000, 19500, 20000, 22000];
T_isa = [288.15, 284.9, 281.7, 278.4, 275.2, 271.9, 268.7, 265.4, ...
         262.2, 258.9, 255.7, 252.4, 249.2, 245.9, 242.7, 239.5, ...
         236.2, 233, 229.7, 226.5, 223.3, 220, 216.8, 216.7, ...
         216.7, 216.7, 216.7, 216.7, 216.7, 216.7, 216.7, 216.7, ...
         216.7, 216.7, 216.7, 216.7, 216.7, 216.7, 216.7, 216.7, 216.7 218.6];
p_isa = [101.325, 95.46, 89.88, 84.56, 79.5, 74.69, 70.12, 65.78, ...
         61.66, 57.75, 54.05, 50.54, 47.22, 44.08, 41.11, 38.3, ...
         35.65, 33.15, 30.8, 28.58, 26.5, 24.54, 22.7, 20.98, ...
         19.4, 17.93, 16.58, 15.33, 14.17, 13.1, 12.11, 11.2, ...
         10.35, 9.572, 8.85, 8.182, 7.565, 6.995, 6.467, 5.98, 5.529, 4.047];

% Plot T vs z
figure;
plot(T1, z, 'b-', 'LineWidth', 2); % Modeled data
hold on;
plot(T_isa, z_isa, 'r--', 'LineWidth', 2); % ISA data
ylabel('Altitude, z (m)');
zlabel('Static Temperature, T (K)');
title('Temperature vs Altitude');
legend('Modeled Data', 'ISA Data');
grid on;

% Plot P vs z
figure;
plot(p1,z, 'b-', 'LineWidth', 2); % Modeled data
hold on;
plot(p_isa, z_isa, 'r--', 'LineWidth', 2); % ISA data
ylabel('Altitude, z (m)');
xlabel('Static Pressure, P (kPa)');
title('Pressure vs Altitude');
legend('Modeled Data', 'ISA Data');
grid on;
