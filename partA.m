clc; clear; close all;

% Constants
a = 986;                 % Coefficient for cp(T) fit
b = 0.179;                   % Constant coefficient for cp(T) fit
cp_4to1 = 1004;            % Constant cp for state 4->1 (J/kg-K)

% Static temperature (K) and entropy (J/(kg·K)) at each state
%T = [245.9, 512.8, 891, 891, 635]; 
%s = [0, 43.95, 658.95, 670.95, 930.95]; 

T = [245.9, 526.8, 2354, 2087, 1558]; 
s = [0, 43.95, 2152, 2164, 2360]; % 

% Extract points for states
T1 = T(1); T2 = T(2); T3 = T(3); Te = T(4); T4 = T(5);
s1 = s(1); s2 = s(2); s3 = s(3); se = s(4); s4 = s(5);

% Step sizes for entropy
ds_23 = 5; % Smaller step size for 2->3
ds_41 = 5; % Smaller step size for 4->1

% 1 -> 2 (straight line)
T_12 = linspace(T1, T2, 100);
s_12 = linspace(s1, s2, 100);

% Initialize entropy and temperature arrays
s_23 = s2:ds_23:s3; % Generate entropy steps
T_23 = zeros(size(s_23)); % Preallocate temperature array
T_23(1) = T2; % Set initial temperature

% Iteratively compute T values
for i = 2:length(s_23)
    % Compute change in entropy
    ds_actual = s_23(i) - s_23(i-1);
    
    % Update temperature using the relation
    dT = (T_23(i-1) * ds_actual) / (a + b * T_23(i-1));   
    T_23(i) = T_23(i-1) + dT;
end



% 3 -> e (straight line)
T_3e = linspace(T3, Te, 100);
s_3e = linspace(s3, se, 100);

% e -> 4 (straight line)
T_e4 = linspace(Te, T4, 100);
s_e4 = linspace(se, s4, 100);

% 4 -> 1 (constant-pressure curve)
T_41 = T4; % Initialize with T4
s_41 = s4:-ds_41:s1; % Entropy steps
for i = 2:length(s_41)
    T_41(i) = T_41(i-1) - (T_41(i-1) * ds_41) / cp_4to1;
end

% Plot T-s diagram
figure;
hold on;

% Add lines and curves
plot(s_12, T_12, 'b', 'LineWidth', 2, 'DisplayName', '1->2 (Straight)');
plot(s_23, T_23, 'r', 'LineWidth', 2, 'DisplayName', '2->3 (Constant-p)');
plot(s_3e, T_3e, 'g', 'LineWidth', 2, 'DisplayName', '3->e (Straight)');
plot(s_e4, T_e4, 'm', 'LineWidth', 2, 'DisplayName', 'e->4 (Straight)');
plot(s_41, T_41, 'c', 'LineWidth', 2, 'DisplayName', '4->1 (Constant-p)');


% Set axes and labels
xlabel('Entropy, s - s1 (J/(kg·K))');
ylabel('Temperature, T (K)');
title('T-s Diagram: 1->2->3->e->4->1');
legend('show');
grid on;
% hold off;
