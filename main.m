% MEKF code - Samuel Wahba (ICARUS)

clc;
close all;
clear;
init;

%% Data import
Ts = 1/100;

data = readtable("static_raw.csv");
data = table2array(data);

n = size(data, 1);

baro_data = data(:, 1);
baro_data = 44307.69.*(1-(baro_data.*(1/1013.25)).^0.190284);
baro_data = movmean(baro_data, 10);

acc_data = data(:, 2:4);
acc_data = acc_data.*9.81;
acc_data(:, 3) = acc_data(:, 3) - ones(n, 1).*9.81;
acc_data = movmean(acc_data, 10);

gyro_data = data(:, 5:7);
gyro_data = movmean(gyro_data, 10);
writematrix([baro_data, acc_data, gyro_data], "../Data/static/cleansed.csv");
writematrix([baro_data, acc_data, gyro_data], "../Data/static/raw_shai.csv");


%% MEKF instance initialization
profile on;

model = TranslationRotationModel(Ts);
ekf = EKFEngine(model);

%% Simulation loop
system_states = zeros(n, ekf.state_size()); % used for graphing data

for i = 1:n
    ekf = ekf.predict_step(get_u(i, acc_data, gyro_data));
    ekf = ekf.update_step(get_z(i, baro_data));
    system_states(i, :) = ekf.x';
end

%profile viewer;
writematrix(system_states, "./data/static_out.csv");
writematrix(system_states, "../Data/static/run.csv");


%% Data handling functions
function z = get_z(i, baro_data)
z = baro_data(i);
end

function u = get_u(i, acc_data, gyro_data)
u = [acc_data(i, :), gyro_data(i, :)];
end