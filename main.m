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

acc_data = data(:, 2:4);
acc_data = acc_data.*9.81;
acc_data(:, 3) = acc_data(:, 3) - ones(n, 1).*9.81;

gyro_data = data(:, 5:7);

%% MEKF instance initialization
profile on;

model = TranslationRotationModel(Ts);
mekf = MEKF(model);

%% Simulation loop
system_states = zeros(n, mekf.state_size()); % used for graphing data

for i = 1:n
    mekf = mekf.predict_step(get_u(i, acc_data, gyro_data));
    mekf = mekf.update_step(get_z(i, baro_data));
    system_states(i, :) = mekf.x';
end

profile viewer;
writematrix(system_states, "./data/static_out.csv");

%% Data handling functions
function z = get_z(i, baro_data)
z = baro_data(i);
end

function u = get_u(i, acc_data, gyro_data)
u = [acc_data(i, :), gyro_data(i, :)];
end