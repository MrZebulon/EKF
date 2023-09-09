% MEKF code - Samuel Wahba (ICARUS)

clc;
close all;
clear;
init;

debug_mode = false;

%% Data import

Fs = 100;
F_mag = 20;

Ts = 1/Fs;

data = readtable("raw_in_8.csv");
data = table2array(data);

n = size(data, 1);

baro_data = data(:, 1);
baro_data = 44307.69.*(1-(baro_data.*(1/1013.25)).^0.190284);

acc_data = data(:, 2:4);

gyro_data = data(:, 5:7);
magneto_data = data(:, 8:10);

%% MEKF instance initialization
if debug_mode
    profile on;
end

model = TranslationRotationModelEKF(Ts);
kf = EKFEngine(model);

%% Simulation loop
system_states = zeros(n, kf.state_size()); % used for graphing data

for i = 1:n
    kf = kf.predict_step([acc_data(i, :), gyro_data(i, :)]);

    kf = kf.update_step_baro([baro_data(i)]);
    if mod(i, Fs/F_mag)
          kf = kf.update_step_mag([magneto_data(i, :)]);
    end

    system_states(i, :) = kf.x';
end

if debug_mode
    profile viewer;
end
writematrix(system_states, "./data/out_8.csv");