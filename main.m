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

data = readtable("half_in_11_2.csv");
data_input_2;

%% MEKF instance initialization
if debug_mode
    profile on;
end

calibration = CalibrationData_11_2();
model = TranslationRotationModelEKF(Ts);
model.set_calibration_data(calibration)

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
writematrix(system_states, "./data/records/out_half_11_2_b.csv");