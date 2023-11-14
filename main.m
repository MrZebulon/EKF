% MEKF code - Samuel Wahba (ICARUS)

clc;
close all;
clear;
init;

debug_mode = false;

%% Data import

Fs = 100;
F_baro = 100;
F_mag = 20;
F_opti = 200;

Ts = 1/Fs;

data = readtable("LA3.csv");
data_input_LA3;

%% MEKF instance initialization
if debug_mode
    profile on;
end

calibration = LA3_Calibration();
model = OptiTrackModelEKF(Ts);
model.set_calibration_data(calibration)

kf = EKFEngine(model);

%% Simulation loop
system_states = zeros(n, kf.state_size()); % used for graphing data

for i = 1:n
    kf = kf.predict_step([acc_data(i, :), gyro_data(i, :)]);

    if mod(i, Fs/F_baro)
          kf = kf.update_step_baro([baro_data(i)]);
    end

    if mod(i, Fs/F_mag)
          kf = kf.update_step_mag([magneto_data(i, :)]);
    end

    if mod(i, Fs/F_opti)
          kf = kf.update_step_optitrack([optitrack_data(i, :)]);
    end
    system_states(i, :) = kf.x';
end

if debug_mode
    profile viewer;
end
writematrix(system_states, "./data/records/out_LA3.csv");