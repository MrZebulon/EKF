clc;
close all;
clear;

addpath models\;
%% Data import
data = readtable("empty_data.csv");
Ts = 1/160;

get_u = @(i) table2array(data(i, ["accel_x", "accel_y", "accel_z", "gyro_x", "gyro_y", "gyro_z"]));
get_z = @(i) table2array(data(i, ["baro"]));
n = size(data, 1); % n = #rows
%% MEKF instance initialization
model = TranslationRotationModel(Ts);
mekf = MEKF(model);

%% Simulation loop
system_states = zeros(n, mekf.state_size()); % used for graphing data

for i = 1:n
    mekf = mekf.predict_step(get_u(i));
    mekf = mekf.update_step(get_z(i));
    system_states(i, :) = (mekf.x);
end