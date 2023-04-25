clc;
close all;
clear;

%% Data import
data = readtable("static.csv");
Ts = 1/160;

n = size(data, 1); % n = #rows
%% MEKF instance initialization
model = TranslationRotationModel(Ts);
mekf = MEKF(model);

%% Simulation loop
system_states = zeros(n, mekf.state_size() + 6); % used for graphing data

for i = 1:n
    mekf = mekf.predict_step(get_u(data, i));
    mekf = mekf.update_step(get_z(data, i));
    system_states(i, :) = [mekf.x' get_u(data, i)];
end

%% Data handling functions
function z = get_z(data, i)
    z = table2array(data(i, ["baro"]));
end

function u = get_u(data, i)
    u = table2array(data(i, ["accel_x", "accel_y", "accel_z", "gyro_x", "gyro_y", "gyro_z"]));
    u = -u;
end