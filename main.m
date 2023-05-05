% MEKF code - Samuel Wahba (ICARUS)

clc;
close all;
clear;

%% Data import
data = readtable("static.csv");
Ts = 1/100;

n = size(data, 1); % n = #rows
%% MEKF instance initialization
%profile on;

model = TranslationRotationModel(Ts);
mekf = MEKF(model);

%% Simulation loop
system_states = zeros(n, mekf.state_size() + 3); % used for graphing data

for i = 1:n
    mekf = mekf.predict_step(get_u(data, i));
    mekf = mekf.update_step(get_z(data, i));
    system_states(i, :) = [mekf.x' mekf.P(4, 4) mekf.P(5, 5) mekf.P(6, 6)];
end

%profile viewer;

%% Data handling functions
function z = get_z(data, i)
    z = table2array(data(i, ["baro"]));
end

function u = get_u(data, i)
    u = table2array(data(i, ["acc_x", "acc_y", "acc_z", "gyro_x", "gyro_y", "gyro_z"]));
end