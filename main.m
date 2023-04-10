clc;
close all;
clear;
%% Data import
data = readtable("empty_data.csv");
Ts = 1/160;

get_u = @(i) table2array(data(i, ["accel_x", "accel_y", "accel_z"]));
get_z = @(i) table2array(data(i, ["baro"]));
n = size(data, 1); % n = #rows
%% MEKF instance initialization
model = SimplifiedModel(Ts);
mekf = MEKF(zeros(10, 1), diag(1e-9 * ones(1, 10)), model);

%% Simulation loop
system_states = zeros(n, mekf.state_size());

for i = 1:n
    u = get_u(i);
    z = get_z(i);
    mekf = mekf.predict_step(u);
    mekf = mekf.update_step(z);
    system_states(i, :) = (mekf.x);
end
