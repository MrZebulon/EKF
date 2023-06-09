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
%writematrix([baro_data, acc_data, gyro_data], "../Data/static/cleansed.csv");
%writematrix([baro_data, acc_data, gyro_data], "../Data/static/raw_shai.csv");


%% MEKF instance initialization
profile on;

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
%writematrix([baro_data, acc_data, gyro_data], "../Data/static/cleansed.csv");
%writematrix([baro_data, acc_data, gyro_data], "../Data/static/raw_shai.csv");


%% MEKF instance initialization

ref = EKFEngine(TranslationRotationModel(Ts));
test = EKFEngine(TranslationRotationModelV2(Ts));


%% Simulation loop
drift = zeros(n, 2); % used for graphing data

for i = 1:n
    ref = ref.predict_step(get_u(i, acc_data, gyro_data));
    ref = ref.update_step(get_z(i, baro_data));

    test = test.predict_step(get_u(i, acc_data, gyro_data));
    test = test.update_step(get_z(i, baro_data));
    drift(i) = norm(test.D - ref.D, "fro");
end

writematrix(drift, "../Data/static/drift.csv");


%% Data handling functions
function z = get_z(i, baro_data)
z = baro_data(i);
end

function u = get_u(i, acc_data, gyro_data)
u = [acc_data(i, :), gyro_data(i, :)];
end