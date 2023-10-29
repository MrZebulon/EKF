data = table2array(data);

n = size(data, 1);

baro_data = data(:, 1);
baro_data = 44307.69.*(1-(baro_data.*(1/1013.25)).^0.190284);

acc_data = data(:, 2:4);

gyro_data = data(:, 5:7);
magneto_data = data(:, 8:10);