data = table2array(data);

n = size(data, 1);

baro_data = data(:, 10);
%baro_data = 44307.69.*(1-(baro_data.*(1/101325)).^0.190284);

acc_data = data(:, 1:3);

gyro_data = data(:, 4:6);
magneto_data = data(:, 7:9);