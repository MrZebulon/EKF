classdef CalibrationData_11_2_fv < handle
    properties(Constant)
        % units = m/s
        accel_bias = [ 4.08000000e-05, -3.83833333e-04, -1.19450000e-03];
        accel_noise = [2.91324295852940e-08 3.18839440550918e-08 3.00457034526850e-08];
        accel_drift = 1e-09;

        % units = rad
        gyro_bias = [ 1.12933333e-06, -1.01733333e-06, 1.59133333e-06];
        gyro_noise = [3.21835802e-10 1.45008582e-09 1.65409894e-10];
        gyro_drift = 1e-09;

        % units = m
        baro_bias = 435.734225664339;
        baro_noise = 0.0056192747617165705;

        % units = uT
        magneto_earth_field = [-19.22854167  10.38547917   9.1331875];
        magneto_earth_field_drift = 1e-15;
        
        magneto_bias = [0, 0, 0];
        magneto_noise = [0.23051571 0.18192492 0.3947644];
        magneto_drift =  1e-10;
    end
end