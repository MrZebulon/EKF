classdef LA3_Calibration < handle
    properties(Constant)
        % units = m/s
        accel_bias = [ 4.08000000e-05, -3.83833333e-04, -1.19450000e-03];
        accel_noise = 8.109744803823443e-08 * ones(1,3);
        accel_drift = 1e-09;

        % units = rad
        gyro_bias = [ 1.12933333e-06, -1.01733333e-06, 1.59133333e-06];
        gyro_noise = 1.0439025981994016e-09 * ones(1,3);
        gyro_drift = 1e-09;

        % units = m
        baro_bias = 435.734225664339;
        baro_noise = 0.0056192747617165705;

        % units = uT
        magneto_earth_field = [-19.22854167  10.38547917   9.1331875];
        magneto_earth_field_drift =  1e-15;
        
        magneto_bias = [0, 0, 0];
        magneto_noise = 0.7998584779996147 * ones(1,3);
        magneto_drift =  1e-10;
    end
end