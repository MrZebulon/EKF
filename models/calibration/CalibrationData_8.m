classdef CalibrationData_8 < handle
    properties(Constant)
        % units = m/s
        accel_bias = [-1.4316165171268747e-05, -8.43600000000002e-05, -0.003618680000000058];
        accel_noise = 1.717235470188011e-07 * ones(1,3);
        accel_drift = 1e-09;

        % units = rad
        gyro_bias = [1.1856937490705716e-06, -2.845387398998638e-08, 1.589500817925012e-06];
        gyro_noise = 1.066263535014038e-09 * ones(1,3);
        gyro_drift = 1e-09;

        % units = m
        baro_bias = 387.5800978850858;
        baro_noise = 0.014769875002697693;

        % units = uT
        magneto_earth_field = [5.312602240618649, 24.403329945967382, -90.09273224111436];
        magneto_earth_field_drift =  1e-15;

        magneto_bias = [0, 0, 0];
        magneto_noise = 0.8871292885904403 * ones(1,3);
        magneto_drift =  1e-10;
    end
end