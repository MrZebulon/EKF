classdef EKFEngine < BaseKFEngine
    methods
        function obj = predict_step(obj, u)
            F = obj.model.get_F_matrix(obj.x, u);
            A = eye(size(F)) + F; % adding identity matrix to ensure that process noise isn't mistakenly removed
            [Qs, w] = obj.model.generate_noise();
            Q = obj.model.get_Q_matrix(obj.x, u, w);
            obj.P = A*obj.P*(A') + Q + Qs;

            obj.x = obj.model.compute_x_new(obj.x, u);
            obj.P = 0.5 *(obj.P + (obj.P)');
            % P has to be symmetric. We could use P + P' instead
            % (as it will always be symmetric). Since Pij = Pji
            % correspond to a covariance, we divide by two (so as to not
            % "count twice" the covariance)
        end

        function obj = update_step(obj, z, h, H, R)
            inov = (z-h)';

            S = H*obj.P*(H')+R;
            K = obj.P*(H')*inv(S);

            obj.x = obj.x + K*inov;
            obj.P = (eye(obj.state_size())-K*H)*obj.P;
        end
        
        function obj = update_step_baro(obj, z)
            [h, H, R] = obj.model.baro_obs_model(obj.x);
            obj.update_step(z, h, H, R);
        end

        function obj = update_step_mag(obj, z)
            [h, H, R] = obj.model.mag_obs_model(obj.x);
            obj.update_step(z, h, H, R);
        end
        
        function Nx = state_size(obj)
            Nx = obj.model.Nx;
        end
    end
end
