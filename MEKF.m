classdef MEKF
    %MEKF Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        % x: state vector
        % P: covariance martix
        %
        % f: state transition model
        % gen_F: generator for df/dx (Jacobian Matrix w/ respect to state vector)
        % gen_G: generator for df/dw (Jacobian matrix w/ respect to noise vector)
        %
        % h: observation model
        % gen_H: generator for dh/dx (Jacobian Matrix w/ respect to state vector)
        %
        % noise_gen: generator for noise vector and process noise matrix
        % R: measurement uncertainty matrix
        x;
        P;
        f;
        gen_F;
        gen_G;
        h;
        gen_H;
        gen_noise;
        R;
    end
    
    methods
        function obj = MKEF(x_init, P_init, f, gen_F, gen_G, h, gen_H, gen_noise, R)
            obj.x = x_init;
            obj.P = P_init;
            obj.f = f;
            obj.gen_F = gen_F;
            obj.gen_G = gen_G;
            obj.h = h;
            obj.gen_H = gen_H;
            obj.gen_noise = gen_noise;
            obj.R = R;
        end

        function obj = predict_step(obj, u, Ts)
            F = obj.gen_F(obj.x, u);
            A = eye(size(F)) + F * Ts;
            x_new = obj.f(obj.x, u) * Ts + obj.x; % x = delta_x + x_ref

            [Qs, w] = obj.gen_noise(Ts);
            G = obj.gen_G(obj.x, u, w);
            P_new = A*obj.P*(A') + G*Qs*(G');

            obj.x = x_new;
            obj.P = P_new;
        end
        
        function obj = update_step(obj, z)
            H = obj.gen_H();
            nx = size(obj.x,1);
            inov = z - obj.h(obj.x);
            S = H*obj.P*(H') + obj.R;
            K = obj.P*(H')*inv(S);

            obj.x = obj.x + K*inov;
            obj.P = (eye(nx)-K*H)*obj.P;
        end
    end
end

