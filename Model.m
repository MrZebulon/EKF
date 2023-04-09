classdef Model
    methods
         function x_new = get_next_state(obj, x, u, dt)
              x_new = 0;
            end
        
            function F = get_F_matrix(obj, x, u, dt)
                F = 0;
            end
            
            function G = get_G_matrix(obj, x, u, w, dt)    
                G = 0;
            end
            
            function z_hat  = get_measurement_estimate(obj, x)
                z_hat = 0;
            end
            
            function H = get_H_matrix(obj)
                H = 0;
            end
            
            function [Qs, w] = generate_noise(obj, dt)
                Qs = 0;
                w = 0;
            end
            
            function R = get_R_matrix(obj)
                R = 0;
            end
    end

end