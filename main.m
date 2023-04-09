clc;
close all;
clear all;
%% MEKF instance initialization
instance = MEKF(zeros(10, 1), diag(1e-9 * ones(1, 10), state_transition, F_generator, G_generator, observation, H_generator, noise_generator, uncertainty));

%% Model definitions
function x_dot = state_transition(x, u, w)
end

function F = F_generator(x, u)
end

function G = G_generator(x, u, w)
end

function z_hat  = observation(x)
end

function H = H_generator()
end

function [Qs, w] = noise_generator(Ts)
end

function R = uncertainty()
end