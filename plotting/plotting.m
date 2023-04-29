close all
clc
clear

%% Import
data = load("./data/static_out.csv");
%data = system_states;

N = size(data, 1);
dt = 1/100;

pos = data(:, 1:3);
vel = data(:, 4:6);
ori = data(:, 7:10);

progressive = true;

%% VELOCITY PLOT
figure
t = linspace(0, N * dt, N);

subplot(3,1,1);
vx = plot(t(1), vel(1, 1));
axis([min(t), max(t), min(vel(:, 1)), max(vel(:, 1))]);
ylabel("v_x (m/s)")
grid on

subplot(3,1,2);
vy = plot(t(1), vel(1, 2));
axis([min(t), max(t), min(vel(:, 2)), max(vel(:, 2))]);
ylabel("v_y (m/s)")
grid on

subplot(3,1,3);
vz = plot(t(1), vel(1, 3));
axis([min(t), max(t), min(vel(:, 3)), max(vel(:, 3))]);
ylabel("v_z (m/s)")
grid on

%% ORIENTATION PLOT
figure
t = linspace(0, N * dt, N);
[yaw, pitch, roll] = quat2angle(ori);

subplot(3,1,1);
ay = plot(t(1), yaw(1));
axis([min(t), max(t), min(yaw(:)), max(yaw(:))]);
ylabel("Yaw (rad)")
grid on

subplot(3,1,2);
ap = plot(t(1), pitch(1));
axis([min(t), max(t), min(pitch(:)), max(pitch(:))]);
ylabel("Pitch (rad)")

grid on

subplot(3,1,3);
ar = plot(t(1), roll(1));
axis([min(t), max(t), min(roll(:)), max(roll(:))]);
ylabel("Roll (rad)")
grid on

%% POSITION PLOT
figure

point = scatter3(pos(1, 1),pos(1, 2),pos(1,3), 'red', 'diamond', 'LineWidth', 1);
hold on

trace = plot3(pos(1, 1),pos(1, 2),pos(1,3), "b--");
hold on

axis([-max(1, max(pos(:, 1))), max(1, max(pos(:, 1))), -max(1, max(pos(:, 2))), max(1, max(pos(:, 2))), -max(1, max(pos(:, 3))), max(1, max(pos(:, 3)))]);
xlabel("x (m)")
ylabel("y (m)")
zlabel("z (m)")

grid on

%% RENDER

if progressive
    for i =  1:N
        set(vx, 'XData', t(1:i), 'YData', vel(1:i, 1))
        set(vy, 'XData', t(1:i), 'YData', vel(1:i, 2))
        set(vz, 'XData', t(1:i), 'YData', vel(1:i, 3))
        drawnow limitrate
    
    
        set(ay, 'XData', t(1:i), 'YData', yaw(1:i))
        set(ap, 'XData', t(1:i), 'YData', pitch(1:i))
        set(ar, 'XData', t(1:i), 'YData', roll(1:i))
        drawnow limitrate
    
        set(trace,'Xdata',pos(1:i, 1),'Ydata',pos(1:i, 2),'Zdata',pos(1:i, 3));
        hold on
        set(point, 'Xdata',pos(i, 1),'Ydata',pos(i, 2),'Zdata',pos(i, 3))
        drawnow limitrate
    end
else
    set(vx, 'XData', t(1:N), 'YData', vel(1:N, 1))
    set(vy, 'XData', t(1:N), 'YData', vel(1:N, 2))
    set(vz, 'XData', t(1:N), 'YData', vel(1:N, 3))
    drawnow limitrate


    set(ay, 'XData', t(1:N), 'YData', yaw(1:N))
    set(ap, 'XData', t(1:N), 'YData', pNtch(1:N))
    set(ar, 'XData', t(1:N), 'YData', roll(1:N))
    drawnow limitrate

    set(trace,'Xdata',pos(1:N, 1),'Ydata',pos(1:N, 2),'Zdata',pos(1:N, 3));
    hold on
    set(point, 'Xdata',pos(N, 1),'Ydata',pos(N, 2),'Zdata',pos(N, 3))
    drawnow limitrate
end