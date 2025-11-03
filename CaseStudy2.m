%% Case Study 2, RC Circut Code and RL Circut Code
% *Name:
%% Instructions
%|publish('Lab_CS2.m','pdf')| in the _Command Window_
clear;
close all;  % uncomment this line if you do not want all figure windows to close when running this code
%% Part 1 - RC Equations
%%%
% Implement Equations (8) and (10) in MATLAB to simulate circuit A (Figure
% 1) with R = 1 kΩ and C = 1 µF. Simulate charging the capacitor using a
% step input: set V_C = 0 V at t = 0 and V_in = 1 V for t > 0. Choose a
% suitable h to model the charging process accurately. Plot V_in and V_C
% vs. time to show the charging of the capacitor. You should observe a
% charging curve similar to the Figure 2 in CS2.
%%%
% Set up simulation parameters:
% TODO: *************************************************************
%Setting up R and C
R = 1e3;
C = 1e-6;
%Choosing an intial value for h, I used that as a step value for 0.005
%seconds
h = 0.00001;
t = 0:h:0.005;
V_in = ones(size(t));
V_in(1) = 0;
V_C = zeros(size(t));
% *******************************************************************
% (please finish the function 'simRCvoltages' in the end and use it here.)
% TODO: *************************************************************
for i = 1:length(V_C);
   V_C(i+1) = (1 - (h/(R*C)))*(V_C(i)) + (h/(R*C))*V_in(i);
end
% *******************************************************************
%%%
% Plot the figure:
% TODO: *************************************************************
plot(V_C);
hold on;
plot(V_in);
ylabel("Voltage");
xlabel("Time (s)");
legend('V_C','V_in');
title("Voltage coming In Shown with the Voltage of the Capacitor");
% *******************************************************************
%% Part 2: Comparison between sampling intervals
% Now, run several versions of your simulation for various temporal sampling intervals h. As h gets
% larger or smaller, how does your simulation's prediction change? Why is this happening? Does
% the charging behavior of a "real" capacitor change as a function of your choice of h?
% Plot the predicted V_out using
% 1) an "accurate" choice of h and
% 2) an "inaccurate" choice of h and
% 3) the theoretical charging curve
% V_C(t) = 1 − exp(−t/RC) [Volts].
% Discuss what happens for the "inaccurate" choice of h.
% Note: Be careful to compare the three curves using correct time axes.
%%%
% Simulate with large h:
% TODO: *************************************************************
%h = 0.1;
h = 0.001;
%h = 0.0001;
%h = 0.00001;
t = 0:h:0.005;
V_C2 = zeros(size(t));
V_in2 = ones(size(t));
for i = 1:length(V_C2);
   V_C2(i+1) = (1 - (h/(R*C)))*(V_C2(i)) + (h/(R*C))*V_in2(i);
end
figure;
plot(V_in2);
hold on;
plot(V_C2);
%plot(V_C);
title("Same Graph with Changing H Values");
legend("V_in","V_C2");
xlabel("Time (s)");
ylabel("Volage");
h = 0.00001;
t = 0:h:0.005;
V_in3 = ones(size(t));
V_theoretical = zeros(size(t));
time = zeros(size(t));
for i = 1:(length(time) - 1);
   t(i+1) = t(i) + h;
end
for i = 1:length(V_theoretical);
       V_theoretical(i) = 1 - exp(-t(i)/(R*C));
end;
figure;
plot(V_theoretical);
hold on;
plot(V_in3);
title("Theoretical Graphs");
ylabel("Voltage");
xlabel("Time");
legend ("V_theoretical", "V_in3");
% *******************************************************************
%% Helper functions "simRCvoltages": [V_C, V_R] = simRCvoltages(V_in,V_C0,R,C,h)
% TODO: *************************************************************
% (Replace this comment with code)
% (Add response here)
%function [V_C, V_R] = simRCvoltages(V_in,V_C0,R,C,h)
%end
%function [V_C, V_R] = simRCvoltages(V_in,V_C,R,C,h)
%    V_C(i+1) = (1-(h/R*C))*(V_C(i)) + (h/(R*C))*V_in(i);
%end
% *******************************************************************
%% *******************************************************************
%Setting Up the RL Circuit
R = 100;
L = 100e-3;
h = 0.000001;
k = 0:h:0.005;
%Setting up V_in, V_C, and ik
V_in = ones(size(k));
V_in(1) = 0;
i = zeros(size(k));
for j = 1:length(i)
   i(j+1) = ((1 - (h*R)/L) * i(j)) + ((h/L) * (V_in(j)));
end;
V_L = zeros(size(k));
for j = 1:length(V_L)
   V_L(j) = V_in(j) - R*(i(j));
end
%% **********************************************************************
%Establishing X
figure;
plot(V_in);
hold on;
plot(V_L);
hold off;


%% ***********************************************************************
%% ***********************************************************************
%% ***********************************************************************
%Running the Function first testing for different R Values
%R = 100;
R_values = [300, 100, -15];
%L = 100e-3;
L_values = [100e-3, 1000e-3, 100e-3];
%C = 1e-6;
C_values = [1e-7, 1e-7, 1e-6];
for n = 1:3
   R = R_values(n);
   L = L_values(n);
   C = C_values(n);
   [k, V_R, V_L, V_C, i] = simulateRLC(R, L, C);
   eval(sprintf('V_R%d = V_R;', n));
end
figure;
plot(k, V_in,'DisplayName', 'V_{in}');
hold on;
for n = 1:3
   eval(sprintf('plot(k, V_R%d, ''DisplayName'', ''V_R%d'');', n, n));
end
xlabel('Time (s)');
ylabel('Voltage (V)');
legend show;
%% ***********************************************************************
soundsc(V_R1,fs);
pause(2);
soundsc(V_R2, fs);
pause(2);
soundsc(V_R3, fs);
%% ***********************************************************************
%Defining a function, because I realized I didn't want to rerun the code
%for each minor change (when experimenting)
function [k, V_R, V_L, V_C, i] = simulateRLC(R, L, C)
%Setting up the Constants
fs = 192e3;
h = 1/fs;
k = 0:h:0.015;
V_in = ones(size(k));
V_in(1) = 0;
A = [[1, h/C]; [-h/L, (1 - (h*R)/L)]];
B = [0, h/L]';
x = zeros(2, length(k));
V_C = zeros(size(k));
i = zeros(size(k));
V_L = zeros(size(k));
for n = 2:length(x)
   x(:, n) = A * x(:, n-1) + B * V_in(n);
   V_C(n) = x(1, n);
   i(n) = x(2, n);
   V_L(n) = L * (i(n) - i(n-1))/h;
end
%Establishing R (Which will be graphed)
V_R = zeros(size(k));
%V_R1 = zeros(size(k));
for m = 2:length(V_R)
   V_R(m) = R * i(m);
   %V_R1(m) = V_in(m) - V_L(m) - V_C(m);
end
end


