% An example file to use the SiSy Class for SiSy1 ZHAW module
%   Created by Sebastian Grünwald, sebiscodes@gmail.com
%   14.12.2023, Winterthur
%   Github: https://github.com/SebisCodes/

%% Simple signal setup by wav
clear; close all; clc;

sisyObj = SiSy; % Init SiSy Object
sisyObj = sisyObj.addWav("dtmf_signal.wav", 1/20); % Second param is the period length in seconds
[t,s,f,N] = sisyObj.getSignal(); % Get the signal and its frequency and sammple amount

disp(sisyObj); % Show values of the sisy object

subplot(1,1,1), plot(t,s); % Plot signal
grid; xlabel('t / s'); ylabel('x(t)');



%% Simple signal setup by own
clear; close all; clc;

t = [0:0.01:1]; % Time from 0 to 1s
s = sin(2*pi*t); % Sin signal over t
fs = 100; % Frequency: 100 samples per second

sisyObj = SiSy; % Init SiSy Object
sisyObj = sisyObj.setSignal(s, fs); % Second param is the period length in seconds
[t2,s2,f,N] = sisyObj.getSignal(); % Get the signal and its frequency and sammple amount

disp(sisyObj); % Show values of the sisy object

subplot(1,1,1), plot(t2,s2); % Plot signal
grid; xlabel('t / s'); ylabel('x(t)');



%% FFT using SiSy
clear; close all; clc;

sisyObj = SiSy; % Init SiSy Object
sisyObj = sisyObj.addWav("dtmf_signal.wav", 1/20); % Second param is the period length in seconds
[t,s,f,N] = sisyObj.getSignal(); % Get the signal and its frequency and sammple amount
[sisyObj, t_p, s_p, fft_f, fft_s] = sisyObj.getFFT(N/21*19); % Get the fft-transformed function

disp(sisyObj); % Show values of the sisy object

subplot(3,1,1), plot(t,s); % Plot signal
subplot(3,1,2), plot(t_p,s_p); % Plot part signal
subplot(3,1,3), plot(fft_f,fft_s); % Plot fft signal

grid; xlabel('t / s'); ylabel('x(t)');



%% Integrals
clear; close all; clc;

t = [0:0.01:1]; % Time from 0 to 1s
s = sin(2*pi*t); % Sin signal over t
fs = 100; % Frequency: 100 samples per second

sisyObj = SiSy; % Init SiSy Object
sisyObj = sisyObj.setSignal(s, fs); % Second param is the period length in seconds
[t2,s2,f,N] = sisyObj.getSignal(); % Get the signal and its frequency and sammple amount
[sisyObj, i_t, i_s, area, absarea] = sisyObj.getIntegral();

disp(sisyObj); % Show values of the sisy object

subplot(2,1,1), plot(t2,s2); % Plot signal
subplot(2,1,2), plot(i_t,i_s); % Plot integral
grid; xlabel('t / s'); ylabel('x(t)');



%% Convolution (Faltung)
clear; close all; clc;

[x,fs] = audioread('aufgabe_4.wav');
x = x'; % Signal x ist ein Zeilenvektor

sisyObj = SiSy; % Init SiSy Object
sisyObj = sisyObj.addWav("aufgabe_4.wav"); % Second param is the period length in seconds
[t,s,f,N] = sisyObj.getSignal(); % Get the signal and its frequency and sammple amount

tau = 1.6e-3; % Timeconstant in seconds
h = (sisyObj.o_Ts/tau) * exp(-t/tau);% Impulse response

[sisyObj, conv_t, conv_s] = sisyObj.getConvolution(h); % Get the convolution using the impulse reaction


disp(sisyObj); % Show values of the sisy object

subplot(2,1,1), plot(t,s,conv_t, conv_s); % Plot signal and convolution

grid; xlabel('t / s'); ylabel('x(t)');