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

subplot(1,1,1), plot(t,s);grid; % Plot signal
grid;



%% Simple signal setup by own
clear; close all; clc;

t = [0:0.01:1]; % Time from 0 to 1s
s = sin(2*pi*t); % Sin signal over t
fs = 100; % Frequency: 100 samples per second

sisyObj = SiSy; % Init SiSy Object
sisyObj = sisyObj.setSignal(s, fs); % Second param is the period length in seconds
[t2,s2,f,N] = sisyObj.getSignal(); % Get the signal and its frequency and sammple amount

disp(sisyObj); % Show values of the sisy object

subplot(1,1,1), plot(t2,s2);grid; % Plot signal
grid;



%% FFT using SiSy
clear; close all; clc;

sisyObj = SiSy; % Init SiSy Object
sisyObj = sisyObj.addWav("dtmf_signal.wav", 1/20); % Second param is the period length in seconds
[t,s,f,N] = sisyObj.getSignal(); % Get the signal and its frequency and sammple amount
[sisyObj, t_p, s_p, fft_f, fft_s] = sisyObj.getFFT(N/21*19); % Get the fft-transformed function

disp(sisyObj); % Show values of the sisy object

subplot(3,1,1), plot(t,s);grid; % Plot signal
subplot(3,1,2), plot(t_p,s_p);grid; % Plot part signal
subplot(3,1,3), plot(fft_f,fft_s);grid; % Plot fft signal

grid;



%% Integrals
clear; close all; clc;

t = [0:0.01:1]; % Time from 0 to 1s
s = sin(2*pi*t); % Sin signal over t
fs = 100; % Frequency: 100 samples per second

sisyObj = SiSy; % Init SiSy Object
sisyObj = sisyObj.setSignal(s, fs, 50, 25); % Second param is the period length in seconds
%sisyObj = sisyObj.setSignal(s, fs); % Integrate the whole function
[t2,s2,f,N] = sisyObj.getSignal(); % Get the signal and its frequency and sammple amount
[sisyObj, i_t, i_s, area, absarea] = sisyObj.getIntegral();

disp(sisyObj); % Show values of the sisy object

subplot(2,1,1), plot(t2,s2);grid; % Plot signal
subplot(2,1,2), plot(i_t,i_s);grid; % Plot integral
grid;



%% Convolution (Faltung)
clear; close all; clc;

[x,fs] = audioread('aufgabe_4.wav');
x = x'; % Column to row

sisyObj = SiSy; % Init SiSy Object
sisyObj = sisyObj.addWav("aufgabe_4.wav"); % Second param is the period length in seconds
[t,s,f,N] = sisyObj.getSignal(); % Get the signal and its frequency and sammple amount

tau = 1.6e-3; % Timeconstant in seconds
h = (sisyObj.o_Ts/tau) * exp(-t/tau);% Impulse response

[sisyObj, conv_t, conv_s] = sisyObj.getConvolution(h); % Get the convolution using the impulse reaction


disp(sisyObj); % Show values of the sisy object

subplot(2,1,1), plot(t,s,conv_t, conv_s);grid; % Plot signal and convolution




%% Get Ak & Bk coefficients
clear; close all; clc;

% Parameter
f0 = 1;     % Frequency
N = 2000;   % Amount of samples
fs = 1000;  % Sampling Frequency

t = [0:N-1]*(2/2000);   % Time 

s=sign(sin(2*pi*f0*(t-0.4)))+sin(t); %Rechtecksignal

N_coeff = 32; % Amount of Ak, Bk and Mk coefficients

sisyObj = SiSy; % Init SiSy Object
sisyObj = sisyObj.setSignal(s, fs); % Second param is the period length in seconds
%sisyObj = sisyObj.setSignal(s, fs, 1000, 500); % ATTENTION! MIGHT BE BUGGY! Second param is the period length in seconds
[t2,s2,fs,N] = sisyObj.getSignal(); % Get the signal and its frequency and sammple amount

[sisyObj, ct, cs, aS, aSA, aSB, Ak, Bk, Mk] = sisyObj.getFourierCoefficients(N_coeff);

disp(sisyObj); % Show values of the sisy object

subplot(4,2,1), plot(t2,s2); % Plot signal
grid; xlabel('t / s'); ylabel('s(t)');
subplot(4,2,2), plot(ct,cs); % Plot signal
grid; xlabel('t / s'); ylabel('part of s(t)');
subplot(4,2,3), plot(ct,aS); % Plot approximated signal
grid; xlabel('t / s'); ylabel('Approximated part of s(t)');
subplot(4,2,5), plot(ct,aSA); % Plot A-part signal
grid; xlabel('t / s'); ylabel('A-part');
subplot(4,2,7), plot(ct,aSB); % Plot B-part signal
grid; xlabel('t / s'); ylabel('B-part');
subplot(4,2,4), stem([1:sisyObj.coeff_N],Ak); hold on; % Plot B-part signal
grid; xlabel('t / s'); ylabel('Ak');
subplot(4,2,6), stem([1:sisyObj.coeff_N],Bk); % Plot B-part signal
grid; xlabel('t / s'); ylabel('Bk');
subplot(4,2,8), stem([1:sisyObj.coeff_N],Mk); % Plot B-part signal
grid; xlabel('t / s'); ylabel('Mk');