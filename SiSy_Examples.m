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



%% FFT
clear; close all; clc;

sisyObj = SiSy; % Init SiSy Object
sisyObj = sisyObj.addWav("aufgabe_4.wav"); % Second param is the period length in seconds
[t,s,f,N] = sisyObj.getSignal(); % Get the signal and its frequency and sample amount
%[sisyObj, t_p, s_p, fft_f, fft_s] = sisyObj.getFFT(N/10, N/10*3.5); % Get the fft-transformed function
[sisyObj, t_p, s_p, fft_f, fft_s] = sisyObj.getFFT(); % Get the fft-transformed function

disp(sisyObj); % Show values of the sisy object

subplot(3,1,1), plot(t,s);grid; % Plot signal
xlabel('t'); ylabel('s(t)');
subplot(3,1,2), plot(t_p,s_p);grid; % Plot part signal
xlabel('t'); ylabel('s(t)');
subplot(3,1,3), plot(fft_f,fft_s);grid; % Plot fft signal
xlabel('f'); ylabel('fft');



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
N = 1000;
t = [0:N-1]*(f0/N);   % Time 
fs = f0*N;

%s=sign(sin(2.0*pi*f0*(t-0.4)))+sin(t); %Rechtecksignal
%s = cos(2*pi*f0*t);
s = sin(2*pi*f0*t);

N_coeff = 10; % Amount of Ak, Bk and Mk coefficients

sisyObj = SiSy; % Init SiSy Object
sisyObj = sisyObj.setSignal(s, fs,N); % Second param is the period length in seconds
%sisyObj = sisyObj.setSignal(s, fs, 1000, 500); % ATTENTION! MIGHT BE BUGGY! Second param is the period length in seconds
[t2,s2,fs,N] = sisyObj.getSignal(); % Get the signal and its frequency and sammple amount

[sisyObj, ct, cs, aS, aSA, aSB, Ak, Bk, Mk, ck, Pk] = sisyObj.getFourierCoefficients(N_coeff);

disp(sisyObj); % Show values of the sisy object

subplot(5,2,1), plot(t2,s2); % Plot signal
grid; xlabel('t / s'); ylabel('s(t)');
subplot(5,2,3), plot(ct,cs); % Plot signal
grid; xlabel('t / s'); ylabel('part of s(t)');
subplot(5,2,5), plot(ct,aS); % Plot approximated signal
grid; xlabel('t / s'); ylabel('Approximated part of s(t)');
subplot(5,2,7), plot(ct,aSA, ct, aSB); % Plot A-part signal
grid; xlabel('t / s'); ylabel('A-part / B-part'); legend('Cos (A)','Sin (B)');
subplot(5,2,2), stem([1:sisyObj.coeff_N]-1,Ak); % Plot Ak values
grid; xlabel('Index * f0'); ylabel('Ak');
subplot(5,2,4), stem([1:sisyObj.coeff_N]-1,Bk); % Plot Bk values
grid; xlabel('Index * f0'); ylabel('Bk');
subplot(5,2,6), stem([1:sisyObj.coeff_N]-1,Mk); % Plot Mk values
grid; xlabel('Index * f0'); ylabel('Mk');
subplot(5,2,8), stem([1:sisyObj.coeff_N]-1,Pk); % Plot Pk values
grid; xlabel('Index * f0'); ylabel('Pk');
% ATTENTION!!! ck may have wrong values!!!
subplot(5,2,9), stem([1:sisyObj.coeff_N]-1,abs(ck)); % Plot abs(ck) values
grid; xlabel('Index * f0'); ylabel('|ck|');
subplot(5,2,10), stem([1:sisyObj.coeff_N]-1,angle(ck)); % Plot arg(ck) values
grid; xlabel('Index * f0'); ylabel('arg(ck)');



%% Bode plot (IN WORK)

tau = [-15:0.1:15];
fg = 1/(2*pi*tau);



% Plot Bode magnitude
subplot(2,1,1);
semilogx(w, 20*log10(mag)); % Convert to dB
grid on;
title('Bode Magnitude Plot');
xlabel('Frequency (rad/s)');
ylabel('Magnitude (dB)');

% Plot Bode phase
subplot(2,1,2);
semilogx(w, phase);
grid on;
title('Bode Phase Plot');
xlabel('Frequency (rad/s)');
ylabel('Phase (degrees)');