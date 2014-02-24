% bayespec.m
%
% Bayesian spectral estimation following Bretthorst
%
% James clark james.clark@ligo.org

%% Signal 

fs=8192;
inj.freq=1500;
inj.tau=0.1;
inj.t0=0.25;
inj.amp=.1;

inj.time = 0:1/fs:1 - 1/fs;

inj.name='ringdown';

signal = bayesignals(inj);

% function gendata(B1,alpha_inj)
% 
% %% General
% fs=1;
% 
% %% Signal
% % fc=200;
% % bw=0.05;
% % signal = gauspuls(time,fc,bw);
% 
% time = 1:512;
% fc=0.3;
% % B1 = 10000;
% % alpha_inj = 0.03;
% signal = B1*cos(fc*time).*exp(-time*alpha_inj);
% 
% signal(time<256) = 0;

%% Noise

noise = randn(1,length(inj.time));

% nbar = 1/length(noise) * sum(noise.^2)

%% Data

data=signal+noise;

%% Diagnostic Plots 

figure
subplot(121)
plot(inj.time,data)
hold on
plot(inj.time,signal,'r')
% ylim([-7,7]);

% Create spectral estimator
h = spectrum.periodogram;
Hpsd = psd(h, data, 'Fs', fs, 'nfft', length(data));
F = Hpsd.Frequencies;
P = Hpsd.Data;
subplot(122)
plot(F,P)

