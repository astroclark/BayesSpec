% bayespec.m
%
% Bayesian spectral estimation following Bretthorst
%
% James clark james.clark@ligo.org

clear
close all

%% General
fs=1024;

%% Signal
% fc=200;
% bw=0.05;
% signal = gauspuls(time,fc,bw);

time = -5:1/fs:5;
% time = 0:512;%fs;
fc = 200;
B1 = 5;
alpha = -1/0.25;
signal = B1*cos(2*pi*fc*time).*exp(time*alpha) + ...
    B1*sin(2*2*pi*fc*time+1).*exp(time*alpha);

signal(time<0) = 0;

%% Get PMNS signal
% filename='waveforms/waveforms.hdf5';
% 
% signal=h5read(filename,['/' 'apr_135135']);
% 
% signal=signal/max(signal);
% 
% time=0:1/fs:length(signal)/fs - 1/fs;


%% Noise

noise = randn(1,length(time));

%% Data

data=signal+noise;

figure
subplot(121)
plot(time,data)
hold on
plot(time,signal,'r')
ylim([-7,7]);

Fmin=1500;
data_f = fft(data);
Cw = abs(data_f).^2;
Cw = 1/length(data) * Cw;
f = [ 0:length(data)/2 , -length(data)/2+1:-1 ]'/ (length(data)/fs);
f = f(f>=0 & f<=fs/2);
Cw = Cw(f>=0 & f<=fs/2);

subplot(122)
semilogy(f,Cw)
% xlim([0,0.5])
% ylim([0,100])

% return

%% Analysis


alpha=-0.1:0.001:0.1;

% Preallocate
[prob, omega] = lorentzspec(data,-5,fs);
[Omega, Alpha] = meshgrid(omega,alpha);
prob_spec = zeros(length(omega),length(alpha));

for a=1:length(alpha)
    [prob, omega] = lorentzspec(data,alpha(a),fs);
    prob_spec(:,a) = prob;
end

prob_freq = sum(prob_spec,2);
prob_alpha = sum(prob_spec);

% prob_freq = lorentzspec(data,0,fs);

% Create spectral estimator
h = spectrum.periodogram;
Hpsd = psd(h, data, 'Fs', fs);
F = Hpsd.Frequencies;
P = Hpsd.Data;

%% Plot results

figure('Position', [100, 100, 800, 400]);
subplot(1,2,1)
plot(alpha,prob_alpha)
subplot(1,2,2)
semilogy(f,prob_freq)


% figure
% surf(Omega/(2*pi),Alpha,log10(prob_spec'),'EdgeColor','none')

