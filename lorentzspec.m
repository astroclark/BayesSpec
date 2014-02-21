% lorentzspec.m
% 
% This is the joint posterior on omega and alpha from section 6.2 (page 88)
% in Bretthorst - single frequency with Lorentzian decay
%
% James Clark, james.clark@ligo.org

function [prob, omega] = lorentzspec(data,alpha,fs)

Ndata=length(data);

%% Mean-squared average of data
% data = data - mean(data);

d2bar=mean(data.^2);

%% gij matrix

if alpha==0
    csum = 0.5;
else
    csum = 0.5*( (1-exp(-2*Ndata*alpha )) / (exp(2*alpha)-1));
end
    

csum=0.5;

%% FFT / Periodogram

data_f = fft(data);
Cw = abs(data_f).^2;
Cw = 1/length(data) * Cw;
f = [ 0:length(data)/2 , -length(data)/2+1:-1 ]'/ (length(data)/fs);
f = f(f>=0 & f<=fs/2);
Cw = Cw(f>=0 & f<=fs/2);
omega = 2*pi*f;

%% Posterior calculation
% See eqn: 6.7

prob = (1 - Cw ./ (Ndata*csum*d2bar)).^((2-Ndata)/2);

prob = prob / trapz(prob);



end

