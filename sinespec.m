% sinespec.m
% 
% This is the  posterior on omega and from Bretthorst
%
% This version computes the projections of the basis functions onto the
% data explicitly for demonstration purposes.  We can, of course, use the
% FFT.
%
% James Clark, james.clark@ligo.org

function [logprob, w] = sinespec(data,fs)

Ndata=length(data);

%% Mean-squared average of data

data = data-mean(data);
time = 1:length(data)/fs;
d2bar=1/Ndata * sum(data.^2);

%% FFT / Periodogram

% NFFT = 4096;
% data_f = fft(data,NFFT);
% R = real(data_f).^2;
% I = imag(data_f).^2;
% f = fs/2*linspace(0,1,NFFT/2+1);
% R = R(1:NFFT/2+1);
% I = I(1:NFFT/2+1);

deltaw=0.00005;
w=0.2:deltaw:0.4;
logprob=zeros(1,length(w));
logmass=-1e10;

% orthonormal model functions
H1 = cos(w'*time);
H2 = sin(w'*time);

% eigenvalues
c = Ndata/2 + sin(Ndata*w) ./ (2*sin(w));
s = Ndata/2 - sin(Ndata*w) ./ (2*sin(w));

% Projetions
Data = repmat(data,[length(w),1]);
R2 = sum(Data.*H1,2).^2;
I2 = sum(Data.*H2,2).^2;

R2=R2(:);
I2=I2(:);
c=c(:);
s=s(:);

logprob = ((2-Ndata)/2) * ...
    log( 1 - (R2./c + I2./s) / (Ndata*d2bar) );
logmass=logsumexp(sum(logprob), log(0.00005));
logprob=logprob-logmass;

function c = logsumexp(a, b)
% Computes log(exp(a) + exp(b)) for matrices a and b while guarding against
% overfow if exp(a) and or exp(b) are very large

d = max(a, b); 
d(d == -Inf) = 0;

c = log(exp(a - d) + exp(b - d)) + d;

end

end