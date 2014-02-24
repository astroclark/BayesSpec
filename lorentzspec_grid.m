% lorentzspec.m
% 
% This is the joint posterior on omega and alpha from section 6.2 (page 88)
% in Bretthorst - single frequency with Lorentzian decay
%
% James Clark, james.clark@ligo.org

function [logprob, Omega, Alpha] = lorentzspec(data,fs)

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
% R = R(1:NFFT/2+1)/Ndata;
% I = I(1:NFFT/2+1)/Ndata;

w=0.2:0.0005:0.4;
alpha=-0.1:0.005:0.1;
% alpha=0.01;
logprob=zeros(length(w),length(alpha));
logmass=-1e10;
for i=1:length(w)
    for j=1:length(alpha)
        
        % orthonormal model functions
        H1 = cos(w(i)*time).*exp(-alpha(j)*time);
        H2 = sin(w(i)*time).*exp(-alpha(j)*time);
        
        % eigenvalues
        if alpha(j)==0
            c = Ndata/2;
            s = c;
        else
            c = 0.5*( (1-exp(-2*Ndata*alpha(j))) ./ (exp(2*alpha(j))-1));
            s = c;
        end
        
        R2 = sum(data.*H1).^2;
        I2 = sum(data.*H2).^2;

        logprob(i,j) = ((2-Ndata)/2) * ...
            log( 1 - (R2/c + I2/s) / (Ndata*d2bar) );
        
        logmass = logsumexp(logmass, logprob(i,j));
        
    end
end
logprob = logprob-logmass+log(0.0005)+log(0.005);

[Omega, Alpha] = meshgrid(w,alpha);

function c = logsumexp(a, b)
% Computes log(exp(a) + exp(b)) for matrices a and b while guarding against
% overfow if exp(a) and or exp(b) are very large

d = max(a, b); 
d(d == -Inf) = 0;

c = log(exp(a - d) + exp(b - d)) + d;

end

end

% function [logprob, W, Alpha] = lorentzspec2(data,fs)


%% Grid-based Approach

% Ndata=length(data);
% 
% %% Mean-squared average of data
% 
% data = data-mean(data);
% time = 1:length(data)/fs;
% d2bar=1/Ndata * sum(data.^2);
% 
% %% FFT / Periodogram
% 
% % set up grids
% dw=0.0001;
% w=0.0:dw:0.5;
% dalpha=0.001;
% alpha=-0.1:dalpha:0.1;
% 
% % construct meshed versions
% [W, Alpha] = meshgrid(w,alpha);
% Data = repmat(data,[length(w),1]);
% 
% logprob=zeros(length(w),length(alpha));
% 
% disp('looping over alpha')
% for i=1:length(alpha)
% 
%     % orthonormal model functions
%     exppart=exp(-alpha(i)*time);
%     Exppart=repmat(exppart,[length(w),1]);
%     
%     H1 = cos(w'*time).*Exppart;
%     H2 = sin(w'*time).*Exppart;
% 
%     % eigenvalues
%     if alpha(i)==0
%         c = Ndata/2;
%         s = c;
%     else
%         c = 0.5*( (1-exp(-2*Ndata*alpha(i))) ./ (exp(2*alpha(i))-1));
%         s = c;
%     end
% 
%     % Project model functions onto data
%     Data = repmat(data,[length(w),1]);
%     R2 = sum(Data.*H1,2).^2;
%     I2 = sum(Data.*H2,2).^2;
% 
%     logprob(:,i) = logposterior(R2, I2, c, s, Ndata, d2bar);
% 
% %     logmass = logsumexp(logmass, logprob(:,j));
% 
% end