% lorentzspec.m
% 
% This is the joint posterior on omega and alpha from section 6.2 (page 88)
% in Bretthorst - single frequency with Lorentzian decay
%
% James Clark, james.clark@ligo.org

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

function [samples,acc] = lorentzspec2(data,fs,nsamples)
%% MCMC Version

% Target distribution
targetpdf = @(x) logposterior(x, data, fs);

% Proposal distribution
proppdf = @(x,y) proposal(x,y);

% Proposal random sampler
proprnd = @(x) proprand(x);

start=proprand(0.25);

% Generate MH samples
[samples,acc] = mhsample(start,nsamples,'logpdf',targetpdf,...
    'proprnd', proprnd, 'proppdf', proppdf);

% figure
% hist3(samples,[500,500])
% set(gcf,'renderer','opengl');
% set(get(gca,'child'),'FaceColor','interp','CDataMode','auto');

% figure
% hist(samples(:,1),500)
% xlabel('w')

figure('Position', [100, 100, 1000, 300]);
subplot(1,3,1)
xi=0:0.0005:0.5;
f=ksdensity(samples,xi,'support','positive');
% f=ksdensity(samples,xi);
plot(xi,f)
subplot(1,3,2)
hist(samples,500)
subplot(1,3,3)
plot(samples)


%% Subroutines

    % log posterior
    function logprob = logposterior(x, data, fs)
        
        % Sample vals
        w=x(1);
        
        % enforce priors
        if w<0
            logprob=-1e10;
        else
    %         alpha=x(2);

            Ndata=length(data);

            % Mean-squared average of data
            data = data-mean(data);
            time = 1:length(data)/fs;
            d2bar=1/Ndata * sum(data.^2);

            % orthonormal model functions
    %         exppart=exp(-alpha*time);

            H1 = cos(w*time);%.*exppart;
            H2 = sin(w*time);%.*exppart;

            % eigenvalues
            c = Ndata/2 + sin(Ndata*w) ./ (2*sin(w));
            s = Ndata/2 - sin(Ndata*w) ./ (2*sin(w));

    %         % eigenvalues
    %         if alpha==0
    %             c = Ndata/2;
    %             s = c;
    %         else
    %             c = 0.5*( (1-exp(-2*Ndata*alpha)) ./ (exp(2*alpha)-1));
    %             s = c;
    %         end

            % Project model functions onto data
            R2 = sum(data.*H1,2).^2;
            I2 = sum(data.*H2,2).^2;

            logprob = ((2-Ndata)/2) * ...
                log( 1 - (R2/c + I2/s) / (Ndata*d2bar) );
        end
        
    end

    % proposal distribution
    function propprob = proposal(x,y)
        
        w_mean=y(1);
        w_sigma=.01;

        propprob=normpdf(x, w_mean, w_sigma);
        
%         w_mean=y(1);
%         w_sigma=.1;
%         
%         alpha_mean=y(2);
%         alpha_sigma=0.1;
%         
%         mu = [w_mean,alpha_mean];
%         sigma = [w_sigma, 0; 0, alpha_sigma];
%         
%         propprob=mvnpdf(x, mu, sigma);

%         propprob=unifpdf(x,0,0.5);
        
    end

    % RNG for proposal
    function x = proprand(y)
        
%         x = 0.5*rand;
        
        w_mean=y(1);
        w_sigma=.01;
%         
        x=normrnd(w_mean,w_sigma);
        
%         w_mean=y(1);
%         w_sigma=.1;
%         
%         alpha_mean=y(2);
%         alpha_sigma=.1;
%         
%         mu = [w_mean,alpha_mean];
%         sigma = [w_sigma, 0; 0, alpha_sigma];
%         
%         x = mvnrnd(mu, sigma, 1);

        
    end
        
        
       

%     % add in logs
%     function c = logsumexp(a, b)
%         % Computes log(exp(a) + exp(b)) for matrices a and b while guarding against
%         % overfow if exp(a) and or exp(b) are very large
% 
%         d = max(a, b); 
%         d(d == -Inf) = 0;
% 
%         c = log(exp(a - d) + exp(b - d)) + d;
%     end

    end