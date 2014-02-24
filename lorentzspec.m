% lorentzspec.m
% 
% This is the joint posterior on omega and alpha from section 6.2 (page 88)
% in Bretthorst - single frequency with Lorentzian decay
%
% This version uses Metropolis-Hastings or Slice sampling.  Beware the
% hard-coded proposal widths for the MH-implementation.  Slice sampling
% handles this automatically.
%
% Currently, this generates samples for a sinusoid with Lorentzian decay
% which starts at time zero and it seems limited by numerical overflows if
% i move to more physicaly intereresting units.
%
% TODO: 
%
% 0) Include an arbitrary start time
%
% 1) Compute the Bayes factor for the model.  This could come from
% the harmonic mean of the posterior at the sample points if i had access
% to that and it was normalised.  An easier alternative might be a
% Laplacian approximation to the posterior? 
% E.g., https://www.stat.washington.edu/~raftery/Research/PDF/kass1995.pdf
%
% 2) Generalise to arbitrary models following Bretthorst's eigenvalue stuff
%
% James Clark, jclark.astro@gmail.com


function [samples,acc] = lorentzspec(data,fs,nsamples,sampling_algorithm)

%% Generate Samples

start=[0.25,0.01, 200];

% Target distribution
targetpdf = @(x) logposterior(x, data, fs);

% Sampling algorithm to use
if nargin==3
    sampling_algorithm='slice'; % slice-sampling is default
end

if strcmp(sampling_algorithm,'mh')
    disp('Performining Metropolis-Hastings sampling')

    % initialise proposal struct (idea is to have something we can update)
    prop_struct.w_sigma=0.001;
    prop_struct.alpha_sigma=0.001;
    prop_struct.t0_sigma=0.01;

    % Proposal distribution
    proppdf = @(x,y) proposal_pdf(x,y);

    % Proposal random sampler
    proprnd = @(x) proposal_sampler(x);

    % Generate MH samples
    [samples,acc] = mhsample(start,nsamples,'logpdf',targetpdf,...
        'proprnd', proprnd, 'proppdf', proppdf, 'symmetric', 1);
else
    disp('Performining slice sampling')
    
    % Generate slice samples
    samples = slicesample(start,nsamples,'logpdf',targetpdf,'thin',5,'burnin',1000);
    acc=nan;
end

%% Compute Bayes Factor
disp('computing Bayes factor from harmonic mean')

% --- Harmonic Mean
logprobs=zeros(1,nsamples);
for i=1:nsamples
    logprobs(i)=targetpdf(samples(i,:));
end


    %% log posterior - Bretthorst
    function logprob = logposterior(x, data, fs)
        
        % Sample vals
        w=x(1);
        alpha=x(2);
        t0=x(3);
        
        % enforce priors
        if w<0 || w>0.5 || alpha<0 || alpha>.1 || t0<0 || t0>length(data)
            logprob=-1e10;
        else
            
            % Derived quantities
            Ndata =length(data);
            data  = data-mean(data);
            time  = 1:length(data)/fs;
            d2bar = 1/Ndata * sum(data.^2);

            % orthonormal model functions
            exppart = exp(-alpha*time);
            H1 = cos(w*(time-t0)).*exppart;
            H2 = sin(w*(time-t0)).*exppart;
            
            H1(time<t0)=0;
            H2(time<t0)=0;

            % eigenvalues
            if alpha==0
                c = Ndata/2;
                s = c;
            else
                c = 0.5*( (1-exp(-2*Ndata*alpha)) ./ (exp(2*alpha)-1));
                s = c;
            end

            % Project data onto model basis
            R2 = sum(data.*H1,2).^2;
            I2 = sum(data.*H2,2).^2;

            logprob = ((2-Ndata)/2) * ...
                log( 1 - (R2/c + I2/s) / (Ndata*d2bar) );
        end
        
    end

    %% proposal distribution
    
    function p = proposal_pdf(x,y)
        
        w_mean=y(1);
        w_sigma=prop_struct.w_sigma;
        
        alpha_mean=y(2);
        alpha_sigma=prop_struct.alpha_sigma;
        
        t0_mean=y(3);
        t0_sigma=prop_struct.t0_sigma;
        
        mu = [w_mean,alpha_mean,t0_mean];
        sigma = [w_sigma, 0, 0; 0, alpha_sigma, 0; 0, 0, t0_sigma];
        
        p=mvnpdf(x, mu, sigma);
 
    end

    %% RNG for proposal
    
    function new_sample = proposal_sampler(old_sample)
        
        w_mean=old_sample(1);
        w_sigma=prop_struct.w_sigma;
        
        alpha_mean=old_sample(2);
        alpha_sigma=prop_struct.alpha_sigma;
        
        t0_mean=old_sample(3);
        t0_sigma=prop_struct.t0_sigma;
        
        mu = [w_mean,alpha_mean,t0_mean];
        sigma = [w_sigma, 0, 0; 0, alpha_sigma, 0; 0, 0, t0_sigma];
        
        new_sample = mvnrnd(mu, sigma, 1);
        
    end
        
end