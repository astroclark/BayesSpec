% lorentzspec2.m
% 
%
% James Clark, jclark.astro@gmail.com


function [Data] = lorentzspec2(data,fs,nsamples)

%% Data struct

Data.data_td=data; % change to FD here if you want
Data.time=0:1/fs:length(data)/fs - 1/fs;
Data.sigma=1; % Noise estimation!
Data.nsamples=nsamples;
Data.fs=fs;

%% Initialise prior
% TODO: add Prior to Data

% --- Ranges
Prior.freqmin=1000;
Prior.freqmax=2000;
Prior.taumin=0;
Prior.taumax=0.5;
Prior.ampmin=0;
Prior.ampmax=5;
Prior.t0min=0;
Prior.t0max=Data.time(end);

% --- Prior PDFs
Prior.freq = @(x) unifpdf(x,Prior.freqmin,Prior.freqmax);
Prior.tau  = @(x) unifpdf(x,Prior.taumin,Prior.taumax);
Prior.amp  = @(x) unifpdf(x,Prior.ampmin,Prior.ampmax);
Prior.t0   = @(x) unifpdf(x,Prior.t0min,Prior.t0max);

Prior.freqrnd = @(~) unifrnd(Prior.freqmin,Prior.freqmax);
Prior.taurnd  = @(~) unifrnd(Prior.taumin,Prior.taumax);
Prior.amprnd  = @(~) unifrnd(Prior.ampmin,Prior.ampmax);
Prior.t0rnd   = @(~) unifrnd(Prior.t0min,Prior.t0max);

% --- Running history for covariance calculation
Data.samples.freq  = zeros(1,nsamples);
Data.samples.tau = zeros(1,nsamples);
Data.samples.amp   = zeros(1,nsamples);
Data.samples.t0    = zeros(1,nsamples);

%% Generate Samples

initial_params=[Prior.freqrnd(),Prior.taurnd(),Prior.amprnd(),Prior.t0rnd()];
Data.counter=1;


% Target distribution
targetpdf_h = @(params) logtargetpdf(params);

mh=1;
if mh
    % Proposal distribution
    proppdf = @(x,y) proposal_pdf(x,y);

    % Proposal random sampler
    proprnd = @(x) proposal_smplr(x);

    % Generate MH samples
    [samples,acc] = mhsample(initial_params,nsamples,'logpdf',targetpdf_h,...
        'proprnd', proprnd, 'logproppdf', proppdf, 'nchain', 1);
else
    
    slicesample(initial_params,nsamples,'logpdf',targetpdf_h);

end

%% Insert samples into Data structure & finalise

Data.samples = params2struct(samples);
Data.samples.acc = acc;

    
    %% Probability Distribution Functions
    
    % --- Target distribution: likelihood x prior
    function logp = logtargetpdf(params)
       
        % Convert params to Params structure
        Params = params2struct(params);
        
        % Add to the sample history
        % FIXME URGENT: ONLY RETAIN ACCEPTED SAMPLES! THERE WAS SOMETHING I
        % SAW...
        Data.samples.freq(Data.counter) = Params.freq;
        Data.samples.tau(Data.counter)  = Params.tau;
        Data.samples.amp(Data.counter)  = Params.amp;
        Data.samples.t0(Data.counter)   = Params.t0;
        
        % Increment the iteration counter
        Data.counter=Data.counter+1;
        
        % Posterior = Likelihood x Prior
        logp = loglikelihood(Params) + logprior(Params);
        
    end

    % --- Likelihood function: Gaussian
    function logl = loglikelihood(Params)
        
        
        % generate waveform
        % TODO: don't do this if only params.amp changes
        waveform = ringdown_td(Params,Data);
        
        % compute likelihood
        residuals = Data.data_td - waveform;
        
        % The actual likelihood
        logl = -0.5*sum( (residuals ./ Data.sigma).^2 );
        
    end
        
    % --- Prior pdf: assume logically independent params
    function logpi = logprior(params)
        
        logpi = ...
            log(Prior.freq(params.freq) * ...
            Prior.tau(params.tau) * ...
            Prior.amp(params.amp) * ...
            Prior.t0(params.t0));
        
        logpi(isnan(logpi)) = -1e100;
        
    end

    % --- Proposal distribution
    function q = proposal_pdf(new_smpl,old_smpl)
        % Slice sample from target distribution
%         q = logtargetpdf(new_smpl);
         
%          Data.counter
%          shape(Data.samples.freq)

        if Data.counter<100
             q = logprior(params2struct(new_smpl));
        else
            % TODO: move the covariance stuff somewhere more 'central'
            if Data.counter<=Data.nsamples
                sample_matrix = [...
                    Data.samples.freq(1:Data.counter); ...
                    Data.samples.tau(1:Data.counter); ...
                    Data.samples.amp(1:Data.counter); ...
                    Data.samples.t0(1:Data.counter)...
                    ];
            else
                sample_matrix = [...
                    Data.samples.freq; ...
                    Data.samples.tau; ...
                    Data.samples.amp; ...
                    Data.samples.t0...
                    ];
            end
            % Multivariate normal based on current sample and covariance
%             S = cholcov(cov(sample_matrix'));
            q = log(mvnpdf(new_smpl,old_smpl,cov(sample_matrix')));
        end

    end

    % --- Proposal sampler
    function smpl = proposal_smplr(old_smpl)
        % Slice sample from target distribution
        %         smpl = slicesample(old_smpl,1,'logpdf',targetpdf_h);
        
        if Data.counter<100
            smpl=[Prior.freqrnd(),...
                Prior.taurnd(),...
                Prior.amprnd(),...
                Prior.t0rnd()];
        else
            % Multivariate normal based on current sample and covariance
            if Data.counter<=Data.nsamples
                sample_matrix = [...
                    Data.samples.freq(1:Data.counter); ...
                    Data.samples.tau(1:Data.counter); ...
                    Data.samples.amp(1:Data.counter); ...
                    Data.samples.t0(1:Data.counter)...
                    ];
            else
                sample_matrix = [...
                    Data.samples.freq; ...
                    Data.samples.tau; ...
                    Data.samples.amp; ...
                    Data.samples.t0...
                    ];
            end
            
            % Multivariate normal based on current sample and covariance
%             S = cholcov(cov(sample_matrix'));
            smpl = mvnrnd(old_smpl,cov(sample_matrix'));
        end
        
    end

    %% Waveform model
    
    % --- Time domain ringdown
    function waveform = ringdown_td(Params,Data)
        
        time=Data.time;
        
        nonzeroidx=Data.time>Params.t0;
        waveform = zeros(1,length(time));
        waveform(nonzeroidx) = Params.amp * ...
            cos(2*pi*Params.freq*(Data.time(nonzeroidx) - Params.t0)) .* ...
            exp(-(Data.time(nonzeroidx)-Params.t0)/Params.tau);
        
    end

    %% Helper Functions
    
    % Convert params array (for use with samplers) to structure (for use
    % with user-defined functions).  Not sure this is actually necessary.
    function Struct = params2struct(params)
        
        if size(params,1)==1 
            Struct.freq = params(1);
            Struct.tau  = params(2);
            Struct.amp  = params(3);
            Struct.t0   = params(4);
        else
            Struct.freq = params(:,1);
            Struct.tau  = params(:,2);
            Struct.amp  = params(:,3);
            Struct.t0   = params(:,4);
        end
    end


end