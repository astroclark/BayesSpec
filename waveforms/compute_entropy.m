% compute_entropy.m
%
% Loads PMNS waveforms from hdf5 and computes their entropy
%
% James Clark, james.clark@ligo.org

filename='waveforms.hdf5';
snr_filename='snr.hdf5';

wf_info = h5info(filename);
wf_data_info = wf_info.Datasets;

snr_info = h5info(snr_filename);
snr_data_info = snr_info.Datasets;

% Preallocate
wf_name=cell(1,length(wf_data_info));
wfs=cell(1,length(wf_data_info));
wf_entropy=zeros(1,length(wf_data_info));
rho_optimal=zeros(1,length(wf_data_info));


%% Precomputed vals

rho50=struct;

names={'apr_135135',...
'dd2_135135_lessvisc',...
'dd2_165165',...
'nl3_135135_lessvisc',...
'nl3_1919_lessvisc',...
'sfho_135135_lessvisc',...
'sfho_1616',...
'sfhx_135135_lessvisc',...
'shen_135135_lessvisc',...
'tm1_135135_lessvisc',...
'tma_135135_lessvisc',...
'hybrid_apr_135135_lessvisc_0_0p5',...
'hybrid_apr_135135_lessvisc_0p05_0p5',...
'hybrid_dd2_135135_lessvisc_0_0p5',...
'hybrid_dd2_135135_lessvisc_0p05_0p5'};

D50=[2.79, 3.21,2.91,3.77,3.88,2.54,3.29,2.75,3.54,3.13,3.20,5.05,...
    3.95,7.76,5.96];

Dstar=20/sqrt(5);

%% Compute entropy

for i=1:length(wf_data_info)
    current_wf = wf_data_info(i);
    wf_name{i} = current_wf.Name;
    wf=h5read(filename,['/' wf_name{i}]);
    
    rho_optimal(i) = h5read(snr_filename,['/' wf_name{i}]);
    
    % comute shannon entropy of waveform, normalised to unity peak
    % amplitude
    wf_normed = wf/max(wf);
    wfs{i} = wf_normed;
    wf_entropy(i) = wentropy(wf_normed,'logenergy');
    
end

%% Match everything up

d50=zeros(1,length(wf_name));
for i=1:length(wf_name)
    current_wf=wf_name{i};
    idx = find(strcmp(names, wf_name{i}));
    
    d50(i) = D50(idx);
    
end

rho50 = Dstar*rho_optimal./d50;



