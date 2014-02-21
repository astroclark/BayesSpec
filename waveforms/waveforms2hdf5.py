#!/usr/bin/env python
"""
waveforms2hdf5.py loops over the list of waveforms defined in this script and
dumps out an hdf5 file for the plus polarisation.  The idea is to then compute
the Shannon entropy of the waveforms using Matlab's wentropy.m function.
"""

import h5py
import numpy as np
import pmns_utils

wfs='/Users/jclark/hmns_repo/results/penultimate_waveforms.txt'
waveform_list=np.loadtxt(wfs,dtype=str)
#waveform_list=['shen_135135_lessvisc','apr_135135']

h5_file=h5py.File('waveforms.hdf5','w')
h5_snr_file=h5py.File('snr.hdf5','w')

for waveform in waveform_list:
    # Generate waveform instance
    wf=pmns_utils.Waveform(waveform)

    # Compute the time series & SNR
    wf.make_wf_timeseries()
    wf.compute_characteristics()

    # Zoom in on signal
    peak_idx=np.argmax(wf.hplus.data.data)
    wf_start_idx=np.argwhere(abs(wf.hplus.data.data)>0)[0]
    wf_end_idx=np.argwhere(abs(wf.hplus.data.data)>0)[-1]
    wf_reduced = wf.hplus.data.data[wf_start_idx:wf_end_idx]

    h5_file[waveform]     = wf_reduced
    h5_snr_file[waveform] = wf.snr_plus
    #h5_file[waveform]=wf_reduced
    #h5_file[waveform+'_snr']=wf.snr_plus

h5_file.close()
h5_snr_file.close()
