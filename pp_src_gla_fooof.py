from fooof import FOOOFGroup
import numpy as np
import scipy.io
import os

# isubj = 3
v=4

for isubj in range (0, 24):

    if os.path.exists('/home/tpfeffer/pp/proc/src/pp_gla_src_fooof_exp_s%d_v%d_proc.txt' % (isubj+1,v)) == True:
        continue

    os.system('touch /home/tpfeffer/pp/proc/src/pp_gla_src_fooof_exp_s%d_v%d_proc.txt' % (isubj+1,v))

    dat = scipy.io.loadmat('/home/tpfeffer/pp/proc/src/pp_gla_src_fooof_s%s_b1_v%d.mat' % (isubj+1,v))

    not_nan_idx  = np.isnan(dat['pxx'][1,1,:])==False
    # not_nan_idx2 = np.isnan(dat['pup'])==False
    # not_nan_idx = (not_nan_idx1) & (not_nan_idx2)

    dat['pxx'] = dat['pxx'][:,:,not_nan_idx]
    dat['pup'] = dat['pup'][:,not_nan_idx]
    dat['pup_df'] = dat['pup_df'][:,not_nan_idx]

    fm = FOOOFGroup(peak_width_limits=[1, 8], min_peak_height=0.05, max_n_peaks=6)                

    freq_range = [3, 40]

    freqs = np.squeeze(dat['fxx'])
    slp = np.zeros([ np.shape(dat['pxx'])[1] , np.shape(dat['pxx'])[2] ])

    for iseg in range(0,np.shape(dat['pxx'])[2]):
        print('Processing segment %d / %d' % (iseg,np.shape(dat['pxx'])[2]))
        power_spectrum = np.squeeze(dat['pxx'][:,:,iseg])
        fm.fit(freqs, np.transpose(power_spectrum), freq_range)
        tmp=fm.get_results()
        for isens in range(0,np.shape(dat['pxx'])[1]):
            slp[isens][iseg] = tmp[isens][0][1]

    scipy.io.savemat('/home/tpfeffer/pp/proc/src/pp_gla_src_fooof_slp_s%d_v%d.mat' % (isubj+1,v), {'slp': slp})

    r = np.zeros([np.shape(dat['pxx'])[1],1])
    for iseg in range(0,np.shape(dat['pxx'])[1]):
        r[iseg] = np.corrcoef(slp[iseg,:],dat['pup'],rowvar=True)[0][1]

    scipy.io.savemat('/home/tpfeffer/pp/proc/src/pp_gla_src_fooof_exp_s%d_v%d.mat' % (isubj+1,v), {'r': r})
