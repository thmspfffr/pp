from fooof import FOOOFGroup
import numpy as np
import scipy.io
import os

# isubj = 4
v=3

for isubj in range (0, 24):

    if os.path.exists('/home/tpfeffer/pp/proc/src/pp_sens_gla_fooof_exp_s%d_proc.txt' % (isubj+1)) == True:
        continue

    os.system('touch /home/tpfeffer/pp/proc/src/pp_sens_gla_fooof_exp_s%d_proc.txt' % (isubj+1))

    dat = scipy.io.loadmat('/home/tpfeffer/pp/proc/src/pp_sens_gla_fooof_s%s_b1_v%d.mat' % (isubj+1,v))

    fm = FOOOFGroup(peak_width_limits=[1, 8], min_peak_height=0.05, max_n_peaks=6)                

    freq_range = [3, 40]

    freqs = np.squeeze(dat['fxx'])
    slp = np.zeros([ np.shape(dat['pxx'])[1] , np.shape(dat['pxx'])[2] ])

    for iseg in range(0,np.shape(dat['pxx'])[2]):
        print('Processing segment %d' % iseg)
        power_spectrum = np.squeeze(dat['pxx'][:,:,iseg])
        spectrum = np.transpose(power_spectrum)
        fm.fit(freqs, spectrum, freq_range)
        tmp=fm.get_results()
        for isens in range(0,np.shape(dat['pxx'])[1]):
            slp[isens][iseg] = tmp[isens][0][1]

    r = np.zeros([np.shape(dat['pxx'])[1],1])
    for iseg in range(0,np.shape(dat['pxx'])[1]):
        r[iseg] = np.corrcoef(slp[iseg,:],dat['pup'],rowvar=True)[0][1]

    scipy.io.savemat('/home/tpfeffer/pp/proc/src/pp_sens_gla_fooof_exp_s%d.mat' % (isubj+1), {'r': r})
