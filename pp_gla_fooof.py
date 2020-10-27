from fooof import FOOOFGroup
import numpy as np
import scipy.io
import os

# isubj = 4
v=3

files = [f for f in os.listdir('/home/tpfeffer/pp/proc/src/') if (f.startswith('pp_sens_gla_fooof')) and (f.endswith('b1_v%d.mat' % v))]

for isubj in range (0, len(files)):

    if os.path.exists('/home/tpfeffer/pp/proc/src/pp_sens_gla_fooof_exp_s%d_proc.txt' % isubj) == True:
        continue

    os.system('touch /home/tpfeffer/pp/proc/src/pp_sens_gla_fooof_exp_s%d_proc.txt' % isubj)

    fn = files[isubj] 

    dat = scipy.io.loadmat('/home/tpfeffer/pp/proc/src/%s' % fn)

    fm = FOOOFGroup(peak_width_limits=[1, 8], min_peak_height=0.05, max_n_peaks=6)                

    freq_range = [3, 40]

    freqs = np.squeeze(dat['fxx'])
    slp = np.zeros([ np.shape(dat['pxx'])[1] , np.shape(dat['pxx'])[2] ])

    for iseg in range(0,np.shape(dat['pxx'])[2]):
        print('Processing segment %d' % iseg)
        power_spectrum = np.squeeze(dat['pxx'][:,:,iseg])
        spectrum = np.transpose(power_spectrum)
        fm.fit(np.squeeze(freqs), spectrum, freq_range)
        tmp=fm.get_results()
        for isens in range(0,np.shape(dat['pxx'])[1]):
            slp[isens][iseg] = tmp[isens][0][1]

    r = np.zeros([np.shape(dat['pxx'])[1],1])
    for iseg in range(0,np.shape(dat['pxx'])[1]):
        r[iseg] = np.corrcoef(slp[iseg,:],dat['pup'],rowvar=True)[0][1]

    scipy.io.savemat('/home/tpfeffer/pp/proc/src/pp_sens_gla_fooof_exp_s%d.mat' % isubj, {'r': r})
