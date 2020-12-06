from fooof import FOOOFGroup
import numpy as np
import scipy.io
import os

v=2

for isubj in range (0, 24):
  for iblock in range(1,3):
        if os.path.exists('/home/tpfeffer/pp/proc/src/pp_gla_src_run_fooof_s%d_b%d_v%d_proc.txt' % (isubj,iblock,v)) == True:
            continue

        os.system('touch /home/tpfeffer/pp/proc/src/pp_gla_src_run_fooof_s%d_b%d_v%d_proc.txt' % (isubj,iblock,v))
        try:
            dat = scipy.io.loadmat('/home/tpfeffer/pp/proc/src/pp_gla_src_powerspectra_s%d_b%d_v%d.mat' % (isubj,iblock,v))
        except:
            print("Error: File not found!")
            continue

        print('Processing S%d B%d ...' % (isubj,iblock))
    
        not_nan_idx  = np.isnan(dat['pxx'][1,1,:])==False

        dat['pxx'] = dat['pxx'][:,:,not_nan_idx]
        dat['pup'] = dat['pup'][:,not_nan_idx]
        dat['pup_df'] = dat['pup_df'][:,not_nan_idx]

        sorted = np.argsort(dat['pup'])[0]

        pup_sorted = dat['pup'][0][sorted]
        pup_means = np.zeros([3,1])
        pup_means[0] = np.mean(pup_sorted[0:int(np.floor(sorted.size/3))])
        pup_means[1] = np.mean(pup_sorted[int(np.floor(sorted.size/3))+1:sorted.size-int(np.floor(sorted.size/3))])
        pup_means[2] = np.mean(pup_sorted[sorted.size-int(np.floor(sorted.size/3))+1:])
        
        np.save('/home/tpfeffer/pp/proc/src/pp_src_gla_fooof_result_pupilmeans_s%d_b%d_v%d.npy' % (isubj,iblock,v),'pup_means')

        first_half  = np.mean(dat['pxx'][:,:,0:int(np.floor(sorted.size/3))],axis=2)
        second_half = np.mean(dat['pxx'][:,:,int(np.floor(sorted.size/3))+1:sorted.size-int(np.floor(sorted.size/3))],axis=2)
        third_half  = np.mean(dat['pxx'][:,:,sorted.size-int(np.floor(sorted.size/3))+1:],axis=2)

        fm = FOOOFGroup(peak_width_limits=[1, 8], min_peak_height=0.05, max_n_peaks=6)                
        fm._maxfev = 30000

        freq_range = [3, 40]

        freqs = np.squeeze(dat['fxx'])
        aperiodic = np.zeros([2,np.shape(dat['pxx'])[1],2])

        fm.fit(freqs, np.transpose(first_half), freq_range)
        fm.save('/home/tpfeffer/pp/proc/src/pp_gla_fooof_result_lo_s%d_b%d_v%d' % (isubj,iblock,v),save_results=True, save_settings=False,save_data=True);

        fm.fit(freqs, np.transpose(second_half), freq_range)
        fm.save('/home/tpfeffer/pp/proc/src/pp_gla_fooof_result_me_s%d_b%d_v%d' % (isubj,iblock,v),save_results=True, save_settings=False,save_data=True);
  
        fm.fit(freqs, np.transpose(third_half), freq_range)
        fm.save('/home/tpfeffer/pp/proc/src/pp_gla_fooof_result_hi_s%d_b%d_v%d' % (isubj,iblock,v),save_results=True, save_settings=False,save_data=True);

