from fooof import FOOOFGroup
import numpy as np
import scipy.io
import os
import time

v=2

for isubj in range(1,42):
    for iblock in range(1,2):

        if os.path.exists('/home/tpfeffer/pp/proc/src/pp_mue_fooof_s%d_b%d_v%d_proc.txt' % (isubj,iblock,v)) == True:
            continue

        os.system('touch /home/tpfeffer/pp/proc/src/pp_mue_fooof_s%d_b%d_v%d_proc.txt' % (isubj,iblock,v))
        try:
            dat = scipy.io.loadmat('/home/tpfeffer/pp/proc/src/pp_mue_src_powerspectra_s%d_b%d_v%d.mat' % (isubj,iblock,v))
        except:
            print("Error: File not found!")
            continue

        print('Processing S%d B%d ...' % (isubj,iblock))
        
        
        not_nan_idx  = np.isnan(dat['pxx'][1,1,:])==False

        dat['pxx'] = dat['pxx'][:,:,not_nan_idx]
        pup = dat['pup'][:,not_nan_idx]; 
        pup_df = dat['pup_df'][:,not_nan_idx]

        if v > 33: 

          sorted = np.argsort(pup)[0]
          pup_sorted = pup[0][sorted]
          pup_means = np.zeros([3,1])
          pup_means[0] = np.mean(pup_sorted[0:int(np.floor(sorted.size/3))])
          pup_means[1] = np.mean(pup_sorted[int(np.floor(sorted.size/3))+1:sorted.size-int(np.floor(sorted.size/3))])
          pup_means[2] = np.mean(pup_sorted[sorted.size-int(np.floor(sorted.size/3))+1:])

          np.save('/home/tpfeffer/pp/proc/src/pp_src_mue_fooof_result_pupilmeans_s%d_b%d_v%d.npy' % (isubj,iblock,v),'pup_means')

          first_half  = np.mean(dat['pxx'][:,:,0:int(np.floor(sorted.size/3))],axis=2)
          second_half = np.mean(dat['pxx'][:,:,int(np.floor(sorted.size/3))+1:sorted.size-int(np.floor(sorted.size/3))],axis=2)
          third_half  = np.mean(dat['pxx'][:,:,sorted.size-int(np.floor(sorted.size/3))+1:],axis=2)

          fm = FOOOFGroup(peak_width_limits=[1, 8], min_peak_height=0.05, max_n_peaks=6)                
          fm._maxfev = 30000

          freq_range = [3, 40]

          freqs = np.squeeze(dat['fxx'])
          aperiodic = np.zeros([2,np.shape(dat['pxx'])[1],2])

          fm.fit(freqs, np.transpose(first_half), freq_range)
          fm.save('/home/tpfeffer/pp/proc/src/pp_mue_fooof_result_lo_s%d_b%d_v%d' % (isubj,iblock,v),save_results=True, save_settings=False,save_data=True);

          fm.fit(freqs, np.transpose(second_half), freq_range)
          fm.save('/home/tpfeffer/pp/proc/src/pp_mue_fooof_result_me_s%d_b%d_v%d' % (isubj,iblock,v),save_results=True, save_settings=False,save_data=True);

          fm.fit(freqs, np.transpose(third_half), freq_range)
          fm.save('/home/tpfeffer/pp/proc/src/pp_mue_fooof_result_hi_s%d_b%d_v%d' % (isubj,iblock,v),save_results=True, save_settings=False,save_data=True);

        else:

          freqs = np.squeeze(dat['fxx'])
          aper = np.empty([2,dat['pxx'].shape[1],dat['pxx'].shape[2]])  
          g = np.zeros([dat['pxx'].shape[1],75,dat['pxx'].shape[2]])
          gg = np.zeros([dat['pxx'].shape[1],75,dat['pxx'].shape[2]])

          for iseg in range(0,dat['pxx'].shape[2]):
            print('%d / %d' % (iseg,dat['pxx'].shape[2]))
            fm = FOOOFGroup(peak_width_limits=[1, 8], min_peak_height=0.05, max_n_peaks=6)                
            fm._maxfev = 30000
            freq_range = [3, 40]
            fm.fit(freqs, np.transpose(dat['pxx'][:,:,iseg]), freq_range)
            tmp = fm.get_results()

            for isens in range(0,dat['pxx'].shape[1]):
              aper[:,isens,iseg] = tmp[isens].aperiodic_params

            F = fm.freqs
            for isens in range(0,dat['pxx'].shape[1]):
              for i in range(0,len(tmp[isens].gaussian_params)):
                  c = tmp[isens].gaussian_params[i][0]
                  w = tmp[isens].gaussian_params[i][2]
                  a = tmp[isens].gaussian_params[i][1]
                  g[isens,:,iseg] = g[isens,:,iseg] + a * np.exp ((-(F-c)**2)/(2*pow(w,2))) 
              b = tmp[isens].aperiodic_params[0]
              e = tmp[isens].aperiodic_params[1]
              gg[isens,:,iseg] = g[isens,:,iseg] + b - np.log10(pow(F,e))

          scipy.io.savemat('/home/tpfeffer/pp/proc/src/pp_mue_collected_fooof_s%d_b%d_v%d.mat' % (isubj,iblock,v), {'g': g, 'gg': gg,  'aper':  aper, 'pup': pup, 'pup_df': pup_df})








        