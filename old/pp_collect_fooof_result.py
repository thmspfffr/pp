from fooof import FOOOFGroup
import numpy as np
import scipy.io
import os
import matplotlib.pyplot as plt 

v = 2
all_subj  = [4,5,6,7,8,9,10,11,12,13,15,16,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34]

for isubj in all_subj:
    for iblock in range(1,3):

        try:    
            fm = FOOOFGroup()
            fm.load('/home/tpfeffer/pp/proc/src/pp_hh_fooof_result_lo_s%d_b%d_v%d.json' % (isubj, iblock, v))
            res_lo = fm.get_results()
            fm.load('/home/tpfeffer/pp/proc/src/pp_hh_fooof_result_me_s%d_b%d_v%d.json' % (isubj, iblock, v))
            res_me = fm.get_results()
            fm.load('/home/tpfeffer/pp/proc/src/pp_hh_fooof_result_hi_s%d_b%d_v%d.json' % (isubj, iblock, v))
            res_hi = fm.get_results()
        except:
            print('Could not find data...')
            continue
        
        # Frequencies: 0 - 40 Hz in steps of 0.25 Hz
        F = fm.freqs
        # g_lo will contain the complete PS fit, i.e., aperiodic + periodic components
        g_lo = np.zeros([len(res_lo),len(F)])
        # aper_lo will contain the aperiodic parameters, i.e., offset and slope
        aper_lo = np.zeros([2,len(res_lo)])
        for isrc in range(0,len(res_lo)):
            print('%d' % isrc)
            aper_lo[:,isrc] = res_lo[isrc].aperiodic_params
            for i in range(0,len(res_lo[isrc].gaussian_params)):
                c = res_lo[isrc].gaussian_params[i][0]
                w = res_lo[isrc].gaussian_params[i][2]
                a = res_lo[isrc].gaussian_params[i][1]
                g_lo[isrc,:] = g_lo[isrc,:] + a * np.exp ((-(F-c)**2)/(2*pow(w,2))) 
            b = res_lo[isrc].aperiodic_params[0]
            e = res_lo[isrc].aperiodic_params[1]
            g_lo[isrc,:] = g_lo[isrc,:] + b - np.log10(pow(F,e)) 

        # g_me: ame as g_lo, but for high pupil diameter segments
        g_me = np.zeros([len(res_me),len(F)])
        # aper_me: same as aper_lo, but for high pupil diameter segments
        aper_me = np.zeros([2,len(res_me)])
        for isrc in range(0,len(res_me)):
            print('%d' % isrc)
            aper_me[:,isrc] = res_me[isrc].aperiodic_params
            for i in range(0,len(res_me[isrc].gaussian_params)):
                c = res_me[isrc].gaussian_params[i][0]
                w = res_me[isrc].gaussian_params[i][2]
                a = res_me[isrc].gaussian_params[i][1]
                g_me[isrc,:] = g_me[isrc,:] + a * np.exp ((-(F-c)**2)/(2*pow(w,2))) 
            b = res_me[isrc].aperiodic_params[0]
            e = res_me[isrc].aperiodic_params[1]
            g_me[isrc,:] = g_me[isrc,:] + b - np.log10(pow(F,e)) 

        # g_hi: ame as g_lo, but for high pupil diameter segments
        g_hi = np.zeros([len(res_hi),len(F)])
        # aper_hi: same as aper_lo, but for high pupil diameter segments
        aper_hi = np.zeros([2,len(res_hi)])
        for isrc in range(0,len(res_hi)):
            print('%d' % isrc)
            aper_hi[:,isrc] = res_hi[isrc].aperiodic_params
            for i in range(0,len(res_hi[isrc].gaussian_params)):
                c = res_hi[isrc].gaussian_params[i][0]
                w = res_hi[isrc].gaussian_params[i][2]
                a = res_hi[isrc].gaussian_params[i][1]
                g_hi[isrc,:] = g_hi[isrc,:] + a * np.exp ((-(F-c)**2)/(2*pow(w,2))) 
            b = res_hi[isrc].aperiodic_params[0]
            e = res_hi[isrc].aperiodic_params[1]
            g_hi[isrc,:] = g_hi[isrc,:] + b - np.log10(pow(F,e)) 

        scipy.io.savemat('/home/tpfeffer/pp/proc/src/pp_hh_collected_fooof_s%d_b%d_v%d.mat' % (isubj,iblock,v), {'g_lo': g_lo, 'g_me': g_me, 'g_hi': g_hi, 'aper_lo': aper_lo, 'aper_me': aper_me, 'aper_hi': aper_hi})

all_subj  = range(1,25)

for isubj in all_subj:
    for iblock in range(1,2):

        try:    
            fm = FOOOFGroup()
            fm.load('/home/tpfeffer/pp/proc/src/pp_gla_fooof_result_lo_s%d_b%d_v%d.json' % (isubj, iblock, v))
            res_lo = fm.get_results()
            fm.load('/home/tpfeffer/pp/proc/src/pp_gla_fooof_result_me_s%d_b%d_v%d.json' % (isubj, iblock, v))
            res_me = fm.get_results()
            fm.load('/home/tpfeffer/pp/proc/src/pp_gla_fooof_result_hi_s%d_b%d_v%d.json' % (isubj, iblock, v))
            res_hi = fm.get_results()
        except:
            print('Could not find data...')
            continue
        
        # Frequencies: 0 - 40 Hz in steps of 0.25 Hz
        F = fm.freqs
        # g_lo will contain the complete PS fit, i.e., aperiodic + periodic components
        g_lo = np.zeros([len(res_lo),len(F)])
        # aper_lo will contain the aperiodic parameters, i.e., offset and slope
        aper_lo = np.zeros([2,len(res_lo)])
        for isrc in range(0,len(res_lo)):
            print('%d' % isrc)
            aper_lo[:,isrc] = res_lo[isrc].aperiodic_params
            for i in range(0,len(res_lo[isrc].gaussian_params)):
                c = res_lo[isrc].gaussian_params[i][0]
                w = res_lo[isrc].gaussian_params[i][2]
                a = res_lo[isrc].gaussian_params[i][1]
                g_lo[isrc,:] = g_lo[isrc,:] + a * np.exp ((-(F-c)**2)/(2*pow(w,2))) 
            b = res_lo[isrc].aperiodic_params[0]
            e = res_lo[isrc].aperiodic_params[1]
            g_lo[isrc,:] = g_lo[isrc,:] + b - np.log10(pow(F,e)) 

        # g_me: ame as g_lo, but for high pupil diameter segments
        g_me = np.zeros([len(res_me),len(F)])
        # aper_me: same as aper_lo, but for high pupil diameter segments
        aper_me = np.zeros([2,len(res_me)])
        for isrc in range(0,len(res_me)):
            print('%d' % isrc)
            aper_me[:,isrc] = res_me[isrc].aperiodic_params
            for i in range(0,len(res_me[isrc].gaussian_params)):
                c = res_me[isrc].gaussian_params[i][0]
                w = res_me[isrc].gaussian_params[i][2]
                a = res_me[isrc].gaussian_params[i][1]
                g_me[isrc,:] = g_me[isrc,:] + a * np.exp ((-(F-c)**2)/(2*pow(w,2))) 
            b = res_me[isrc].aperiodic_params[0]
            e = res_me[isrc].aperiodic_params[1]
            g_me[isrc,:] = g_me[isrc,:] + b - np.log10(pow(F,e)) 

        # g_hi: ame as g_lo, but for high pupil diameter segments
        g_hi = np.zeros([len(res_hi),len(F)])
        # aper_hi: same as aper_lo, but for high pupil diameter segments
        aper_hi = np.zeros([2,len(res_hi)])
        for isrc in range(0,len(res_hi)):
            print('%d' % isrc)
            aper_hi[:,isrc] = res_hi[isrc].aperiodic_params
            for i in range(0,len(res_hi[isrc].gaussian_params)):
                c = res_hi[isrc].gaussian_params[i][0]
                w = res_hi[isrc].gaussian_params[i][2]
                a = res_hi[isrc].gaussian_params[i][1]
                g_hi[isrc,:] = g_hi[isrc,:] + a * np.exp ((-(F-c)**2)/(2*pow(w,2))) 
            b = res_hi[isrc].aperiodic_params[0]
            e = res_hi[isrc].aperiodic_params[1]
            g_hi[isrc,:] = g_hi[isrc,:] + b - np.log10(pow(F,e)) 

        scipy.io.savemat('/home/tpfeffer/pp/proc/src/pp_gla_collected_fooof_s%d_b%d_v%d.mat' % (isubj,iblock,v), {'g_lo': g_lo, 'g_me': g_me, 'g_hi': g_hi, 'aper_lo': aper_lo, 'aper_me': aper_me, 'aper_hi': aper_hi})


all_subj  = range(1,42)
all_subj = np.delete(all_subj,[9,11,16,18,21,26,34,37,38,39])

for isubj in all_subj:
    print('%s' % isubj)
    for iblock in range(1,2):

        try:    
            fm = FOOOFGroup()
            fm.load('/home/tpfeffer/pp/proc/src/pp_mue_fooof_result_lo_s%d_b%d_v%d.json' % (isubj, iblock, v))
            res_lo = fm.get_results()
            fm.load('/home/tpfeffer/pp/proc/src/pp_mue_fooof_result_me_s%d_b%d_v%d.json' % (isubj, iblock, v))
            res_me = fm.get_results()
            fm.load('/home/tpfeffer/pp/proc/src/pp_mue_fooof_result_hi_s%d_b%d_v%d.json' % (isubj, iblock, v))
            res_hi = fm.get_results()
        except:
            print('Could not find data...')
            continue
        
        # Frequencies: 0 - 40 Hz in steps of 0.25 Hz
        F = fm.freqs
        # g_lo will contain the complete PS fit, i.e., aperiodic + periodic components
        g_lo = np.zeros([len(res_lo),len(F)])
        # aper_lo will contain the aperiodic parameters, i.e., offset and slope
        aper_lo = np.zeros([2,len(res_lo)])
        for isrc in range(0,len(res_lo)):
            print('%d' % isrc)
            aper_lo[:,isrc] = res_lo[isrc].aperiodic_params
            for i in range(0,len(res_lo[isrc].gaussian_params)):
                c = res_lo[isrc].gaussian_params[i][0]
                w = res_lo[isrc].gaussian_params[i][2]
                a = res_lo[isrc].gaussian_params[i][1]
                g_lo[isrc,:] = g_lo[isrc,:] + a * np.exp ((-(F-c)**2)/(2*pow(w,2))) 
            b = res_lo[isrc].aperiodic_params[0]
            e = res_lo[isrc].aperiodic_params[1]
            g_lo[isrc,:] = g_lo[isrc,:] + b - np.log10(pow(F,e)) 

        # g_me: ame as g_lo, but for high pupil diameter segments
        g_me = np.zeros([len(res_me),len(F)])
        # aper_me: same as aper_lo, but for high pupil diameter segments
        aper_me = np.zeros([2,len(res_me)])
        for isrc in range(0,len(res_me)):
            print('%d' % isrc)
            aper_me[:,isrc] = res_me[isrc].aperiodic_params
            for i in range(0,len(res_me[isrc].gaussian_params)):
                c = res_me[isrc].gaussian_params[i][0]
                w = res_me[isrc].gaussian_params[i][2]
                a = res_me[isrc].gaussian_params[i][1]
                g_me[isrc,:] = g_me[isrc,:] + a * np.exp ((-(F-c)**2)/(2*pow(w,2))) 
            b = res_me[isrc].aperiodic_params[0]
            e = res_me[isrc].aperiodic_params[1]
            g_me[isrc,:] = g_me[isrc,:] + b - np.log10(pow(F,e)) 

        # g_hi: ame as g_lo, but for high pupil diameter segments
        g_hi = np.zeros([len(res_hi),len(F)])
        # aper_hi: same as aper_lo, but for high pupil diameter segments
        aper_hi = np.zeros([2,len(res_hi)])
        for isrc in range(0,len(res_hi)):
            print('%d' % isrc)
            aper_hi[:,isrc] = res_hi[isrc].aperiodic_params
            for i in range(0,len(res_hi[isrc].gaussian_params)):
                c = res_hi[isrc].gaussian_params[i][0]
                w = res_hi[isrc].gaussian_params[i][2]
                a = res_hi[isrc].gaussian_params[i][1]
                g_hi[isrc,:] = g_hi[isrc,:] + a * np.exp ((-(F-c)**2)/(2*pow(w,2))) 
            b = res_hi[isrc].aperiodic_params[0]
            e = res_hi[isrc].aperiodic_params[1]
            g_hi[isrc,:] = g_hi[isrc,:] + b - np.log10(pow(F,e)) 

        scipy.io.savemat('/home/tpfeffer/pp/proc/src/pp_mue_collected_fooof_s%d_b%d_v%d.mat' % (isubj,iblock,v), {'g_lo': g_lo, 'g_me': g_me, 'g_hi': g_hi, 'aper_lo': aper_lo, 'aper_me': aper_me, 'aper_hi': aper_hi})
