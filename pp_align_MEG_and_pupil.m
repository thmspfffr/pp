%%  
% --------
% cut out and align MEG and pupil recordings
% takes data from ~/pupmod/matlab/pupmod_align_data.m
% is used at input for:
% (1) pp_src_pupil_power_correlations.m
% --------
% Note that Task data is processed separately in 
% pp_cnt_align_MEG_and_pupil.m
% --------

clear
restoredefaultpath

SUBJLIST  = [4 5 6 7 8 9 10 11 12 13 15 16 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34];
v = 1;

addpath ~/Documents/MATLAB/fieldtrip-20160919/
ft_defaults

outdir = '~/pp/proc/sens/';
addpath ~/pconn/matlab/
ord   = pconn_randomization;
%%
for isubj = SUBJLIST
  
  im = find(ord(isubj,:)==1);
  
  for iblock = 1:2
    
    fn = sprintf('pp_align_MEG_and_pupil_s%d_b%d_v%d',isubj,iblock,v);
    if tp_parallel(fn,outdir,1,0)
      continue
    end
    
    fprintf('Processing subj%d block%d ...\n',isubj,iblock);

    try
      % load meg data from pupmod_align_data.m (same as for powcorr)
      load(sprintf('~/pupmod/proc/pupmod_sens_cleandat_s%d_m%d_b%d_v%d.mat',isubj,m,iblock,1))
      % load pupil data
      load(sprintf('~/pp/proc/pup/pp_pupil_diameter_cleaned_s%d_m%d_b%d.mat',isubj,im,iblock))
    catch me
      continue
    end

    data.trial = data.trial(:,1:data.end_of_recording);
    data.trial = data.trial(:,end:-1:1);
    
    if size(pupil,2)>1
      raw_pupil=pupil(:,4);
    else
      raw_pupil=pupil;
    end
    
    pup = resample(raw_pupil,400,1000);
    pup = pup(end:-1:1);
      
    len = min([size(pup,1) size(data.trial,2)]);
    if len/400 > 600
      len = 400*600;
    end

    pup = pup(1:len);
    pup = pup(end:-1:1);

    data.trial = data.trial(:,1:len);
    data.trial = data.trial(:,end:-1:1);
    data.trial(:,isnan(pup))=nan;
    
    dat = data.trial; clear data
    
    save([outdir fn],'dat','pup');
    
    tp_parallel(fn,outdir,0,0)

  end
end
  
