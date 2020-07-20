%% pp_save_osf
% correlate bandpass filtered (or via wavelets) pupil and MEG signals

clear
SUBJLIST  = [4 5 6 7 8 9 10 11 12 13 15 16 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34];
v = 1;
addpath ~/Documents/MATLAB/fieldtrip-20160919/
ft_defaults

outdir = '~/pp/proc/sens/';
addpath ~/pconn/matlab/
ord = pconn_randomization;


for isubj = 8
  
  im = find(ord(isubj,:)==1);
  
  for iblock = 2
    
%     fn = sprintf('pp_sens_pupil_power_correlations_s%d_b%d_v%d',isubj,iblock,v);
%     if tp_parallel(fn,outdir,1,0)
%       continue
%     end
%     
    fprintf('Processing subj%d block%d ...\n',isubj,iblock);
    
    try
      % load meg data
      load(sprintf('~/pp/proc/pp_sens_cleandat_s%d_m%d_b%d_v%d.mat',isubj,im,iblock,v))
      % load pupil data
      load(sprintf('~/pp/proc/pup/pp_pupil_diameter_cleaned_s%d_m%d_b%d.mat',isubj,im,iblock))
      
    catch me
%       r = nan(274,25);
%       save([outdir fn '.mat'],'r')
      continue
    end
       
    pupil = resample(pupil,400,1000);
    
    data.trial = data.trial(:,1:data.end_of_recording);
    data.time  = data.time(:,1:data.end_of_recording);
    data.trial = data.trial(:,end:-1:1);
    data.time  = data.time(:,end:-1:1);
    pupil = pupil(end:-1:1,:);
    
    len = min([size(pupil,1) size(data.trial,2)]);
    if len/400 > 600
      len = 400*600;
    end
    
    data.trial = data.trial(:,1:len);
    data.time = data.time(:,1:len);
    pupil = pupil(1:len,:);
    
    data.trial = data.trial(:,end:-1:1);
    data.time = data.time(:,end:-1:1);
    pupil = pupil(end:-1:1,:);
    
    data=rmfield(data,{'start_of_recording';'end_of_recording';'art';'idx';'trigger'});
    
    save(sprintf('~/pp/proc/pp_osf_s%d_b%d.mat',isubj,iblock),'data','pupil')
    
    
    
  end
end
