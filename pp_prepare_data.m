%% pp_hh_prepare_all_data

clear
restoredefaultpath

SUBJLIST = [4 5 6 7 8 9 10 11 12 13 15 16 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34];

addpath ~/Documents/MATLAB/fieldtrip-20160919/
addpath ~/pconn/matlab/
load(sprintf('~/pp/proc/pp_atlas_BNA.mat'))

outdir = '~/pp/proc/src/';
ord    = pconn_randomization;


%% 200 ms for resting state, 1100 ms for task-counting
% -------------------------
for isubj = SUBJLIST
  
  % identify placebo condition (ord==1)
  im = find(ord(isubj,:)==1);
  
  for iblock = 1:2

    fprintf('Processing subj%d block%d ...\n',isubj,iblock);
    
    try
      % load cleaned meg data
%       load(sprintf('~/pp/data/ham/pupmod_rest_sens_cleandat_s%d_m%d_b%d_v%d.mat',isubj,im,iblock,1))
      load(sprintf('~/pupmod/proc/sens/pupmod_rest_sens_cleandat_s%d_m%d_b%d_v%d.mat',isubj,im,iblock,1))
      load(sprintf('~/pp/proc/pup/pp_pupil_diameter_cleaned_s%d_m%d_b%d.mat',isubj,im,iblock))
    catch me
      continue
    end
    
    load(sprintf('~/pconn/proc/preproc/pconn_preproc_data_s%d_m%d_b%d_v2.mat',isubj,im,iblock))
    label = data.label;
    
    pupil = pupil(1:end-200,:);
    save(sprintf('~/pp/data/ham/pp_rest_s%d_b%d_v%d.mat',isubj,iblock,1),'dat','pupil','label','art')

  end
end


error('!!!')
%%

for isubj = SUBJLIST
  
  % identify placebo condition (ord==1)
  im = find(ord(isubj,:)==1);
  
  for iblock = 1:2

    fprintf('Processing subj%d block%d ...\n',isubj,iblock);
    
    try
      % load cleaned meg data
      load(sprintf('~/pupmod/proc/sens/pupmod_task_sens_cleandat_s%d_m%d_b%d_v%d.mat',isubj,im,iblock,1))
      load(sprintf('~/pp/proc/pup/pp_pupil_diameter_cleaned_counting_s%d_m%d_b%d.mat',isubj,im,iblock))
    catch me
      continue
    end
    
    load(sprintf('~/pconn_cnt/proc/preproc/pconn_cnt_preproc_data_count_s%d_m%d_b%d_v1.mat',isubj,im,iblock))
    label = data.label;
    
%     pupil = pupil(1:end-lag(isubj,iblock),:);
    
    save(sprintf('~/pp/data/ham/pp_task_s%d_b%d_v%d.mat',isubj,iblock,1),'dat','pupil','label')
    
  end
end


%%
% clear siz
% for isubj = SUBJLIST
%     im = find(ord(isubj,:)==1);
% 
%   for iblock = 1 : 2
%     try
% load(sprintf('~/pp/proc/pup/pp_pupil_diameter_cleaned_s%d_m%d_b%d.mat',isubj,im,iblock))
% siz(isubj,1,iblock) = size(pupil,1);
% 
% load(sprintf('~/pp/proc/pup/pp_pupil_diameter_cleaned_counting_s%d_m%d_b%d.mat',isubj,im,iblock))
% siz(isubj,2,iblock) = size(pupil,1);
% 
%     catch me
%     end
% 
% 
%   end
% end
% 
% % siz = nanmean(siz(SUBJLIST,:,:),3);
% siz = siz(SUBJLIST,:,:);
% siz (siz==0) = nan;
for isubj = SUBJLIST
  
  % identify placebo condition (ord==1)
  im = find(ord(isubj,:)==1);
  
  for iblock = 1:2

  load(sprintf('~/pconn_cnt/proc/preproc/pconn_cnt_trig_raw_s%d_m%d_b%d_v1.mat',isubj,im,iblock))
  [p,l]=findpeaks(trig.trial{1});
  lag(isubj,iblock)=1100+1000*(600-(mean(l(end)+find(trig.trial{1}(l(end):end)>0.5)-1)/1200-l(1)/1200))
%   
% sum(p==100)
% (l(find(p==50,1,'last'))-l(find(p==100,1,'first')))/1200
% size(l)
% (l(end)-l(1))/1200
%   l(1)/1200
%   tmp=600-cumsum(diff(l(2:end))/1200); tmp(end)
  end
  
end