%% pp_src_pupil_power_correlations
% correlate bandpass filtered (or via wavelets) pupil and MEG signals

clear
restoredefaultpath

SUBJLIST  = [4 5 6 7 8 9 10 11 12 13 15 16 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34];
% focus on all pupil
v = 1;
% focus on low freq pupil
v = 2
addpath ~/Documents/MATLAB/fieldtrip-20160919/
ft_defaults

outdir = '~/pp/proc/src/';
addpath ~/pconn/matlab/
ord = pconn_randomization;
load(sprintf('~/pp/proc/pp_atlas_BNA.mat'))


for isubj = SUBJLIST
  
  im = find(ord(isubj,:)==1);
  isubj
  
  for iblock = 1:2

    
    fprintf('Processing subj%d block%d ...\n',isubj,iblock);
    
    try
      % load meg data
%       load(sprintf('~/pp/proc/pp_sens_cleandat_s%d_m%d_b%d_v%d.mat',isubj,im,iblock,1))
      % load pupil data
      load(sprintf('~/pp/proc/pup/pp_pupil_diameter_cleaned_s%d_m%d_b%d.mat',isubj,im,iblock))
      
    catch me
      all_pup(:,isubj,iblock) = nan(100000,1);
%       save([outdir fn '.mat'],'src_r')
      continue
    end
    
    k = 2;
    f_sample = 1000;
    fnq = f_sample/2; 
    hil_hi = 0.005;
    hil_lo = 200;
    hil_Wn=[hil_hi/fnq hil_lo/fnq];
    [bhil, ahil] = butter(k, hil_Wn);
    
    tmp = filtfilt(bhil, ahil, pupil(:,4));
    pupil = pupil(1:100000,4)-tmp(1:100000);
    
    all_pup(:,isubj,iblock) = pupil;
    
    

  end
end

all_pup = all_pup(:,SUBJLIST,:);
%%
for isubj = 1 :28
  isubj
%   figure; set(gcf,'color','w');
%   subplot(1,2,1);
%   plot(all_pup(:,isubj,1));

starting_values = 0.0001:0.0001:0.01;
for i = 1 : length(starting_values)
  [tmp_l(i), fval(i)] = tp_fitexpdecay(1:100000,all_pup(:,isubj,1),starting_values(i));
end
[~,i_min]=sort(fval);
lamb(isubj,1) = tmp_l(i_min(1));

% X = [ones(1,length(all_pup(:,isubj,1)))' (1:100000)'./400];
%   Y = (all_pup(:,isubj,1)-all_pup(1,isubj,1))./all_pup(1,isubj,1);
%   tmp = X\Y; slp(isubj,1) = tmp(2);
%   subplot(1,2,2);
%   plot(all_pup(:,isubj,2));
%   X = [ones(1,length(all_pup(:,isubj,2)))' (1:100000)'./400];
%   Y = (all_pup(:,isubj,2)-all_pup(1,isubj,2))./all_pup(1,isubj,2);
%   tmp = X\Y; slp(isubj,2) = tmp(2);
starting_values = 0.0004:0.002:0.1;
for i = 1 : length(starting_values)
  [tmp_l(i), fval(i)] = tp_fitexpdecay(1:100000,all_pup(:,isubj,2),starting_values(i));
end
[~,i_min]=sort(fval);
lamb(isubj,2) = tmp_l(i_min(1));

end
%%
% 
% slp(5,:)=[]; slp(16,:)=[];
% 
figure; set(gcf,'color','w');
scatter(lamb(:,1),lamb(:,2))
% scatter(slp(:,1),slp(:,2))
% % axis([-0.035 0.035 -0.02 0.02])
tp_editplots
lsline; xlabel('Lambda (Block 1)'); ylabel('Lambda Block 2')
[~,p]=corr(lamb(:,1),lamb(:,2))

