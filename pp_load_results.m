function [plt_gla, plt_hh,plt_mue,plt_all]  = pp_load_results()

%% LOAD GLASGOW DATA

outdir = '~/pp/proc/src/';
ord    = pconn_randomization;
load(sprintf('~/pp/proc/pp_atlas_BNA.mat'))
load ~/standard_sourcemodel_BNA_5mm.mat

addpath('~/Documents/MATLAB/fieldtrip-20181231/')
ft_defaults

trans = pp_transfer_gla2hh;

v=3;
SUBJ = 1:24; SUBJ([5 9])=[];

for n_subj = 1: length(SUBJ)
  isubj = SUBJ(n_subj);
  for iblock = 1 : 1
    clear src_r
    try
      load(sprintf([outdir 'pp_gla_src_pupil_power_correlations_s%d_b%d_v%d.mat'],isubj,iblock,v));
      
      plt_gla.pow_sens(:,:,n_subj) = outp.tp_sens_pow;
      plt_gla.corr_sens(:,:,n_subj) = outp.sens_r;
      plt_gla.corr_src(:,:,n_subj) = outp.src_r(trans,:);
      plt_gla.corr_src_df(:,:,n_subj) = outp.src_r_df(trans,:);
      plt_gla.corr_src_filt(:,:,n_subj) = outp.src_r_filt(trans,:);
      plt_gla.corr_src_df_filt(:,:,n_subj) = outp.src_r_df_filt(trans,:);
    catch me
      plt_gla.pow_sens(:,:,n_subj) = nan(248,25);
      plt_gla.corr_sens(:,:,n_subj)= nan(248,25);
      plt_gla.corr_src(:,:,n_subj) = nan(8799,25);
      plt_gla.nai_src(:,:,n_subj) = nan(8799,25);
      plt_gla.corr_src_df(:,:,n_subj) = nan(8799,25);
      plt_gla.corr_src_filt(:,:,n_subj) = nan(8799,25);
      plt_gla.corr_src_df_filt(:,:,n_subj) = nan(8799,25);
      
      warning('!')
    end
  end
end
% 
for igrid = 1 : max(BNA.tissue_5mm(:))
  plt_gla.corr_src_BNA(igrid,:,:) = tanh(mean(atanh(plt_gla.corr_src(BNA.tissue_5mm == igrid,:,:))));
  plt_gla.corr_src_filt_BNA(igrid,:,:) = tanh(mean(atanh(plt_gla.corr_src_filt(BNA.tissue_5mm == igrid,:,:))));
  plt_gla.corr_src_df_BNA(igrid,:,:) = tanh(mean(atanh(plt_gla.corr_src_df(BNA.tissue_5mm == igrid,:,:))));
  plt_gla.corr_src_df_filt_BNA(igrid,:,:) = tanh(mean(atanh(plt_gla.corr_src_df_filt(BNA.tissue_5mm == igrid,:,:))));
end

plt_hh.pow_sens = nan(274,25,28,2);
plt_hh.corr_sens = nan(274,25,28,2);

SUBJLIST = [4 5 6 7 8 9 10 11 12 13 15 16 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34];
for isubj = SUBJLIST
  isubj
  for iblock = 1 : 2
    clear src_r
    try
      load(sprintf([outdir 'pp_src_pupil_power_correlations_s%d_b%d_v%d.mat'],isubj,iblock,v));
      load(sprintf('~/pconn/proc/preproc/pconn_chanidx_s%d_m%d.mat',isubj,find(ord(isubj,:)==1)))
      idx=logical(idx);
      isubj
      plt_hh.pow_sens(idx,:,isubj,iblock) = outp.tp_sens_pow;
      plt_hh.corr_sens(idx,:,isubj,iblock) = outp.sens_r;
      plt_hh.corr_src(:,:,isubj,iblock) = outp.src_r;
      plt_hh.corr_src_df(:,:,isubj,iblock) = outp.src_r_df;
      plt_hh.corr_src_filt(:,:,isubj,iblock) =outp.src_r_filt;
      plt_hh.corr_src_df_filt(:,:,isubj,iblock) = outp.src_r_df_filt;
    catch me
      warning('!!!')
      plt_hh.pow_sens(idx,:,isubj,iblock) = nan(sum(idx),25,1,1);
      plt_hh.corr_sens(idx,:,isubj,iblock) = nan(sum(idx),25,1,1);
      plt_hh.corr_src(:,:,isubj,iblock) = nan(8799,25);
%       plt_hh.nai_src(:,:,isubj,iblock) = nan(8799,25);
      plt_hh.corr_src_df(:,:,isubj) = nan(8799,25);
      plt_hh.corr_src_filt(:,:,isubj) = nan(8799,25);
      plt_hh.corr_src_df_filt(:,:,isubj) = nan(8799,25);
      continue
    end
  end
end

plt_hh.corr_sens = nanmean(plt_hh.corr_sens(:,:,SUBJLIST,:),4);
plt_hh.corr_src = nanmean(plt_hh.corr_src(:,:,SUBJLIST,:),4);
% plt_hh.nai_src = nanmean(plt_hh.nai_src(:,:,SUBJLIST,:),4);
plt_hh.pow_sens = nanmean(plt_hh.pow_sens(:,:,SUBJLIST,:),4);
plt_hh.corr_src_df = nanmean(plt_hh.corr_src_df(:,:,SUBJLIST,:),4);
plt_hh.corr_src_filt= nanmean(plt_hh.corr_src_filt(:,:,SUBJLIST,:),4);
plt_hh.corr_src_df_filt = nanmean(plt_hh.corr_src_df_filt(:,:,SUBJLIST,:),4);
% 
for igrid = 1 : max(BNA.tissue_5mm(:))
  plt_hh.corr_src_BNA(igrid,:,:) = tanh(mean(atanh(plt_hh.corr_src(BNA.tissue_5mm == igrid,:,:))));
  plt_hh.corr_src_df_BNA(igrid,:,:) = tanh(mean(atanh(plt_hh.corr_src_df(BNA.tissue_5mm == igrid,:,:))));
  plt_hh.corr_src_filt_BNA(igrid,:,:) = tanh(mean(atanh(plt_hh.corr_src_filt(BNA.tissue_5mm == igrid,:,:))));
  plt_hh.corr_src_df_filt_BNA(igrid,:,:) = tanh(mean(atanh(plt_hh.corr_src_df_filt(BNA.tissue_5mm == igrid,:,:))));
end

% LOAD MUENSTER DATA
% ----------------------------
for isubj = 1:37
  isubj
  for iblock = 1 : 1
    clear src_r
    try
      load(sprintf([outdir 'pp_mue_src_pupil_power_correlations_s%d_b%d_v%d.mat'],isubj,iblock,v));
%       idx=logical(idx);
      isubj
      plt_mue.pow_sens(:,:,isubj,iblock) = outp.tp_sens_pow;
      plt_mue.corr_sens(:,:,isubj,iblock) = outp.sens_r;
      plt_mue.corr_src(:,:,isubj,iblock) = outp.tp_src_r;
%       plt_mue.corr_src_df(:,:,isubj,iblock) = outp.src_r_df;
%       plt_mue.corr_src_filt(:,:,isubj,iblock) =outp.src_r_filt;
%       plt_mue.corr_src_df_filt(:,:,isubj,iblock) = outp.src_r_df_filt;
    catch me
      warning('!!!')
      plt_mue.pow_sens(:,:,isubj,iblock) = nan(275,25,1,1);
      plt_mue.corr_sens(:,:,isubj,iblock) = nan(275,25,1,1);
      plt_mue.corr_src(:,:,isubj,iblock) = nan(8799,25);
%       plt_hh.nai_src(:,:,isubj,iblock) = nan(8799,25);
%       plt_mue.corr_src_df(:,:,isubj) = nan(8799,25);
%       plt_mue.corr_src_filt(:,:,isubj) = nan(8799,25);
%       plt_mue.corr_src_df_filt(:,:,isubj) = nan(8799,25);
      continue
    end
  end
end

for igrid = 1 : max(BNA.tissue_5mm(:))
  plt_mue.corr_src_BNA(igrid,:,:) = tanh(mean(atanh(plt_mue.corr_src(BNA.tissue_5mm == igrid,:,:))));
%   plt_mue.corr_src_df_BNA(igrid,:,:) = tanh(mean(atanh(plt_mue.corr_src_df(BNA.tissue_5mm == igrid,:,:))));
%   plt_mue.corr_src_filt_BNA(igrid,:,:) = tanh(mean(atanh(plt_mue.corr_src_filt(BNA.tissue_5mm == igrid,:,:))));
%   plt_mue.corr_src_df_filt_BNA(igrid,:,:) = tanh(mean(atanh(plt_mue.corr_src_df_filt(BNA.tissue_5mm == igrid,:,:))));
end

% SENSOR LEVEL: SORTED BY SENSOR POSITION 
% ----------------------------
cfg=[];
cfg.layout='4D248.lay';
lay = ft_prepare_layout(cfg);
minmax_gla=[min(lay.pos(:,2)) max(lay.pos(:,2))];
ser_gla = linspace(minmax_gla(1),minmax_gla(2),40);
plt_gla.corr_sens_ord= zeros(size(ser_gla,2)-1,25,22);

for i = 1 : size(ser_gla,2)-1
  idx = lay.pos(:,2)<ser_gla(i+1) & lay.pos(:,2)>ser_gla(i);
  plt_gla.corr_sens_ord(i,:,:) = nanmean(plt_gla.corr_sens(idx,:,:),1);
end

cfg=[];
cfg.layout='CTF275.lay';
lay = ft_prepare_layout(cfg);
lay.pos([203 276 277],:)=[];
minmax_hh=[min(lay.pos(:,2)) max(lay.pos(:,2))];
ser_hh = linspace(minmax_hh(1),minmax_hh(2),40);
plt_hh.corr_sens_ord= zeros(size(ser_hh,2)-1,25,28);

for i = 1 : size(ser_gla,2)-1
  idx = lay.pos(:,2)<ser_hh(i+1) & lay.pos(:,2)>ser_hh(i);
  plt_hh.corr_sens_ord(i,:,:) = nanmean(plt_hh.corr_sens(idx,:,:),1);
end

for i = 1 : size(ser_gla,2)-1
  idx = lay.pos(:,2)<ser_hh(i+1) & lay.pos(:,2)>ser_hh(i);
  plt_mue.corr_sens_ord(i,:,:) = nanmean(plt_mue.corr_sens(idx,:,:),1);
end

% COLLAPSE ACROSS DATASETS
plt_all.corr_src = cat(3,plt_hh.corr_src,plt_gla.corr_src,plt_mue.corr_src);
plt_all.corr_src_BNA = cat(3,plt_hh.corr_src_BNA,plt_gla.corr_src_BNA,plt_mue.corr_src_BNA);
plt_all.corr_src_df = cat(3,plt_hh.corr_src_df,plt_gla.corr_src_df);
plt_all.corr_src_df_BNA = cat(3,plt_hh.corr_src_df_BNA,plt_gla.corr_src_df_BNA);


% -----------------------------
% LOAD CROSS FREQUENCY
% -----------------------------
% 
% v=3;
% 
% subj_idx=[];
% for isubj = SUBJLIST
%   isubj
%   for iblock = 1 : 2
%     clear src_r
%     try
%       load(sprintf([outdir 'pp_src_pupil_power_correlations_crossfreq_s%d_b%d_v%d.mat'],isubj,iblock,v));
%       subj_idx=[subj_idx isubj];
%       plt_hh.cf_corr(:,:,:,isubj,iblock)=src_r;
%     catch me
%       plt_hh.cf_corr(:,:,:,isubj,iblock)=nan(8799,25,21);
%     end
%   end
% end
% 
% subj_idx = unique(subj_idx);
% plt_hh.cf_corr=nanmean(plt_hh.cf_corr(:,:,:,subj_idx,:),5);

% for igrid = 1 : max(BNA.tissue_5mm(:))
%   plt_hh.cf_corr_BNA(igrid,:,:,:) = mean(plt_hh.cf_corr(BNA.tissue_5mm == igrid,:,:,:));
% end

% -----------------------------
% LOAD CROSS FREQUENCY
% -----------------------------

% v=3;
% 
% subj_idx=[];
% for isubj = SUBJLIST
%   isubj
%   for iblock = 1 : 2
%     clear src_r
%     try
%       load(sprintf([outdir 'pp_cnt_src_pupil_power_correlations_crossfreq_s%d_b%d_v%d.mat'],isubj,iblock,v));
%       subj_idx=[subj_idx isubj];
%       plt_hh.cf_corr_cnt(:,:,:,isubj,iblock)=src_r;
%     catch me
%       plt_hh.cf_corr_cnt(:,:,:,isubj,iblock)=nan(8799,25,21);
%     end
%   end
% end
% 
% subj_idx = unique(subj_idx);
% plt_hh.cf_corr=nanmean(plt_hh.cf_corr_cnt(:,:,:,subj_idx,:),5);
% 
% for igrid = 1 : max(BNA.tissue_5mm(:))
%   plt_hh.cf_corr_BNA_cnt(igrid,:,:,:) = mean(plt_hh.cf_corr_cnt(BNA.tissue_5mm == igrid,:,:,:));
% end
% 

% v=1;
% SUBJ = 1:24; SUBJ([5 9])=[];
% 
% for n_subj = 1: length(SUBJ)
%   isubj = SUBJ(n_subj);
%   for iblock = 1 : 1
%     clear src_r
%     try
%       load(sprintf([outdir 'pp_gla_src_pupil_power_correlations_crossfreq_s%d_b%d_v%d.mat'],isubj,iblock,v));
%       plt_gla.cf_corr(:,:,:,n_subj)=src_r(trans,:,:);
%     catch me
%       plt_gla.cf_corr(:,:,:,n_subj)=nan(8799,25,21);
%     end
%   end
% end
% for igrid = 1 : max(BNA.tissue_5mm(:))
%   plt_gla.cf_corr_BNA(igrid,:,:,:) = mean(plt_gla.cf_corr(BNA.tissue_5mm == igrid,:,:,:));
% end
