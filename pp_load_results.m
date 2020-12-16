function [plt_gla, plt_hh,plt_mue,plt_hh_cnt,plt_all]  = pp_load_results(v)

%% LOAD GLASGOW DATA

outdir = '~/pp/proc/src/';
ord    = pconn_randomization;
load(sprintf('~/pp/proc/pp_atlas_BNA.mat'))
load ~/standard_sourcemodel_BNA_5mm.mat

addpath ~/Documents/MATLAB/fieldtrip-20190224/
ft_defaults

trans = pp_transfer_gla2hh;

SUBJ = 1:24; SUBJ([5 9])=[];

freqoi    = 2.^(1:(1/4):7);
for ifreq = 1 : length(freqoi)
  [~,opt]=tp_mkwavelet(freqoi(ifreq),0.5,400,0.8);
  nlags(ifreq) = floor(10/(opt.n_shift/400))*2+1;
end

for n_subj = 1: length(SUBJ)
  isubj = SUBJ(n_subj);
  fprintf('Glasgow: Subj%d\n',isubj)
  load(sprintf([outdir 'pp_gla_src_pupil_power_correlations_s%d_b1_v%d.mat'],isubj,v));
  plt_gla.pow_sens(:,:,n_subj)      = outp.sens_pow;
  plt_gla.corr_sens(:,:,n_subj)     = outp.sens_r;
  plt_gla.corr_src(:,:,n_subj)      = outp.src_r(trans,:);
  plt_gla.corr_src_df(:,:,n_subj)   = outp.src_r_df(trans,:);
  %       plt_gla.sens_mi(:,:,n_subj)       = outp.sens_mi;
  %       plt_gla.sens_mi0(:,:,n_subj)      = outp.sens_mi0;
  plt_gla.src_mi(:,:,n_subj)        = outp.src_mi;
  for ifreq = 1 : 25
    plt_gla.xcorr{ifreq}(:,:,n_subj)      = outp.xcorr{ifreq};
    plt_gla.xcorr_df{ifreq}(:,:,n_subj)  	= outp.xcorr_df{ifreq};
    plt_gla.xcorr_lags{ifreq}(:,:,n_subj)	= outp.xcorr_lags{ifreq};
  end
  clear outp
end
% 
for igrid = 1 : max(BNA.tissue_5mm(:))
  plt_gla.corr_src_BNA(igrid,:,:) = tanh(mean(atanh(plt_gla.corr_src(BNA.tissue_5mm == igrid,:,:))));
  plt_gla.corr_src_df_BNA(igrid,:,:) = tanh(mean(atanh(plt_gla.corr_src_df(BNA.tissue_5mm == igrid,:,:))));
%   plt_gla.src_mi_BNA(igrid,:,:)    = mean(plt_gla.src_mi(BNA.tissue_5mm == igrid,:,:));
end

SUBJLIST = [4 5 6 7 8 9 10 11 12 13 15 16 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34];

for ifreq = 1: length(freqoi)
  plt_hh.xcorr{ifreq} = nan(nlags(ifreq),275,34,2);
end

for isubj = SUBJLIST
  fprintf('Hamburg: Subj%d\n',isubj)
  for iblock = 1 : 2
    try
      load(sprintf([outdir 'pp_src_pupil_power_correlations_s%d_b%d_v%d.mat'],isubj,iblock,v));
      plt_hh.pow_sens(:,:,isubj,iblock)     = outp.sens_pow;
      plt_hh.corr_sens(:,:,isubj,iblock)    = outp.sens_r;
      plt_hh.corr_src(:,:,isubj,iblock)     = outp.src_r;
      plt_hh.corr_src_df(:,:,isubj,iblock)  = outp.src_r_df;
%       plt_hh.sens_mi(:,:,isubj,iblock)      = outp.sens_mi;
%       plt_hh.sens_mi0(:,:,isubj,iblock)     = outp.sens_mi0;
      plt_hh.src_mi(:,:,isubj,iblock)       = outp.src_mi;
      for ifreq = 1 : 25
        plt_hh.xcorr{ifreq}(:,:,isubj,iblock)       = outp.xcorr{ifreq};
        plt_hh.xcorr_df{ifreq}(:,:,isubj,iblock)  	= outp.xcorr_df{ifreq};
        plt_hh.xcorr_lags{ifreq}(:,:,isubj,iblock)	= outp.xcorr_lags{ifreq};
      end
      clear outp
    catch me
      warning('!!!')
      plt_hh.pow_sens(:,:,isubj,iblock)     = nan(275,25,1,1);
      plt_hh.corr_sens(:,:,isubj,iblock)    = nan(275,25,1,1);
      plt_hh.corr_src(:,:,isubj,iblock)     = nan(8799,25);
      plt_hh.corr_src_df(:,:,isubj)         = nan(8799,25);
%       plt_hh.sens_mi(:,:,isubj,iblock)      = nan(275,25,1,1);
%       plt_hh.sens_mi0(:,:,isubj,iblock)     = nan(275,25,1,1);
%       plt_hh.src_mi(:,:,isubj,iblock)       = nan(8799,25);
      continue
    end
  end
end

plt_hh.pow_sens     = nanmean(plt_hh.pow_sens(:,:,SUBJLIST,:),4);
% plt_hh.sens_mi      = nanmean(plt_hh.sens_mi(:,:,SUBJLIST,:),4);
% plt_hh.sens_mi0     = nanmean(plt_hh.sens_mi0(:,:,SUBJLIST,:),4);
plt_hh.corr_sens    = nanmean(plt_hh.corr_sens(:,:,SUBJLIST,:),4);
plt_hh.corr_src     = nanmean(plt_hh.corr_src(:,:,SUBJLIST,:),4);
plt_hh.corr_src_df  = nanmean(plt_hh.corr_src_df(:,:,SUBJLIST,:),4);
% plt_hh.src_mi       = nanmean(plt_hh.src_mi(:,:,SUBJLIST,:),4);

for ifreq = 1: length(freqoi)
  plt_hh.xcorr{ifreq} = nanmean(plt_hh.xcorr{ifreq}(:,:,SUBJLIST,:),4);
end

for igrid = 1 : max(BNA.tissue_5mm(:))
  plt_hh.corr_src_BNA(igrid,:,:) = tanh(mean(atanh(plt_hh.corr_src(BNA.tissue_5mm == igrid,:,:))));
  plt_hh.corr_src_df_BNA(igrid,:,:) = tanh(mean(atanh(plt_hh.corr_src_df(BNA.tissue_5mm == igrid,:,:))));
  plt_hh.src_mi_BNA(igrid,:,:)    = mean(plt_hh.src_mi(BNA.tissue_5mm == igrid,:,:));
end

for isubj = SUBJLIST
  fprintf('Hamburg CNT: Subj%d\n',isubj)
  for iblock = 1 : 2
    try
      load(sprintf([outdir 'pp_cnt_src_pupil_power_correlations_s%d_b%d_v%d.mat'],isubj,iblock,v));
      plt_hh_cnt.pow_sens(:,:,isubj,iblock)     = outp.sens_pow;
      plt_hh_cnt.corr_sens(:,:,isubj,iblock)    = outp.sens_r;
      plt_hh_cnt.corr_src(:,:,isubj,iblock)     = outp.src_r;
      plt_hh_cnt.corr_src_df(:,:,isubj,iblock)  = outp.src_r_df;
%       plt_hh_cnt.sens_mi(:,:,isubj,iblock)      = outp.sens_mi;
%       plt_hh_cnt.sens_mi0(:,:,isubj,iblock)     = outp.sens_mi0;
%       plt_hh_cnt.src_mi(:,:,isubj,iblock)       = outp.src_mi;
      for ifreq = 1 : 25
        plt_hh_cnt.xcorr{ifreq}(:,:,isubj,iblock)       = outp.xcorr{ifreq};
        plt_hh_cnt.xcorr_df{ifreq}(:,:,isubj,iblock)    = outp.xcorr_df{ifreq};
        plt_hh_cnt.xcorr_lags{ifreq}(:,:,isubj,iblock)	= outp.xcorr_lags{ifreq};
      end
      clear outp
    catch me
      warning('!!!')
      plt_hh_cnt.pow_sens(:,:,isubj,iblock)     = nan(275,25,1,1);
      plt_hh_cnt.corr_sens(:,:,isubj,iblock)    = nan(275,25,1,1);
      plt_hh_cnt.corr_src(:,:,isubj,iblock)     = nan(8799,25);
      plt_hh_cnt.corr_src_df(:,:,isubj)         = nan(8799,25);
%       plt_hh_cnt.sens_mi(:,:,isubj,iblock)      = nan(275,25,1,1);
%       plt_hh_cnt.sens_mi0(:,:,isubj,iblock)     = nan(275,25,1,1);
%       plt_hh_cnt.src_mi(:,:,isubj,iblock)       = nan(8799,25);
      continue
    end
  end
end

for ifreq = 1: length(freqoi)
  plt_hh_cnt.xcorr{ifreq} = nanmean(plt_hh_cnt.xcorr{ifreq}(:,:,SUBJLIST,:),4);
end

plt_hh_cnt.pow_sens     = nanmean(plt_hh_cnt.pow_sens(:,:,SUBJLIST,:),4);
% plt_hh_cnt.sens_mi      = nanmean(plt_hh_cnt.sens_mi(:,:,SUBJLIST,:),4);
% plt_hh_cnt.sens_mi0     = nanmean(plt_hh_cnt.sens_mi0(:,:,SUBJLIST,:),4);
plt_hh_cnt.corr_sens    = nanmean(plt_hh_cnt.corr_sens(:,:,SUBJLIST,:),4);
plt_hh_cnt.corr_src     = nanmean(plt_hh_cnt.corr_src(:,:,SUBJLIST,:),4);
plt_hh_cnt.corr_src_df  = nanmean(plt_hh_cnt.corr_src_df(:,:,SUBJLIST,:),4);
% plt_hh_cnt.src_mi       = nanmean(plt_hh_cnt.src_mi(:,:,SUBJLIST,:),4);


for igrid = 1 : max(BNA.tissue_5mm(:))
  plt_hh_cnt.corr_src_BNA(igrid,:,:) = tanh(mean(atanh(plt_hh_cnt.corr_src(BNA.tissue_5mm == igrid,:,:))));
  plt_hh_cnt.corr_src_df_BNA(igrid,:,:) = tanh(mean(atanh(plt_hh_cnt.corr_src_df(BNA.tissue_5mm == igrid,:,:))));
%   plt_hh_cnt.src_mi_BNA(igrid,:,:)    = mean(plt_hh_cnt.src_mi(BNA.tissue_5mm == igrid,:,:));
end

% LOAD MUENSTER DATA
% ----------------------------
% list = [10 17 19 22 27 35  40]

SUBJLIST=1:41; SUBJLIST([10,12,17,19,22,27,35,38,39,40])=[];
for n_subj = 1:length(SUBJLIST)
  isubj = SUBJLIST(n_subj);
  fprintf('Muenster: Subj%d\n',isubj)
  
  load(sprintf([outdir 'pp_mue_src_pupil_power_correlations_s%d_b1_v%d.mat'],isubj,v));
  
  plt_mue.pow_sens(:,:,n_subj)     = outp.sens_pow;
  plt_mue.corr_sens(:,:,n_subj)    = outp.sens_r;
  plt_mue.corr_src(:,:,n_subj)     = outp.src_r;
  plt_mue.corr_src_df(:,:,n_subj)  = outp.src_r_df;
  %       plt_mue.sens_mi(:,:,n_subj)      = outp.sens_mi;
  %       plt_mue.src_mi(:,:,n_subj)       = outp.src_mi;
  %       plt_mue.sens_mi0(:,:,n_subj)     = outp.sens_mi0;
  for ifreq = 1 : 25
    plt_mue.xcorr{ifreq}(:,:,n_subj)       = outp.xcorr{ifreq};
    plt_mue.xcorr_df{ifreq}(:,:,n_subj)    = outp.xcorr_df{ifreq};
    plt_mue.xcorr_lags{ifreq}(:,:,n_subj)  = outp.xcorr_lags{ifreq};
  end
  clear outp
end
end

for igrid = 1 : max(BNA.tissue_5mm(:))
  plt_mue.corr_src_BNA(igrid,:,:) = tanh(mean(atanh(plt_mue.corr_src(BNA.tissue_5mm == igrid,:,:))));
  plt_mue.corr_src_df_BNA(igrid,:,:) = tanh(mean(atanh(plt_mue.corr_src_df(BNA.tissue_5mm == igrid,:,:))));
%   plt_mue.src_mi_BNA(igrid,:,:)    = mean(plt_mue.src_mi(BNA.tissue_5mm == igrid,:,:));
end

% ------------------------
% ORDER SENSORS FROM FRONT TO BACK FOR GLA
% -----------------------
cfg=[];
cfg.layout='4D248.lay';
lay = ft_prepare_layout(cfg);
minmax_gla=[min(lay.pos(1:248,2)) max(lay.pos(1:248,2))];
ser_gla = linspace(minmax_gla(1),minmax_gla(2),40);

for i = 1 : size(ser_gla,2)-1
  idx = lay.pos(:,2)<ser_gla(i+1) & lay.pos(:,2)>ser_gla(i);
  plt_gla.corr_sens_ord(i,:,:) = nanmean(plt_gla.corr_sens(idx,:,:),1);
%   plt_gla.sens_mi_ord(i,:,:) = nanmean(plt_gla.sens_mi(idx,:,:),1);
%   plt_gla.sens_mi_ord0(i,:,:) = nanmean(plt_gla.sens_mi0(idx,:,:),1);
  for ifreq = 1: length(freqoi)
    plt_gla.xcorr_ord{ifreq}(:,i,:) = nanmean(plt_gla.xcorr{ifreq}(:,idx,:),2);
  end
end

% ------------------------
% ORDER SENSORS FROM FRONT TO BACK FOR HH, HH TASK and MUE
% -----------------------
cfg=[];
cfg.layout='CTF275.lay';
lay = ft_prepare_layout(cfg);
minmax_hh=[min(lay.pos(1:275,2)) max(lay.pos(1:275,2))];
ser_hh = linspace(minmax_hh(1),minmax_hh(2),40);

for i = 1 : size(ser_hh,2)-1
  idx = lay.pos(:,2)>=ser_hh(i) & lay.pos(:,2)<=ser_hh(i+1);
  plt_hh_cnt.corr_sens_ord(i,:,:)   = nanmean(plt_hh_cnt.corr_sens(idx,:,:),1);
  plt_mue.corr_sens_ord(i,:,:)      = nanmean(plt_mue.corr_sens(idx,:,:),1);
  plt_hh.corr_sens_ord(i,:,:)   = nanmean(plt_hh.corr_sens(idx,:,:),1);
%   plt_hh_cnt.sens_mi_ord(i,:,:)     = nanmean(plt_hh_cnt.sens_mi(idx,:,:),1);
%   plt_hh_cnt.sens_mi_ord0(i,:,:)    = nanmean(plt_hh_cnt.sens_mi0(idx,:,:),1);
  for ifreq = 1: length(freqoi)
    plt_hh.xcorr_ord{ifreq}(:,i,:) = nanmean(plt_hh.xcorr{ifreq}(:,idx,:),2);
    plt_hh_cnt.xcorr_ord{ifreq}(:,i,:) = nanmean(plt_hh_cnt.xcorr{ifreq}(:,idx,:),2);
    plt_mue.xcorr_ord{ifreq}(:,i,:) = nanmean(plt_mue.xcorr{ifreq}(:,idx,:),2);
  end
end

% COLLAPSE ACROSS DATASETS
plt_all.corr_src = cat(3,plt_hh.corr_src,plt_gla.corr_src,plt_mue.corr_src);
plt_all.corr_src_BNA = cat(3,plt_hh.corr_src_BNA,plt_gla.corr_src_BNA,plt_mue.corr_src_BNA);
plt_all.corr_src_df = cat(3,plt_hh.corr_src_df,plt_gla.corr_src_df,plt_mue.corr_src_df);
plt_all.corr_src_df_BNA = cat(3,plt_hh.corr_src_df_BNA,plt_gla.corr_src_df_BNA,plt_mue.corr_src_df_BNA);




