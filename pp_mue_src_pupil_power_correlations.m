%% pp_mue_src_pupil_power_correlations
% correlate bandpass filtered (or via wavelets) pupil and MEG signals

clear
restoredefaultpath

% -------------------------
% VERSION 1: no pupil lag
% -------------------------
v = 1;
freqoi    = 2.^(1:(1/4):7);
lag = 0;
% -------------------------
% VERSION 3: with pupil lag
% -------------------------
% v = 2;
% freqoi    = 2.^(1:(1/4):7);
% lag = 1;
% -------------------------

addpath('~/Documents/MATLAB/fieldtrip-20181231/')
addpath ~/pconn/matlab/
load(sprintf('~/pp/proc/pp_atlas_BNA.mat'))

ft_defaults

outdir = '~/pp/proc/src/';
ord    = pconn_randomization;

addpath /home/gnolte/meth/highlevel/

d=dir('~/pp/data_gla/fw4bt/osfstorage/data/ms01/meg/*mat');

SUBJLIST = [];
for i = 1 : length(d) 
    SUBJLIST = [SUBJLIST; d(i).name(end-7:end-4)];  
end

%%
% -------------------------
for isubj = 1:size(SUBJLIST,1)
  
  clear data dat pupil pup dataf src_r
  
  for iblock = 1:1
    %
    fn = sprintf('pp_mue_src_pupil_power_correlations_s%d_b%d_v%d',isubj,iblock,v);
    if tp_parallel(fn,outdir,1,0)
      continue
    end
    %
    fprintf('Processing subj%d block%d ...\n',isubj,iblock);
    
    try
      % load pupil data
      load(sprintf('~/pp/data_gla/fw4bt/osfstorage/data/ms01/pupil/%s_pupil_preproc_lp4.mat',SUBJLIST(isubj,:)))
      pupil = reshape(data.trial{1}(4,:),[size(data.trial{1}(4,:),2) 1]);
      
      % load meg data
      load(sprintf('~/pp/data_gla/fw4bt/osfstorage/data/ms01/meg/cleanmeg_%s.mat',SUBJLIST(isubj,:)))
      data = cleanmeg; clear cleanmeg
      f_sample = data.fsample;
   
      cfg=[];
      cfg.layout='CTF275.lay';
      lay = ft_prepare_layout(cfg);
      [~, outp.chanidx] = ismember(lay.label(1:275),data.label);
    
    catch me
      src_r = nan(246,25);
      save([outdir fn '.mat'],'src_r')
      continue
    end
    
    cfg = [];
    cfg.resamplefs = 400;
    data = ft_resampledata(cfg,data);
    f_sample = data.fsample;
    
    k = 2;
    fnq = f_sample/2;
    hil_hi = 0.005;
    hil_lo = 2;
    hil_Wn=[hil_hi/fnq hil_lo/fnq];
    [bhil, ahil] = butter(k, hil_Wn);
    
    pupil = filtfilt(bhil, ahil, pupil);
    pupil = resample(pupil,400,600);
    
    if lag
        pup_shift = round(f_sample*0.93); % 930s from hoeks and levelt (1992?)
        pupil = pupil(pup_shift:end); pupil(end+1:end+pup_shift-1)=nan;
    end
    pupil_df = diff(pupil);
    
    data.trial{1}(:,isnan(pupil))=nan(size(data.trial{1},1),sum(isnan(pupil)));
    
%     tmp = data;
    data.avg = data.trial{1}'; %data.trial{1} = [];

    load(sprintf('~/pp/data_gla/fw4bt/osfstorage/data/ms01/leadfields/lf_%s.mat',SUBJLIST(isubj,:)))
    outp.sens_pow   = nan(275,25);
    outp.sens_r     = nan(275,25);
    outp.sens_r_df  = nan(275,25);
    outp.sens_mi    = nan(275,25);
    outp.sens_mi_df = nan(275,25);
    
    for ifreq=1:numel(freqoi)
      
      fprintf('Freq: %d\n',ifreq)

      % -------------------------------
      % compute csd
      % -------------------------------
      para          = [];
      para.freq     = freqoi(ifreq);
      para.fsample  = f_sample;  
      para.overlap = 0.8;
      [csd, dataf,opt]=tp_compute_csd_wavelets(data.avg',para);
      
      % -------------------------------
      % compute power from csd
      % -------------------------------
      tmp = diag(abs(csd)); 
      outp.sens_pow(outp.chanidx>0,ifreq) = tmp(outp.chanidx(outp.chanidx>0));
      % -------------------------------
      % prepare pupil signal
      % -------------------------------
      nseg=floor((size(data.avg,1)-opt.n_win)/opt.n_shift+1);
      pup = nan(nseg,1);
      pup_df = nan(nseg,1);
      for j=1:nseg
        tmp = pupil((j-1)*opt.n_shift+1:(j-1)*opt.n_shift+opt.n_win);
        pup(j) = mean(tmp.*gausswin(opt.n_win,3));
        if j==nseg
            tmp2 = pupil_df((j-1)*opt.n_shift+1:(j-1)*opt.n_shift+opt.n_win-1);
            pup_df(j) = mean(tmp2.*gausswin(opt.n_win-1,3));
        else
            tmp2 = pupil_df((j-1)*opt.n_shift+1:(j-1)*opt.n_shift+opt.n_win);
            pup_df(j) = mean(tmp2.*gausswin(opt.n_win,3));
        end
      end

      % -------------------------------
      % beamforming
      % -------------------------------
      para      = [];
      para.iscs = 1;
      para.reg  = 0.05;
      tp_filt   = tp_beamformer(real(csd),lf,para);
      % -------------------------------
      % identify artifactual segments
      idx = find(~isnan(dataf(1,:))'&~isnan(pup)&~isnan(pup_df));
      % -------------------------------
      % sensor-level pupil power correlation
      % -------------------------------
      env = (abs(dataf(outp.chanidx(outp.chanidx>0),idx)).^2)';
      outp.sens_r(outp.chanidx>0,ifreq) = corr(pup(idx),env,'type','spearman');
      outp.sens_r_df(outp.chanidx>0,ifreq) = corr(pup_df(idx),env,'type','spearman');
      % -------------------------------
      % sensor-level mutual information
      % -------------------------------
      cnpup     = copnorm(pup(idx));
      cnpup_df  = copnorm(pup_df(idx));
      cnpow     = copnorm(env);
      tmp = []; tmp_df = [];
      for isens = 1 : size(dataf,1)
        tmp(isens)    = mi_gg_dfi_ak(cnpow(:,isens),cnpup,[]);
        tmp_df(isens) = mi_gg_dfi_ak(cnpow(:,isens),cnpup_df,[]);
      end
      outp.sens_mi(outp.chanidx>0,ifreq) = tmp;
      outp.sens_mi_df(outp.chanidx>0,ifreq) = tmp_df;
      % -------------------------------
      % sensor-level cross correlation
      % -------------------------------
      nlags=floor(10/(opt.n_shift/f_sample)); % roughly 10s
      for isens = 1 : size(env,2)
        tmp_pup = pup(idx)-mean(pup(idx));
        tmp_pup_df = pup_df(idx)-mean(pup_df(idx));
        tmp_env = env-nanmean(env,1);
        [outp.xcorr{ifreq}(:,isens),lags] = xcorr(tmp_pup(idx),tmp_env(:,isens),nlags,'coeff');
        [outp.xcorr_df{ifreq}(:,isens),lags] = xcorr(tmp_pup_df(idx),tmp_env(:,isens),nlags,'coeff');
      end
      outp.xcorr{ifreq}(:,outp.chanidx>0) = outp.xcorr{ifreq}(:,outp.chanidx(outp.chanidx>0));
      outp.xcorr{ifreq}(:,outp.chanidx==0)= nan;
      outp.xcorr_df{ifreq}(:,outp.chanidx>0) = outp.xcorr_df{ifreq}(:,outp.chanidx(outp.chanidx>0));
      outp.xcorr_df{ifreq}(:,outp.chanidx==0)= nan;
      
      lags=lags*(opt.n_shift/f_sample);
      outp.xcorr_lags{ifreq} = lags;
      % -------------------------------
      % correlate power with pupil
      % -------------------------------
      src_pow = abs(tp_filt'*dataf(:,idx)).^2; % source level power fluct
      outp.src_r(:,ifreq) = corr(pup(idx),src_pow','type','spearman');
      outp.src_r_df(:,ifreq) = corr(pup_df(idx),src_pow','type','spearman');
      % -------------------------------
      % source level mutual information
      % -------------------------------
      cnpup     = copnorm(pup(idx));
      cnpup_df  = copnorm(pup_df(idx));
      cnpow     = copnorm(src_pow)';
      for isrc = 1 : size(src_pow,1)
        outp.src_mi(isrc,ifreq) = mi_gg_dfi_ak(cnpow(:,isrc),cnpup,[]);
        outp.src_mi_df(isrc,ifreq) = mi_gg_dfi_ak(cnpow(:,isrc),cnpup_df,[]);
      end
      % -------------------------------
      % source-level cross correlation
      % -------------------------------
      
%       nlags=floor(10/(opt.n_shift/f_sample)); % roughly 10s      
%       for isrc = 1 : size(src_pow,1)
%         tmp_pup = pup(idx)-mean(pup(idx));
%         tmp_pup_df = pup_df(idx)-nanmean(pup_df(idx));
%         tmp_env = src_pow(isrc,:)-nanmean(src_pow(isrc,:),2);
%         [outp.src_xcorr{ifreq}(:,isrc)] = xcorr(tmp_pup,tmp_env,nlags,'coeff');
%         [outp.src_xcorr_df{ifreq}(:,isrc)] = xcorr(tmp_pup_df,tmp_env,nlags,'coeff');
%       end
      
      clear src pup csd tp_csd dataf 
      
    end

    save([outdir fn '.mat'],'outp')
    tp_parallel(fn,outdir,0)
    
    clear outp
    
  end
end

exit
