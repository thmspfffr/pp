%% pp_src_pupil_power_correlations
% correlate bandpass filtered (or via wavelets) pupil and MEG signals

clear
restoredefaultpath

% -------------------------
% VERSION 1: no pupil lag
% -------------------------
% v = 1;
% SUBJLIST = [4 5 6 7 8 9 10 11 12 13 15 16 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34];
% lag = 0;
% -------------------------
% VERSION 3: with pupil lag
% -------------------------
v = 2;
SUBJLIST = [4 5 6 7 8 9 10 11 12 13 15 16 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34];
lag = 1;
% -------------------------

addpath ~/Documents/MATLAB/fieldtrip-20160919/
addpath ~/pconn/matlab/
load(sprintf('~/pp/proc/pp_atlas_BNA.mat'))

ft_defaults

outdir = '~/pp/proc/src/';
ord    = pconn_randomization;

freqoi=2.^(1:(1/4):7); % 2-128 Hz as per Hipp et al. (2012) Nat Neurosci

%%
% -------------------------
for isubj = SUBJLIST
  
  % identify placebo condition (ord==1)
  im = find(ord(isubj,:)==1);
  
  for iblock = 1:2
    %
    fn = sprintf('pp_src_pupil_power_correlations_s%d_b%d_v%d',isubj,iblock,v);
    %     if tp_parallel(fn,outdir,1,0)
    %       continue
    %     end
    %
    fprintf('Processing subj%d block%d ...\n',isubj,iblock);
    
    try
      % load cleaned meg data
      load(sprintf('~/pp/data/ham/pupmod_rest_sens_cleandat_s%d_m%d_b%d_v%d.mat',isubj,im,iblock,1))
    catch me
      continue
    end
    
    cfg=[];
    cfg.layout='CTF275.lay';
    lay = ft_prepare_layout(cfg);
    [~, outp.chanidx] = ismember(lay.label(1:275),label);
    
    % bp-filter and resample pupil
    % ------
    k = 2; f_sample = 1000;
    fnq = f_sample/2;
    hil_hi = 0.005; hil_lo = 2;
    hil_Wn=[hil_hi/fnq hil_lo/fnq];
    [bhil, ahil] = butter(k, hil_Wn);
    
    pupil = filtfilt(bhil, ahil, pupil);
    pupil = resample(pupil,400,1000);
    % ------
    f_sample = 400;
    % align pupil and meg (at signal offset)
    % ------

    dat = dat(:,end:-1:1);
    
    len = min([size(pupil,1) size(dat,2)]);
    if len/400 > 600
      len = 400*600;
    end
    
    dat = dat(:,1:len);
    pupil = pupil(1:len);
    
    dat = dat(:,end:-1:1);
    pupil = pupil(end:-1:1);
    % ------
    
    % pupil shift: 930 ms from hoeks & levelt (1992)
    if lag
      pup_shift = round(400*0.93);
      pupil = pupil(pup_shift:end); pupil(end+1:end+pup_shift-1)=nan;
    end
    
    pupil_df = diff(pupil);
    
    dat(:,isnan(pupil))=nan(size(dat,1),sum(isnan(pupil)));
    
    load(['/home/tpfeffer/pp/proc/src/' sprintf('pp_sa_s%d_m%d_b%d_v%d.mat',isubj,im,iblock,1)],'sa');
    
    lf = sa.L_genemaps_aal;
    
    save(sprintf('~/pp/proc/sens/pp_meg_pupil_lfs_s%d_b%d.mat',isubj,iblock),'dat','pupil','lf');
    
    [outp.pxx,outp.fxx]=pwelch(dat(:,~isnan(dat(1,:)))',hanning(400),0,1:1:200,400);
    
    outp.sens_pow    = nan(275,25);
    outp.sens_r      = nan(275,25);
    outp.sens_r_df   = nan(275,25);
    outp.sens_mi_df  = nan(275,25);
    outp.sens_mi     = nan(275,25);
    outp.sens_mi0    = nan(275,25);
    outp.sens_mi_df0 = nan(275,25);
    
    for ifreq=1:length(freqoi)
      ifreq
      
      para          = [];
      para.freq     = freqoi(ifreq);
      para.fsample  = f_sample;
      para.overlap  = 0.8;
      [csd, dataf,opt]=tp_compute_csd_wavelets(dat,para);
      
      tmp = diag(abs(csd));
      outp.sens_pow(outp.chanidx>0,ifreq) = tmp(outp.chanidx(outp.chanidx>0));
      
      % -------------------------------
      % prepare pupil signal
      % -------------------------------
      % take pupil_df as signal is one sample shorter
      % taking pupil or dat will result in an error
      nseg=floor((size(pupil_df,1)-opt.n_win)/opt.n_shift+1);
      
      pup = nan(nseg,1);
      pup_df = nan(nseg,1);
      for j=1:nseg
        tmp = pupil((j-1)*opt.n_shift+1:(j-1)*opt.n_shift+opt.n_win);
        tmp2 = pupil_df((j-1)*opt.n_shift+1:(j-1)*opt.n_shift+opt.n_win);
        pup(j) = mean(tmp.*gausswin(opt.n_win,3));
        pup_df(j) = mean(tmp2.*gausswin(opt.n_win,3));
      end
      
      % find indices of non-NAN segments
      idx_valid = find(~isnan(dataf(1,:)));
      idx_valid_df = find(~isnan(pup_df)');
      idx = intersect(idx_valid,idx_valid_df);
      
      env = abs(dataf).^2;
      
      % correlate pupil with sensor level signal
      outp.sens_r(outp.chanidx>0,ifreq) = corr(pup(idx),env(outp.chanidx(outp.chanidx>0),idx)','type','Spearman');
      outp.sens_r_df(outp.chanidx>0,ifreq) = corr(pup_df(idx),env(outp.chanidx(outp.chanidx>0),idx)','type','Spearman');
      
      % -------------------------------
      % sensor-level mutual information
      % -------------------------------
      cnpup     = copnorm(pup(idx));
      cnpup0    = cnpup(end:-1:1);
      cnpup_df  = copnorm(pup_df(idx));
      cnpup_df0 = cnpup_df(end:-1:1);
      cnpow     = copnorm(env(outp.chanidx(outp.chanidx>0),idx))';
      tmp = []; tmp_df = []; tmp0 = []; tmp_df0 = [];
      for isens = 1 : size(cnpow,2)
        tmp(isens)    = mi_gg_dfi_ak(cnpow(:,isens),cnpup,[]);
        tmp0(isens)   = mi_gg_dfi_ak(cnpow(:,isens),cnpup0,[]);
        tmp_df(isens) = mi_gg_dfi_ak(cnpow(:,isens),cnpup_df,[]);
        tmp_df0(isens) = mi_gg_dfi_ak(cnpow(:,isens),cnpup_df0,[]);
      end
      outp.sens_mi(outp.chanidx>0,ifreq)    = tmp;
      outp.sens_mi0(outp.chanidx>0,ifreq)   = tmp0;
      outp.sens_mi_df(outp.chanidx>0,ifreq) = tmp_df;
      outp.sens_mi_df0(outp.chanidx>0,ifreq) = tmp_df0;
      % -------------------------------
      % sensor-level cross correlation
      % -------------------------------
      nlags=floor(10/(opt.n_shift/f_sample)); % roughly 10s
      for isens = 1 : size(env,1)
        tmp_pup = pup(idx)-mean(pup(idx));
        tmp_pup_df = pup_df(idx)-nanmean(pup_df(idx));
        tmp_env = env(isens,idx)-nanmean(env(isens,idx),2);
        %           prod_std(isens) = std(tmp_pup)*std(tmp_env);
        [outp.xcorr{ifreq}(:,isens),lags] = xcorr(tmp_pup,tmp_env,nlags,'coeff');
        [outp.xcorr_df{ifreq}(:,isens),lags] = xcorr(tmp_pup_df,tmp_env,nlags,'coeff');
      end
      outp.xcorr{ifreq}(:,outp.chanidx>0) = outp.xcorr{ifreq}(:,outp.chanidx(outp.chanidx>0));
      outp.xcorr{ifreq}(:,outp.chanidx==0)= nan;
      outp.xcorr_df{ifreq}(:,outp.chanidx>0) = outp.xcorr_df{ifreq}(:,outp.chanidx(outp.chanidx>0));
      outp.xcorr_df{ifreq}(:,outp.chanidx==0)= nan;
      
      lags = lags*(opt.n_shift/f_sample);
      outp.xcorr_lags{ifreq} = lags;
      % -------------------------------
      % beamforming (ignore name of leadfield)
      % -------------------------------
      para          = [];
      para.reg      = 0.05;
      [filt,pow] = tp_beamformer(real(csd),sa.L_genemaps_aal,para);
      % --------------
      % beamform again with noise to compute "NAI"
      % --------------
      Lr = reshape(sa.L_genemaps_aal,[size(sa.L_genemaps_aal,1) 8799*3]); csd_noise = Lr*Lr';
      [~,noise]      = tp_beamformer(real(csd_noise),sa.L_genemaps_aal,para);
      outp.src_nai(:,ifreq) = pow./noise;
      % -------------------------------
      
      % compute power & lowpass filter
      src_pow = abs(filt'*dataf(:,idx)).^2;
      
      % correlate with pupil
      outp.src_r(:,ifreq) = corr(pup(idx),src_pow','type','Spearman');
      outp.src_r_df(:,ifreq) = corr(pup_df(idx),src_pow','type','Spearman');
      
      % mutual information
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
      
      clear src pup pup_df
      
    end
    
    save([outdir fn '.mat'],'outp')
    tp_parallel(fn,outdir,0)
    
    clear src_r all_nai outp
    
  end
end

exit
