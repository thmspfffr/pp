%% pp_gla_src_pupil_power_correlations
% correlate bandpass filtered (or via wavelets) pupil and MEG signals

clear
restoredefaultpath

% -------------------------
% VERSION 1: no pupil lag
% -------------------------
v = 1;
SUBJLIST  = 1:24;
freqoi    = 2.^(1:(1/4):7);
lag = 0;
% -------------------------
% VERSION 3: with pupil lag
% -------------------------
% v = 2;
% SUBJLIST  = 1:24;
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
%%
% -------------------------
for isubj = SUBJLIST
  
  clear data dat pupil pup dataf src_r
  
  for iblock = 1:1
    %
    fn = sprintf('pp_gla_src_pupil_power_correlations_s%d_b%d_v%d',isubj,iblock,v);
    if tp_parallel(fn,outdir,1,0)
      continue
    end
    %
    fprintf('Processing subj%d block%d ...\n',isubj,iblock);
    
    try
      % load pupil data
      if isubj < 10
        load(sprintf('~/pp/data_gla/fw4bt/osfstorage/data/gla01/pupil/sub0%d_gla_pp.mat',isubj));
      else
        load(sprintf('~/pp/data_gla/fw4bt/osfstorage/data/gla01/pupil/sub%d_gla_pp.mat',isubj));
      end
      
      pupil = data.trial{1}';
      f_sample = data.fsample;
      
      % load meg data
      if isubj < 10
        load(sprintf('~/pp/data_gla/fw4bt/osfstorage/data/gla01/meg/sub0%d_gla_meg.mat',isubj));
      else
        load(sprintf('~/pp/data_gla/fw4bt/osfstorage/data/gla01/meg/sub%d_gla_meg.mat',isubj));
      end
      
      
      artifPnts=data.cfg.artfctdef.visual.artifact;
      
      cfg=[];
      cfg.layout='4D248.lay';
      lay = ft_prepare_layout(cfg);
      [~, outp.chanidx] = ismember(lay.label(1:248),data.label);
      
      chanidx = outp.chanidx;
      save(sprintf('~/pp/proc/src/chanidx_s%d.mat',isubj),'chanidx')
      
%       
      [outp.pxx,outp.fxx]=pwelch(data.trial{1}',hanning(400),0,1:1:200,400);
            
      outp.sens_pow    = nan(248,25);
      outp.sens_r      = nan(248,25);
      outp.sens_r_sp   = nan(248,25);
      outp.sens_mi     = nan(248,25);
      outp.sens_mi_df  = nan(248,25);
    catch me
      src_r = nan(246,25);
      save([outdir fn '.mat'],'src_r')
      continue
    end
    
    k = 2;
    fnq = f_sample/2;
    hil_hi = 0.005;
    hil_lo = 2;
    hil_Wn=[hil_hi/fnq hil_lo/fnq];
    [bhil, ahil] = butter(k, hil_Wn);
    
    pupil = filtfilt(bhil, ahil, pupil(:,4));
    
    
    if lag
        pup_shift = round(f_sample*0.93); % 930s from hoeks and levelt (1992?)
        pupil = pupil(pup_shift:end); pupil(end+1:end+pup_shift-1)=nan;
    end
    
    pupil_df = diff(pupil);
    
    data.trial{1}(:,isnan(pupil))=nan(size(data.trial{1},1),sum(isnan(pupil)));
    data.avg = data.trial{1}'; %data.trial{1} = [];

    if isubj < 10
      load(sprintf('~/pp/data_gla/fw4bt/osfstorage/data/gla01/leadfields/sub0%d_gla_lf_BNA5mm.mat',isubj))
    else
      load(sprintf('~/pp/data_gla/fw4bt/osfstorage/data/gla01/leadfields/sub%d_gla_lf_BNA5mm.mat',isubj))
    end
    
    for iart = 1 : size(artifPnts,1)
        data.avg(artifPnts(iart,1):artifPnts(iart,2),:)=NaN;
    end
      
    for ifreq=1:size(freqoi,2)
      
      fprintf('Freq: %d\n',ifreq)
      
      % -------------------------------
      % compute csd
      % -------------------------------
      para          = [];
      para.freq     = freqoi(ifreq);
      para.fsample  = 400;  
      para.overlap = 0.5;
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
      [filt,pow]      = tp_beamformer(real(csd),lf,para);
      % -------------------------------
      % Project noise
      % -------------------------------
      Lr = reshape(lf,[size(lf,1) 8799*3]); csd_noise = Lr*Lr';
      [~,noise]      = tp_beamformer(real(csd_noise),lf,para);
      outp.src_nai(:,ifreq) = pow./noise;
      % -------------------------------
      idx = find(~isnan(dataf(1,:))'&~isnan(pup)&~isnan(pup_df));
      % -------------------------------
      % sensor-level pupil power correlation
      % -------------------------------
      env = (abs(dataf(outp.chanidx(outp.chanidx>0),idx)).^2)';
      
      outp.sens_r(outp.chanidx>0,ifreq) = corr(pup(idx),env,'type','spearman');
      outp.sens_r_sp(outp.chanidx>0,ifreq) = corr(pup(idx),env,'type','spearman');
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
      nlags=floor(10/(opt.n_shift/400)); % roughly 10s
      for isens = 1 : size(env,2)
        [outp.xcorr{ifreq}(:,isens),lags]=xcorr(pup(idx),env(:,isens),nlags,'normalized');
      end
      lags=lags*(opt.n_shift/f_sample);
      outp.xcorr_lags{ifreq} = lags;
      
      % compute source-level power
      % -------------------------------
      src_pow = abs(filt'*dataf).^2;
%       % -------------------------------
      % correlate power with pupil
      % -------------------------------
      outp.src_r(:,ifreq) = corr(pup(idx),src_pow(:,idx)','type','spearman');
      outp.src_r_df(:,ifreq) = corr(pup_df(idx),src_pow(:,idx)','type','spearman');
      cnpup     = copnorm(pup(idx));
      cnpup_df  = copnorm(pup_df(idx));
      cnpow     = copnorm(src_pow(:,idx))';
      for isrc = 1 : size(src_pow,1)
        outp.src_mi(isrc,ifreq) = mi_gg_dfi_ak(cnpow(:,isrc),cnpup,[]);
        outp.src_mi_df(isrc,ifreq) = mi_gg_dfi_ak(cnpow(:,isrc),cnpup_df,[]);
      end
      % -------------------------------

      clear src_pow pup pup_df csd dataf 
      
    end

    save([outdir fn '.mat'],'outp')
    tp_parallel(fn,outdir,0)
    
    clear outp
    
  end
end

error('!')
%%
