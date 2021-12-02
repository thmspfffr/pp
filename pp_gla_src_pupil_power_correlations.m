%% pp_gla_src_pupil_power_correlations
% correlate bandpass filtered (or via wavelets) pupil and MEG signals

clear
restoredefaultpath

% -------------------------
% VERSION 1: no pupil lag (0 ms)
% -------------------------
% v = 1;
% lag = 0;
% -------------------------
% VERSION 2: with pupil lag (930 ms)
% -------------------------
% v = 2;
% lag = 1;
% -------------------------
% VERSION 3: with pupil lag (500 ms)
% -------------------------
v = 3;
lag = 2;
% -------------------------

addpath('~/Documents/MATLAB/fieldtrip-20181231/')
addpath ~/pconn/matlab/
load(sprintf('~/pp/proc/pp_atlas_BNA.mat'))
freqoi    = 2.^(1:(1/4):7);
SUBJLIST  = 1:24;

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
      
      % load meg data
      if isubj < 10
        load(sprintf('~/pp/data_gla/fw4bt/osfstorage/data/gla01/meg/sub0%d_gla_meg.mat',isubj));
      else
        load(sprintf('~/pp/data_gla/fw4bt/osfstorage/data/gla01/meg/sub%d_gla_meg.mat',isubj));
      end
      
      artifPnts=data.cfg.artfctdef.visual.artifact;
      f_sample = data.fsample;
      
      cfg=[];
      cfg.layout='4D248.lay';
      lay = ft_prepare_layout(cfg);
      [~, outp.chanidx] = ismember(lay.label(1:248),data.label);
      
      chanidx = outp.chanidx;
%       save(sprintf('~/pp/proc/src/chanidx_s%d.mat',isubj),'chanidx')
      
%       
      [outp.pxx,outp.fxx]=pwelch(data.trial{1}',hanning(f_sample),0,1:1:200,400);
            
      outp.sens_pow    = nan(248,25);
      outp.sens_r      = nan(248,25);
      outp.sens_r_sp   = nan(248,25);
      outp.sens_mi     = nan(248,25);
      outp.sens_mi_df  = nan(248,25);
      outp.sens_mi0    = nan(248,25);
      outp.sens_mi_df0 = nan(248,25);
      
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
    
    % pupil shift: 930 ms from hoeks & levelt (1992)
    if lag==1 % 930 ms lag
      pup_shift = round(f_sample*0.93);
      pupil = pupil(pup_shift:end); pupil(end+1:end+pup_shift-1)=nan;
    elseif lag == 2 % 500 ms lag
      pup_shift = round(f_sample*0.5);
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
      para              = [];
      para.freq         = freqoi(ifreq);
      para.fsample      = 400;  
      para.overlap      = 0.8;
      [csd, dataf,opt]  = tp_compute_csd_wavelets(data.avg',para);
      
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
      para       = [];
      para.iscs  = 1;
      para.reg   = 0.05;
      [filt,pow] = tp_beamformer(real(csd),lf,para);
      % -------------------------------
      % Project noise
      % -------------------------------
      Lr = reshape(lf,[size(lf,1) 8799*3]); csd_noise = Lr*Lr';
      [~,noise]             = tp_beamformer(real(csd_noise),lf,para);
      outp.src_nai(:,ifreq) = pow./noise;
      % -------------------------------
      idx = find(~isnan(dataf(1,:))'&~isnan(pup)&~isnan(pup_df));
      % -------------------------------
      % sensor-level pupil power correlation
      % -------------------------------
      env = (abs(dataf(outp.chanidx(outp.chanidx>0),idx)).^2)';
      
      outp.sens_r(outp.chanidx>0,ifreq)     = corr(pup(idx),env,'type','spearman');
      outp.sens_r_sp(outp.chanidx>0,ifreq)  = corr(pup(idx),env,'type','spearman');
      % -------------------------------
      % sensor-level mutual information
      % -------------------------------
      cnpup     = copnorm(pup(idx));
      cnpup0    = cnpup(end:-1:1);
      cnpup_df  = copnorm(pup_df(idx));
      cnpup_df0 = cnpup_df(end:-1:1);
      cnpow     = copnorm(env);
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
      outp.xcorr{ifreq} = nan(nlags*2+1,248);
      outp.xcorr_df{ifreq} = nan(nlags*2+1,248);
      tmp_xcorr = []; tmp_xcorr_df = [];
      for isens = 1 : size(env,2)
        tmp_pup = pup(idx)-mean(pup(idx));
        tmp_pup_df = pup_df(idx)-mean(pup_df(idx));
        tmp_env = env(:,isens)-nanmean(env(:,isens),1);
        [tmp_xcorr(:,isens),lags] = xcorr(tmp_pup,tmp_env,nlags,'coeff');
        [tmp_xcorr_df(:,isens),lags] = xcorr(tmp_pup_df,tmp_env,nlags,'coeff');
      end
      outp.xcorr{ifreq}(:,outp.chanidx>0) = tmp_xcorr;
      outp.xcorr_df{ifreq}(:,outp.chanidx>0) = tmp_xcorr_df;
      lags = lags*(opt.n_shift/f_sample);
      outp.xcorr_lags{ifreq} = lags;
      
      % compute source-level power
      % -------------------------------
      src_pow = abs(filt'*dataf(:,idx)).^2;
%       % -------------------------------
      % correlate power with pupil
      % -------------------------------
      outp.src_r(:,ifreq) = corr(pup(idx),src_pow','type','spearman');
      outp.src_r_df(:,ifreq) = corr(pup_df(idx),src_pow','type','spearman');
      cnpup     = copnorm(pup(idx));
      cnpup_df  = copnorm(pup_df(idx));
      cnpow     = copnorm(src_pow)';
      for isrc = 1 : size(src_pow,1)
        outp.src_mi(isrc,ifreq)     = mi_gg_dfi_ak(cnpow(:,isrc),cnpup,[]);
        outp.src_mi_df(isrc,ifreq)  = mi_gg_dfi_ak(cnpow(:,isrc),cnpup_df,[]);
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

      clear src_pow pup pup_df csd dataf 
      
    end

    save([outdir fn '.mat'],'outp')
    tp_parallel(fn,outdir,0)
    
    clear outp
    
  end
end

% exit

