%% pp_gla_src_pupil_power_correlations
% correlate bandpass filtered (or via wavelets) pupil and MEG signals

clear
restoredefaultpath

% -------------------------
% VERSION 4
% -------------------------
v = 3;
% include 28 subjects, as in pfeffer et al. (2018) plos biology
SUBJLIST  = 1:24;
freqoi    = 2.^(1:(1/4):7);
% -------------------------

addpath ~/pconn/matlab/
load(sprintf('~/pp/proc/pp_atlas_BNA.mat'))


outdir = '~/pp/proc/src/';
ord    = pconn_randomization;

addpath /home/gnolte/meth/highlevel/
%%
% -------------------------
for isubj = 1:24
  addpath('~/Documents/MATLAB/fieldtrip-20181231/')
  ft_defaults
  clear data dat pupil pup dataf src_r
  
  for iblock = 1:1
%     %
%     fn = sprintf('pp_gla_src_pupil_power_correlations_s%d_b%d_v%d',isubj,iblock,v);
%     if tp_parallel(fn,outdir,1,0)
%       continue
%     end
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
      
      rmpath(genpath('~/Documents/MATLAB/fieldtrip-20181231/'))
      
      [outp.pxx,outp.fxx]=pwelch(data.trial{1}',hanning(400),0,1:1:200,400);
            
      outp.tp_sens_pow = nan(248,25);
      outp.sens_r      = nan(248,25);
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
    
%     pup_shift = round(f_sample*0.93); % 930s from hoeks and levelt (1992?)
%     pupil = pupil(pup_shift:end); pupil(end+1:end+pup_shift-1)=nan;
%     pupil_df = diff(pupil);
    
    pupil = pupil(500:end-500);
    len = length(pupil);

    data.avg = filtfilt(bhil, ahil, data.trial{1}(:,500:end-500)'); %data.trial{1} = [];
    
    data.avg = data.avg(1:len,:);
    
    outp.sens_r = corr(pupil,data.avg);
    
    for isens = 1 : size(data.avg,2)
      [outp.xcorr(:,isens),outp.lag]=xcorr(pupil,data.avg(:,isens),1200,'normalized');
    end

    if isubj < 10
      load(sprintf('~/pp/data_gla/fw4bt/osfstorage/data/gla01/leadfields/sub0%d_gla_lf_BNA5mm.mat',isubj))
    else
      load(sprintf('~/pp/data_gla/fw4bt/osfstorage/data/gla01/leadfields/sub%d_gla_lf_BNA5mm.mat',isubj))
    end
         
    for ifreq=1:numel(freqoi)
      
      fprintf('Freq: %d\n',ifreq)
      
      for iart = 1 : size(artifPnts,1)
        data.avg(artifPnts(iart,1):artifPnts(iart,2),:)=NaN;
      end
 
      % -------------------------------
      % compute csd
      % -------------------------------
      para          = [];
      para.freq     = freqoi(ifreq);
      para.fsample  = 400;  
      para.overlap  = 0.5;
      [csd, dataf,opt]=tp_compute_csd_wavelets(data.avg',para);
      
      % -------------------------------
      % compute power from csd
      % -------------------------------
      tmp = diag(abs(csd)); 
      outp.tp_sens_pow(outp.chanidx>0,ifreq) = tmp(outp.chanidx(outp.chanidx>0));

      % -------------------------------
      % prepare pupil signal
      % -------------------------------
      nseg=floor((size(pupil_df,1)-opt.n_win)/opt.n_shift+1);
      
      pup = nan(nseg,1);
      pup_df= nan(nseg,1);
      for j=1:nseg
        tmp = pupil((j-1)*opt.n_shift+1:(j-1)*opt.n_shift+opt.n_win);
        tmp2 = pupil_df((j-1)*opt.n_shift+1:(j-1)*opt.n_shift+opt.n_win);
        pup(j) = mean(tmp.*gausswin(opt.n_win,3));
        pup_df(j) = mean(tmp2.*gausswin(opt.n_win,3));
      end
      
      % find indices of non-NAN segments
      idx_valid = find(~isnan(dataf(1,:)));
      idx_valid_df = find(~isnan(pup_df)');
      idx_valid = intersect(idx_valid,idx_valid_df);
      
      env = abs(dataf(:,idx_valid)).^2;
      f_sample = 400/opt.n_shift;
      env_filt=lowpass(env,hil_lo,f_sample);
      
      % correlate pupil with sensor level signal
      outp.sens_r(outp.chanidx>0,ifreq) = corr(pup(idx_valid),env(outp.chanidx(outp.chanidx>0),:)','type','Spearman');
      outp.sens_r_df(outp.chanidx>0,ifreq) = corr(pup_df(idx_valid),env(outp.chanidx(outp.chanidx>0),:)','type','Spearman');
      outp.sens_r_filt(outp.chanidx>0,ifreq) = corr(pup(idx_valid),env_filt(outp.chanidx(outp.chanidx>0),:)','type','Spearman');
      outp.sens_r_df_filt(outp.chanidx>0,ifreq) = corr(pup_df(idx_valid),env_filt(outp.chanidx(outp.chanidx>0),:)','type','Spearman');
      % -------------------------------
      % beamforming
      % -------------------------------
      para      = [];
      para.iscs = 1;
      para.reg  = 0.05;
      [filt,tp_pow]      = tp_beamformer(real(csd),lf,para);
      % -------------------------------
      % Project noise
      % -------------------------------
      Lr = reshape(lf,[size(lf,1) 8799*3]); csd_noise = Lr*Lr';
      [~,noise]      = tp_beamformer(real(csd_noise),lf,para);
      outp.tp_src_nai(:,ifreq) = tp_pow./noise;
      % -------------------------------

      % -------------------------------
      % compute source-level power
      % -------------------------------
      env = abs(filt'*dataf(:,idx_valid)).^2;
      f_sample = 400/opt.n_shift;
      env_filt=lowpass(env,hil_lo,f_sample);
      
%     % -------------------------------
      % correlate power with pupil
      % -------------------------------
      outp.src_r(:,ifreq) = corr(pup(idx_valid),env','type','Spearman');
      outp.src_r_df(:,ifreq) = corr(pup_df(idx_valid),env','type','Spearman');
      outp.src_r_filt(:,ifreq) = corr(pup(idx_valid),env_filt','type','Spearman');
      outp.src_r_df_filt(:,ifreq) = corr(pup_df(idx_valid),env_filt','type','Spearman');
      % -------------------------------
      
      clear src pup csd tp_csd dataf
      
    end

    save([outdir fn '.mat'],'outp')
    tp_parallel(fn,outdir,0)
    
    clear outp
    
  end
end

error('!')
%%
