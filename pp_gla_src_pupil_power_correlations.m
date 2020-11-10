%% pp_gla_src_pupil_power_correlations
% correlate bandpass filtered (or via wavelets) pupil and MEG signals

clear
restoredefaultpath

% -------------------------
% VERSION 1: no pupil lag
% -------------------------
% v = 1;
% SUBJLIST  = 1:24;
% freqoi    = 2.^(1:(1/4):7);
% lag = 0;
% -------------------------
% VERSION 3: with pupil lag
% -------------------------
v = 2;
SUBJLIST  = 1:24;
freqoi    = 2.^(1:(1/4):7);
lag = 1;
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
for isubj = 1:24
  
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
            
      outp.sens_pow = nan(248,25);
      outp.sens_r      = nan(248,25);
      outp.sens_r_sp   = nan(248,25);
    catch me
      src_r = nan(246,25);
      save([outdir fn '.mat'],'src_r')
      continue
    end
%     len = min([size(data.trial{1},2) size(pupil,1)]);
%     data.trial{1} = data.trial{1}(:,1:len);
%     pupil = pupil(1:len,:);
    
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
   
                
    for ifreq=1:numel(freqoi)
      
      fprintf('Freq: %d\n',ifreq)
    
      % -------------------------------
      % START CK TF-ANALYSIS
      % -------------------------------
%       srate = 400;
%       % freq analysis, pass 1, yields TF Fourier rep ##########################
%       % - this is for multiplication with the ortho-normalised spatial filter
%       tfcfg=[]; % this will be carried along ...
%       tfcfg.method='wavelet';
%       tfcfg.output='fourier';
%       tfcfg.channel={'MEG'};
%       tfcfg.foi=freqoi(ifreq);
%       tfcfg.width=5.83; % again, as per Hipp et al. (2012) Nat Neurosci
%       tempSD=1./(2*pi*(tfcfg.foi./tfcfg.width)); % temporal SD in sec
%       tempWd=round(3*tempSD*srate)/srate; % set to 3 std dev, comp with 1000 Hz sampl rate
%       tfcfg.toi=tempWd.*(1:floor(data.time{1}(end)./tempWd));
%       % keep in mind that fieldtrip uses a proprietary setting of the gwidth
%       % parameter (default = 3 cycles) internally that is independent of the
%       % here requested time axis
%       tfcfg.pad='nextpow2';
%       tf=ft_freqanalysis(tfcfg,data);
%       % reset freq to requested freq
%       tf.freq=freqoi(ifreq);      
      
%       artifPnts=data.cfg.artfctdef.visual.artifact;
      
      for iart = 1 : size(artifPnts,1)
        data.avg(artifPnts(iart,1):artifPnts(iart,2),:)=NaN;
      end
      
%       tfPnts   =tf.time*srate;
% 
%       % only discard those bins that "center" on an artifact
%       critDist=diff(tfPnts([1,2]),[],2); % set critical dist to central 100%
%       keepBins=logical([]);
%       % discard TF bins that overlap with artifacts (set zero in keepBins)
%       for ibin=1:numel(tfPnts)
%           keepBins(ibin)=~any(abs(tfPnts(ibin)-ceil(mean(artifPnts,2)))<critDist);   
%       end

      % additionally omit edge datapoints to excl artif
%       timeAx=tf.time; % original time axis
%       timeKp=dsearchn(timeAx.',[3;timeAx(end)-3]); % excl ~ 1st & last 3 sec
%       keepBins(1:timeKp(1)-1)=false;
%       keepBins(timeKp(2)+1:end)=false;

      % compute cross-specral density by time bin, then average csd's
%       csd=zeros(numel(tf.label)*[1,1]);
%       % csd calc excludes artifact bins
%       csdTime=tf.time(keepBins);
%       csdData=tf.fourierspctrm(:,:,:,keepBins);
%       for itbin=1:numel(csdTime)
%           fspec=squeeze(csdData(:,:,:,itbin)).';
%           for ichan=1:numel(tf.label);
%               csd(:,ichan)=csd(:,ichan)+fspec(ichan)*conj(fspec);
%           end
%       end
%       csd=csd./numel(tf.time(keepBins)); % avg cross-spectral dens matrix
%       csdData=[]; csdTime=[];
      % -------------------------------
      % END CK TF-ANALYSIS
      % -------------------------------
      
      % -------------------------------
      % pupil ck
      % -------------------------------
%       ck_pup = pupil(round(tfPnts));
      
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
%       [ck_filt,ck_pow]      = tp_beamformer(real(csd),lf,para);
      % -------------------------------
      % Project noise
      % -------------------------------
      Lr = reshape(lf,[size(lf,1) 8799*3]); csd_noise = Lr*Lr';
      [~,noise]      = tp_beamformer(real(csd_noise),lf,para);
      outp.src_nai(:,ifreq) = pow./noise;
      % -------------------------------

%       ck_dataf = squeeze(tf.fourierspctrm);
%       
      idx = find(~isnan(dataf(1,:))'&~isnan(pup));
%       ck_idx_valid = find(~isnan(ck_dataf(1,:))'&~isnan(ck_pup)&keepBins');
%       idx=intersect(ck_idx_valid,tp_idx_valid);
      % -------------------------------
      % sensor-level pupil power correlation
      % -------------------------------
      outp.sens_r(outp.chanidx>0,ifreq) = corr(pup(idx),abs(dataf(outp.chanidx(outp.chanidx>0),idx).^2)','type','spearman');
      outp.sens_r_sp(outp.chanidx>0,ifreq) = corr(pup(idx),abs(dataf(outp.chanidx(outp.chanidx>0),idx).^2)','type','spearman');
      % -------------------------------
      % compute source-level power
      % -------------------------------
      tp_src = abs(filt'*dataf).^2;
%       % -------------------------------
      % correlate power with pupil
      % -------------------------------
      outp.src_r(:,ifreq) = corr(pup(idx),tp_src(:,idx)','type','spearman');
      outp.src_r_df(:,ifreq) = corr(pup_df(idx),tp_src(:,idx)','type','spearman');
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
