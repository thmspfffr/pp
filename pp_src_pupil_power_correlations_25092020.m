%% pp_src_pupil_power_correlations
% correlate bandpass filtered (or via wavelets) pupil and MEG signals

clear
restoredefaultpath

% -------------------------
% VERSION 4
% -------------------------
v = 3;
% include 28 subjects, as in pfeffer et al. (2018) plos biology
SUBJLIST = [4 5 6 7 8 9 10 11 12 13 15 16 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34];
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
    if tp_parallel(fn,outdir,1,0)
      continue
    end
% %     
    fprintf('Processing subj%d block%d ...\n',isubj,iblock);
    
    try
      % load cleaned meg data
      load(sprintf('~/pupmod/proc/sens/pupmod_rest_sens_cleandat_s%d_m%d_b%d_v%d.mat',isubj,im,iblock,1))
      % load cleanted pupil data
      load(sprintf('~/pp/proc/pup/pp_pupil_diameter_cleaned_s%d_m%d_b%d.mat',isubj,im,iblock))
    catch me
      continue
    end
    
    % bp-filter and resample pupil
    % ------
    k = 2; f_sample = 1000;
    fnq = f_sample/2;
    hil_hi = 0.005; hil_lo = 2;
    hil_Wn=[hil_hi/fnq hil_lo/fnq];
    [bhil, ahil] = butter(k, hil_Wn);
 
    pupil = filtfilt(bhil, ahil, pupil(:,4));
    pupil = resample(pupil,400,1000);
    % ------
    
    % align pupil and meg (at signal offset)
    % ------
    if size(pupil,2)>3
      pupil = pupil(end:-1:1,4);
    else
      pupil = pupil(end:-1:1);
    end
    
%     data.trial = data.trial(:,1:data.end_of_recording);
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
    pup_shift = round(400*0.93);   
    pupil = pupil(pup_shift:end); pupil(end+1:end+pup_shift-1)=nan;
    
    pupil_df = diff(pupil);
      
    dat(:,isnan(pupil))=nan(size(dat,1),sum(isnan(pupil)));
    
    load(['/home/tpfeffer/pp/proc/src/' sprintf('pp_sa_s%d_m%d_b%d_v%d.mat',isubj,im,iblock,1)],'sa');  
    
    lf = sa.L_genemaps_aal;
    
    save(sprintf('~/pp/proc/sens/pp_meg_pupil_lfs_s%d_b%d.mat',isubj,iblock),'dat','pupil','lf');
    
    [outp.pxx,outp.fxx]=pwelch(dat(:,~isnan(dat(1,:)))',hanning(400),0,1:1:200,400);
    
    for ifreq=1:numel(freqoi)
      ifreq

      para          = [];
      para.freq     = freqoi(ifreq);
      para.fsample  = 400;  
      [csd, dataf,opt]=tp_compute_csd_wavelets(dat,para);
            
      outp.tp_sens_pow(:,ifreq) = diag(abs(csd));

      % -------------------------------
      % prepare pupil signal
      % -------------------------------
      nseg=floor((size(dat,2)-opt.n_win)/opt.n_shift+1);
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
      
      env = abs(dataf(:,idx_valid).^2);
      
      k = 2; f_sample = 400/opt.n_shift;
      fnq = f_sample/2;
      hil_hi = 0.005; hil_lo = 2;
      hil_Wn=[hil_hi/fnq hil_lo/fnq];
      [bhil, ahil] = butter(k, hil_Wn);

      env = filtfilt(bhil, ahil, env);
 
      % correlate pupil with sensor level signal
      outp.sens_r(:,ifreq) = corr(pup(idx_valid),abs(dataf(:,idx_valid).^2)');
      outp.sens_r_df(:,ifreq) = corr(pup_df(idx_valid),abs(dataf(:,idx_valid).^2)');
      
      % beamforming (ignore name of leadfield)
      % --------------
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
      % compute power
      src = abs(filt'*dataf).^2;
      
      % correlate with pupila
      outp.src_r(:,ifreq) = corr(pup(idx_valid),src(:,idx_valid)','type','Spearman');
      outp.src_r_df(:,ifreq) = corr(pup_df(idx_valid),src(:,idx_valid)','type','Spearman');
      
      clear src pup

    end

    save([outdir fn '.mat'],'outp')
    tp_parallel(fn,outdir,0)
    
    clear src_r all_nai outp
    
  end
end

error('!')

