
%% pp_hh_src_powerspectra
% correlate bandpass filtered (or via wavelets) pupil and MEG signals

clear
restoredefaultpath

% -------------------------
% VERSION 1: no pupil lag
% -------------------------
% v = 1;
% SUBJLIST = [4 5 6 7 8 9 10 11 12 13 15 16 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34];
% lag = 0;
% win_len = 800;
% overlap = 2; % 50% overlap
% -------------------------
% VERSION 2: with pupil lag
% -------------------------
% v = 2;
% SUBJLIST = [4 5 6 7 8 9 10 11 12 13 15 16 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34];
% lag = 1;
% win_len = 800;
% overlap = 2; % 50% overlap
% -------------------------
% VERSION 11: no pupil lag, less overlap
% -------------------------
% v = 11;
% SUBJLIST = [4 5 6 7 8 9 10 11 12 13 15 16 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34];
% lag = 0;
% win_len = 800;
% overlap = 1; % 0% overlap
% -------------------------
% VERSION 2: with pupil lag
% -------------------------
v = 22;
SUBJLIST = [4 5 6 7 8 9 10 11 12 13 15 16 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34];
lag = 1;
win_len = 800;
overlap = 1; % 0% overlap
% -------------------------

addpath ~/Documents/MATLAB/fieldtrip-20160919/
addpath ~/pconn/matlab/
load(sprintf('~/pp/proc/pp_atlas_BNA.mat'))

% ft_defaults

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
    fn = sprintf('pp_hh_src_powerspectra_s%d_b%d_v%d',isubj,iblock,v);
    if tp_parallel(fn,outdir,1,0)
      continue
    end
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
    
    clear csd
    
    for ifreq=1:length(freqoi)
      ifreq
      
      para          = [];
      para.freq     = freqoi(ifreq);
      para.fsample  = 400;
      para.overlap  = 0.5;
      [csd(:,:,ifreq)]=tp_compute_csd_wavelets(dat,para);
      
    end
    
    csd = nanmean(csd,3);
    
    para          = [];
    para.reg      = 0.05;
    filt = tp_beamformer(real(csd),sa.L_genemaps_aal,para);
    % --------------
    
    idx = isnan(dat(1,:))' | isnan(pupil);
    
    dat(:,idx) = nan;
    pupil(idx) = nan;
    pupil_df(idx) = nan;
    
    opt.n_win = win_len; % 10s segment length, i.e., 0.1:0.1:100
    opt.n_shift = win_len/overlap; 
    
    nseg=floor((size(dat,2)-opt.n_win)/opt.n_shift+1);
    
    clear pxx fxx pup pup_df
    ff = 2:1/(opt.n_win/400):128;
    
    pxx = nan(size(ff,2),max(BNA.tissue_5mm(:)),nseg);
    for iseg = 1 : nseg
      
      fprintf('%d / %d\n',iseg,nseg)
      seg_dat = dat(:,(iseg-1)*opt.n_shift+1:(iseg-1)*opt.n_shift+opt.n_win)'*filt;
      
      if any(isnan(seg_dat(:,1)))
        pup(iseg) = nan;
        pup_df(iseg)=nan;
        continue
      end
      
      [tmp_pxx,fxx]=pwelch(seg_dat,hanning(opt.n_win),0.5,ff,400,'power');
      
      for igrid = 1 : max(BNA.tissue_5mm(:))
        pxx(:,igrid,iseg) = mean(tmp_pxx(:,BNA.tissue_5mm == igrid),2);
      end
      
      pup(iseg)  = nanmean(pupil((iseg-1)*opt.n_shift+1:(iseg-1)*opt.n_shift+opt.n_win));
      if iseg~=nseg
        pup_df(iseg) = nanmean(pupil_df((iseg-1)*opt.n_shift+1:(iseg-1)*opt.n_shift+opt.n_win));
      else
        pup_df(iseg) = nanmean(pupil_df((iseg-1)*opt.n_shift+1:(iseg-1)*opt.n_shift+opt.n_win-1));
      end
    end
    
    pxx=single(pxx);
    save([outdir fn '.mat'],'pxx','fxx','pup','pup_df')
    tp_parallel(fn,outdir,0)
    
    clear pxx fxx pup pup_df
    
  end
end

% exit






