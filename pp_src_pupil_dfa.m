%% pp_src_pupil_dfa
% correlate bandpass filtered (or via wavelets) pupil and MEG signals

clear
restoredefaultpath

SUBJLIST  = [4 5 6 7 8 9 10 11 12 13 15 16 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34];
v = 1;
addpath ~/Documents/MATLAB/fieldtrip-20160919/
ft_defaults

outdir = '~/pp/proc/dfa/';
addpath ~/pconn/matlab/
ord = pconn_randomization;
load(sprintf('~/pp/proc/pp_atlas_BNA.mat'))

calc_interv = [3 75];
overlap = 0.9;

for isubj = SUBJLIST
  
  im = find(ord(isubj,:)==1);
  
  for iblock = 1:2
    
    fn = sprintf('pp_src_dfa_s%d_b%d_v%d',isubj,iblock,v);
    if tp_parallel(fn,outdir,1,0)
      continue
    end
    
    fprintf('Processing subj%d block%d ...\n',isubj,iblock);
    
    try
      % load meg data
      load(sprintf('~/pp/proc/pp_sens_cleandat_s%d_m%d_b%d_v%d.mat',isubj,im,iblock,v))
      % load pupil data
      load(sprintf('~/pp/proc/pup/pp_pupil_diameter_cleaned_s%d_m%d_b%d.mat',isubj,im,iblock))
      
    catch me
      seg.dfa = nan;
      seg.pupil = nan;
      save([outdir fn '.mat'],'seg')
      continue
    end
    
    k = 2;
    f_sample = 1000;
    fnq = f_sample/2;
    hil_hi = 0.005;
    hil_lo = 2;
    hil_Wn=[hil_hi/fnq hil_lo/fnq];
    [bhil, ahil] = butter(k, hil_Wn);
    
    pupil = filtfilt(bhil, ahil, pupil(:,4));
    
    pupil = resample(pupil,400,1000);
    
    data.trial = data.trial(:,1:data.end_of_recording);
    data.trial = data.trial(:,end:-1:1);
    pupil = pupil(end:-1:1);
    
    len = min([size(pupil,1) size(data.trial,2)]);
    if len/400 > 600
      len = 400*600;
    end
    
    data.trial = data.trial(:,1:len);
    pupil = pupil(1:len);
    
    data.trial = data.trial(:,end:-1:1);
    pupil = pupil(end:-1:1);
    
    % replace this later
    pupil(isnan(data.trial(1,:)))=[];
    data.trial(:,isnan(data.trial(1,:))) = [];
    
    flp = 8;           % lowpass frequency of filter
    fhi = 12;
    
    para.ord = 4;
    k=4;                  % 2nd order butterworth filter
    fnq=400/2;       % Nyquist frequency
    Wn=[flp/fnq fhi/fnq]; % butterworth bandpass non-dimensional frequency
    [bfilt,afilt]=butter(k,Wn);
    
    data_f = abs(hilbert(single(filtfilt(bfilt,afilt,data.trial'))));
    
    load(['/home/tpfeffer/pp/proc/src/' sprintf('pp_sa_s%d_m%d_b%d_v%d.mat',isubj,im,iblock,v)],'sa');
    
    eplen    = size(data_f,1);
    seglen   = floor(90*400);
    segshift = floor(overlap*seglen);
    nseg     = floor((eplen-seglen)/segshift+1);
    
    
    if isubj >= 32
      load(sprintf('~/pconn/proc/preproc/pconn_chanidx_s%d_m%d.mat',31,3))
    elseif isubj < 4
      load(sprintf('~/pconn/proc/preproc/pconn_chanidx_s%d_m%d.mat',4,1))
    elseif isubj == 17
      load(sprintf('~/pconn/proc/preproc/pconn_chanidx_s%d_m%d.mat',16,1))
    else
      load(sprintf('~/pconn/proc/preproc/pconn_chanidx_s%d_m%d.mat',isubj,im))
    end
    clear seg
    for itime = 1 : nseg
      
      timerange = [(itime-1)*segshift+1 (itime-1)*segshift+seglen];
      
      pup_dat          = pupil((itime-1)*segshift+1:(itime-1)*segshift+seglen);
      dat              = data_f((itime-1)*segshift+1:(itime-1)*segshift+seglen,:);
      
      seg.pupil(itime) = mean(pup_dat);
      tmp = tp_dfa(dat,calc_interv,400,0.5,15);
      
      seg.dfa(:,itime) = pconn_sens_interp274(idx,tmp.exp);
      
    end
    save([outdir fn '.mat'],'seg')
    tp_parallel(fn,outdir,0)
    
  end
end

error('!')
%%
%
clear dfa
i = 0
for isubj = SUBJLIST(1:end-1)
  for iblock = 1 : 2
    clear src_r
    try
      load(sprintf([outdir 'pp_src_dfa_s%d_b%d_v%d.mat'],isubj,iblock,v));
      isubj
      
      if ~isnan(seg.dfa)
       dfa.hi_dfa(:,isubj,iblock) = mean(seg.dfa(:,seg.pupil>median(seg.pupil)),2);
       dfa.lo_dfa(:,isubj,iblock) = mean(seg.dfa(:,seg.pupil<median(seg.pupil)),2);
       dfa.hi_pup(isubj,iblock) = mean(seg.pupil(seg.pupil>median(seg.pupil)));
       dfa.lo_pup(isubj,iblock) = mean(seg.pupil(seg.pupil<median(seg.pupil)));
      else
       dfa.hi_dfa(:,isubj,iblock) = nan(274,1);
       dfa.lo_dfa(:,isubj,iblock) = nan(274,1);
       dfa.hi_pup(isubj,iblock) = nan;
       dfa.lo_pup(isubj,iblock) = nan;
      end
    catch me
       dfa.hi_dfa(:,isubj,iblock) = nan(274,1);
       dfa.lo_dfa(:,isubj,iblock) = nan(274,1);
       dfa.hi_pup(isubj,iblock) = nan;
       dfa.lo_pup(isubj,iblock) = nan;
      continue
    end
    %     all_corr(:,:,isubj,iblock) = r;
    
  end
end

dfa.hi_dfa = nanmean(dfa.hi_dfa(:,SUBJLIST(1:end-1),:),3);
dfa.lo_dfa = nanmean(dfa.lo_dfa(:,SUBJLIST(1:end-1),:),3);
% dfa.hi_pup(isubj,iblock) = nan;
% dfa.lo_pup(isubj,iblock) = nan;
%% PLOT

load /home/tpfeffer/pconn/proc/src/pconn_sa_s26_m2_b1_v4.mat
pars            = [];
pars.markersize = 0;
pars.linewidth  = 9;
pars.cbar       = 0;
pars.scale      = [-0.02 0.02]
pars.cmap       = jet;
pars.resolution = 600;

par = mean(dfa.hi_dfa,2)-mean(dfa.lo_dfa,2)

figure; set (gcf,'color','w')

showfield_colormap(par,sa.locs_2D,pars);
drawnow

