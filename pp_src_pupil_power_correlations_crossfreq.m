%% pp_src_pupil_power_correlations_crossfreq
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

pupil_freqoi(:,1) = 2.^(-9:(0.5):1);
pupil_freqoi(:,2) = 2.^(-8:(0.5):2);
%%
% -------------------------
for isubj = SUBJLIST
  
  im = find(ord(isubj,:)==1);
  
  for iblock = 1:2
%     
    fn = sprintf('pp_src_pupil_power_correlations_crossfreq_s%d_b%d_v%d',isubj,iblock,v);
    if tp_parallel(fn,outdir,1,0)
      continue
    end
% %     
    fprintf('Processing subj%d block%d ...\n',isubj,iblock);
    
    try
      % load meg data
      load(sprintf('~/pupmod/proc/sens/pupmod_rest_sens_cleandat_s%d_m%d_b%d_v%d.mat',isubj,im,iblock,1))
      % load pupil data
      load(sprintf('~/pp/proc/pup/pp_pupil_diameter_cleaned_s%d_m%d_b%d.mat',isubj,im,iblock))
      
    catch me
      outp.src_r = nan(246,25);
      save([outdir fn '.mat'],'outp')
      continue
    end
    
    if size(pupil,2)>3
      raw_pupil = pupil(:,4);
    else
      raw_pupil = pupil(:,1);
    end
      
    clear pupil
    
    raw_pupil = resample(raw_pupil,400,1000);

    raw_pupil = raw_pupil(end:-1:1);

    len = min([size(raw_pupil,1) size(dat,2)]);
    if len/400 > 600
      len = 400*600;
    end
    
    raw_pupil = raw_pupil(1:len);
    raw_pupil = raw_pupil(end:-1:1);
    
    pup_shift = round(400*0.93); % 930s from hoeks and levelt (1992?)
    raw_pupil = raw_pupil(pup_shift:end);
    raw_pupil(end+1:end+pup_shift-1)=nan;
    
    dat = dat(:,1:len);
    dat(:,isnan(raw_pupil))=[];
    raw_pupil(isnan(raw_pupil))=[];
  
    for ipupfreq = 1 : size(pupil_freqoi,1) 
      pupil(:,ipupfreq) = ft_preproc_bandpassfilter(raw_pupil', 400, [pupil_freqoi(ipupfreq,1) pupil_freqoi(ipupfreq,2)], 1, 'but', [], 'no');
    end
        
    load(['/home/tpfeffer/pp/proc/src/' sprintf('pp_sa_s%d_m%d_b%d_v%d.mat',isubj,im,iblock,1)],'sa');
        
    for ifreq=1:numel(freqoi)
      
      para             = [];
      para.freq        = freqoi(ifreq);
      para.fsample     = 400;  
      [csd, dataf,opt] = tp_compute_csd_wavelets(dat,para);
                  
       % -------------------------------
      % prepare pupil signal
      % -------------------------------
      nseg=floor((size(dat,2)-opt.n_win)/opt.n_shift+1);
      
      for ipupfreq = 1 : size(pupil_freqoi,1)
        ipupfreq
        for j=1:nseg
          dloc2=dat(:,(j-1)*opt.n_shift+1:(j-1)*opt.n_shift+opt.n_win)';
          if any(isnan(dloc2(:,1)))
            pup{ipupfreq}{ifreq}(:,j) = nan;
            continue
          else
            tmp = pupil((j-1)*opt.n_shift+1:(j-1)*opt.n_shift+opt.n_win,ipupfreq);
            pup{ipupfreq}{ifreq}(:,j) = mean(tmp.*gausswin(opt.n_win,3));                 
          end
        end
      end
      
      % Compute cross spectrum & find NaN segments
      idx_valid = find(~isnan(dataf(1,:)));
      nanidx{ifreq} = idx_valid;
      % -----------
      % beamforming (ignore name of leadfield)
      % --------------
      para     	= [];
      para.reg  = 0.05;
      [filt]    = tp_beamformer(real(csd),sa.L_genemaps_aal,para);
      % --------------
      % compute power
      % --------------
      src_pow{ifreq} = abs(filt'*dataf).^2;
     
    end
    
    % compute pupil diameter for segments
    % ---------------
    src_r = nan(size(sa.L_genemaps_aal,2),size(src_pow,2),size(pupil,2));
    
    for ipup = 1 : size(pupil,2)
      ipup
      for ifreq = 1 : size(src_pow,2)
        src_r(:,ifreq,ipup) = corr(src_pow{ifreq}(:,nanidx{ifreq})',pup{ipup}{ifreq}(nanidx{ifreq})');
      end
    end
    
    save([outdir fn '.mat'],'src_r')
    tp_parallel(fn,outdir,0)
    
    
    clear src_r all_nai outp
    
  end
end

error('!')

