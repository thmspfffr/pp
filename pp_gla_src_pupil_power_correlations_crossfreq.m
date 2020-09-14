%% pp_gla_src_pupil_power_correlations_crossfreq
% correlate bandpass filtered (or via wavelets) pupil and MEG signals

clear
restoredefaultpath

SUBJLIST  = 1:24;
v = 1;

addpath ~/Documents/MATLAB/fieldtrip-20160919/
ft_defaults
pupil_freqoi(:,1) = 2.^(-9:(0.5):1);
pupil_freqoi(:,2) = 2.^(-8:(0.5):2);

freqoi=2.^(1:(1/4):7); % 1-128 Hz as per Hipp et al. (2012) Nat Neurosci

    
outdir = '~/pp/proc/src/';
addpath ~/pconn/matlab/
ord = pconn_randomization;
load(sprintf('~/pp/proc/pp_atlas_BNA.mat'))

%%
for isubj = SUBJLIST
  
  im = find(ord(isubj,:)==1);
  
  for iblock = 1:1
          
    clear src_pow dat pup pupil dataf

    fn = sprintf('pp_gla_src_pupil_power_correlations_crossfreq_s%d_b%d_v%d',isubj,iblock,v);
    if tp_parallel(fn,outdir,1,0)
      continue
    end
    
    fprintf('Processing subj%d block%d ...\n',isubj,iblock);

    try
      % load pupil data
      if isubj < 10
        load(sprintf('~/pp/data_gla/fw4bt/osfstorage/data/gla01/pupil/sub0%d_gla_pp.mat',isubj));
      else
        load(sprintf('~/pp/data_gla/fw4bt/osfstorage/data/gla01/pupil/sub%d_gla_pp.mat',isubj));
      end
      
      pup = data.trial{1}';
      pup = pup(:,4);
      
      % load meg data
      if isubj < 10
        load(sprintf('~/pp/data_gla/fw4bt/osfstorage/data/gla01/meg/sub0%d_gla_meg.mat',isubj));
      else
        load(sprintf('~/pp/data_gla/fw4bt/osfstorage/data/gla01/meg/sub%d_gla_meg.mat',isubj));
      end
      
      
    catch me
      src_r = nan(8799,25);
      save([outdir fn '.mat'],'src_r')
      continue
    end
    
   
    dat = data.trial{1}; clear data
    pup_shift = round(400*0.93); % 930s from hoeks and levelt (1992?)    %     

    for ipupfreq = 1 : length(pupil_freqoi)
      pupil(:,ipupfreq) = ft_preproc_bandpassfilter(pup', 400, [pupil_freqoi(ipupfreq,1) pupil_freqoi(ipupfreq,2)], 1, 'but', [], 'no');
    end
    % ----------------
    pupil(end+1:end+pup_shift-1,:)=nan;
    dat(:,end-pup_shift:end)=nan;
    
    clear pup
    
    if isubj < 10
      load(sprintf('~/pp/data_gla/fw4bt/osfstorage/data/gla01/leadfields/sub0%d_gla_lf_BNA5mm.mat',isubj))
    else
      load(sprintf('~/pp/data_gla/fw4bt/osfstorage/data/gla01/leadfields/sub%d_gla_lf_BNA5mm.mat',isubj))
    end
    
    for ifreq=1:numel(freqoi)
      ifreq

      para          = [];
      para.freq     = freqoi(ifreq);
      para.fsample  = 400;  
      [csd, dataf,opt]=tp_compute_csd_wavelets(dat,para);
      nseg=floor((size(dat,2)-opt.n_win)/opt.n_shift+1);

      for ipupfreq = 1 : size(pupil_freqoi,1)
        
        for j=1:nseg
          dloc2=dat(:,(j-1)*opt.n_shift+1:(j-1)*opt.n_shift+opt.n_win)';
          if any(isnan(dloc2(:,1)))
%             warning('NaN detected')
            pup{ipupfreq}{ifreq}(:,j) = nan;
            continue
          else
            tmp = pupil((j-1)*opt.n_shift+1:(j-1)*opt.n_shift+opt.n_win,ipupfreq);
            pup{ipupfreq}{ifreq}(:,j) = mean(tmp.*gausswin(opt.n_win,3));
          end
        end
      end
      
      % Compute cross spectrum & find NaN segments
      idx_valid     = find(~isnan(dataf(1,:)));
      nanidx{ifreq} = idx_valid;
      % -----------
      % beamforming (ignore name of leadfield)
      % --------------
      para     	= [];
      para.reg  = 0.05;
      filt      = tp_beamformer(real(csd),lf,para);
      % --------------
      src_pow{ifreq} = abs((filt'*dataf).^2); % compute power
      % -----------

      clear src
    end
    
    % compute pupil diameter for segments
    % ---------------
    for ipup = 1 : size(pupil,2)
      ipup
      for ifreq = 1 : size(src_pow,2)
        src_r(:,ifreq,ipup) = corr(src_pow{ifreq}(:,nanidx{ifreq})',pup{ipup}{ifreq}(nanidx{ifreq})');
      end
    end
    
    save([outdir fn '.mat'],'src_r')
    tp_parallel(fn,outdir,0)
    
    clear data pupil pup dat tf src_r nanidx outp
  end
  
end


error('!')
%%