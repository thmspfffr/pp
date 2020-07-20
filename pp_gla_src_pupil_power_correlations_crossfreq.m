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
      src_r = nan(246,25);
      save([outdir fn '.mat'],'src_r')
      continue
    end
    
    len = min([size(data.trial{1},2) size(pup,1)]);
    data.trial{1} = data.trial{1}(:,1:len);
    
    dat = data.trial{1}; clear data
    pup = pup(1:len,:);
    pup_shift = round(400*0.93); % 930s from hoeks and levelt (1992?)    
    pup = pup(pup_shift:end); 
%     

    for ipupfreq = 1 : length(pupil_freqoi)
      pupil(:,ipupfreq) = ft_preproc_bandpassfilter(pup', 400, [pupil_freqoi(ipupfreq,1) pupil_freqoi(ipupfreq,2)], 1, 'but', [], 'no');
    end
    % ----------------
    pupil(end+1:end+pup_shift-1,:)=nan;
    dat(:,end-pup_shift:end)=nan;
    
    clear pup
    
    if isubj < 10
      load(sprintf('~/pp/data_gla/fw4bt/osfstorage/data/gla01/spatfilt/sub0%d_gla_sf_BNA5mm.mat',isubj))
    else
      load(sprintf('~/pp/data_gla/fw4bt/osfstorage/data/gla01/spatfilt/sub%d_gla_sf_BNA5mm.mat',isubj))
    end
    
    spat_filt = data; clear data

    for ifreq=1:numel(freqoi)
      ifreq
      % DEFINE WAVELETS HERE
      % --------------------------------
      [KERNEL,~,opt]=tp_mkwavelet(freqoi(ifreq), .5, 400);
      % --------------------------------
      
      nseg=floor((size(dat,2)-opt.n_win)/opt.n_shift+1);
      clear dataf
      % this can probably be coded more elegantly, but it gets job done
      for j=1:nseg
        
        dloc2=dat(:,(j-1)*opt.n_shift+1:(j-1)*opt.n_shift+opt.n_win)';
         
        if any(isnan(dloc2(:,1)))
%           warning('NaN detected')
          dataf(:,j)=nan(1,size(dloc2,2));
          continue
        else
          dataf(:,j)=dloc2'*KERNEL;
        end
      end
      
      for ipupfreq = 1 : size(pupil_freqoi,1)
        ipupfreq
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
      csd           = (dataf(:,idx_valid)*dataf(:,idx_valid)')/size(dataf(:,idx_valid ),2);
      nanidx{ifreq} = idx_valid;
      % -----------
      filt = spat_filt{ifreq}.spatFilt;
      src       = abs((filt*dataf).^2); % compute power
      % -----------
%       src = conv(src
      % compute average 
      for igrid = 1 : max(BNA.tissue_5mm(:))
        src_pow{ifreq}(:,igrid) = mean(src(BNA.tissue_5mm == igrid,nanidx{ifreq}));
      end
      clear src
    end
    
    % compute pupil diameter for segments
    % ---------------
    for ipup = 1 : size(pupil,2)
      ipup
      for ifreq = 1 : size(src_pow,2)
        src_r(:,ifreq,ipup) = corr(src_pow{ifreq}(nanidx{ifreq},:),pup{ipup}{ifreq}(nanidx{ifreq})');
      end
    end
    
    save([outdir fn '.mat'],'src_r')
    tp_parallel(fn,outdir,0)
    
    clear data pupil pup dat tf src_r nanidx
  end
  
end


error('!')
%%
addpath ~/Documents/MATLAB/cbrewer/cbrewer/
cmap = cbrewer('div', 'RdBu', 256,'pchip'); cmap = cmap(end:-1:1,:);
clear src_r_all

for isubj = SUBJLIST
  
  for iblock = 1 : 1
    
    try
      load(sprintf('/home/tpfeffer/pp/proc/src/pp_gla_src_pupil_power_correlations_crossfreq_s%d_b%d_v%d.mat',isubj,iblock,v))
      src_r_all(:,:,isubj,iblock) = squeeze(mean(src_r(:,:,2:end-3),1)); clear src_r
    catch me
      src_r_all(:,:,isubj,iblock) = nan(246,25);
    end
    
  end
end

src_r_all = nanmean(src_r_all(:,:,SUBJLIST,:),4);

figure_w;
subplot(1,2,1)
imagesc(nanmean(src_r_all,3)',[-0.02 0.02]);
set(gca,'ydir','normal');
pupil_freqoi(:,1) = 2.^(-9:(0.5):1);
pupil_freqoi(:,2) = 2.^(-8:(0.5):2);
set(gca,'ytick',[1:2:17],'yticklabel',round(mean(pupil_freqoi(2:2:end-3,:),2)*1000)/1000)
set(gca,'xtick',[1:4:25],'xticklabel',freqoi(1:4:end))
xlabel('MEG frequency [Hz]')
ylabel('Pupil frequency [Hz]')
tp_editplots

subplot(1,2,2)
h=ttest(src_r_all,zeros(size(src_r_all)),'dim',3)
imagesc(h'.*nanmean(src_r_all,3)',[-0.02 0.02]);
set(gca,'ydir','normal')
set(gca,'ytick',[1:2:17],'yticklabel',round(mean(pupil_freqoi(2:2:end-3,:),2)*1000)/1000)
set(gca,'xtick',[1:4:25],'xticklabel',freqoi(1:4:end))
xlabel('MEG frequency [Hz]')
tp_editplots
colormap(cmap)

print(gcf,'-dpdf',sprintf('~/pp/plots/pp_gla_src_crossfreq_v%d.pdf',v))

%%
clear src_r_all

for isubj = SUBJLIST
  
  for iblock = 1 : 2
    
    try
      load(sprintf('/home/tpfeffer/pp/proc/src/pp_src_pupil_power_correlations_crossfreq_s%d_b%d_v1.mat',isubj,iblock))
      src_r_all(:,:,:,isubj,iblock) = src_r(:,:,2:end-3); clear src_r
    catch me
      src_r_all(:,:,:,isubj,iblock) = nan(246,25,17);
    end
    
  end
end

src_r_all = nanmean(src_r_all(:,:,:,SUBJLIST,:),5);

pup_freq = 4;

[h,~,~,s] = ttest(squeeze(src_r_all(:,:,pup_freq,:)),zeros(size(squeeze(src_r_all(:,:,pup_freq,:)))),'dim',3,'alpha',0.01);
ltz = 100*sum(h>0&s.tstat>0)/246;
stz = 100*sum(h>0&s.tstat<0)/246;

figure_w;
subplot(3,3,[1 4]); 
imagesc(nanmean(src_r_all(:,:,pup_freq,:),4).*h,[-0.02 0.02]);
title(sprintf('Freq: %.3f Hz',pupil_freqoi(pup_freq)))
ylabel('BNA region'); xlabel('MEG freq. [Hz]')
set(gca,'xtick',[1:4:25],'xticklabel',freqoi(1:4:end))
tp_editplots

subplot(3,3,7); hold on
plot(ltz,'color','r')
plot(stz,'color','b')
axis([1 25 0 100])
set(gca,'xtick',[1:4:25],'xticklabel',freqoi(1:4:end))
set(gca,'ytick',[0 50 100],'yticklabel',[0 50 100])
tp_editplots

pup_freq = 8;
[h,~,~,s] = ttest(squeeze(src_r_all(:,:,pup_freq,:)),zeros(size(squeeze(src_r_all(:,:,pup_freq,:)))),'dim',3,'alpha',0.01);
ltz = 100*sum(h>0&s.tstat>0)/246;
stz = 100*sum(h>0&s.tstat<0)/246;

subplot(3,3,[2 5]); 
imagesc(nanmean(src_r_all(:,:,pup_freq,:),4).*h,[-0.02 0.02]);
title(sprintf('Freq: %.3f Hz',pupil_freqoi(pup_freq)))
xlabel('MEG freq. [Hz]')
set(gca,'xtick',[1:4:25],'xticklabel',freqoi(1:4:end))
tp_editplots

subplot(3,3,8); hold on
plot(ltz,'color','r')
plot(stz,'color','b')
axis([1 25 0 100])
set(gca,'xtick',[1:4:25],'xticklabel',freqoi(1:4:end))
set(gca,'ytick',[0 50 100],'yticklabel',[0 50 100])
tp_editplots

pup_freq = 16;
[h,~,~,s] = ttest(squeeze(src_r_all(:,:,pup_freq,:)),zeros(size(squeeze(src_r_all(:,:,pup_freq,:)))),'dim',3,'alpha',0.01);
ltz = 100*sum(h>0&s.tstat>0)/246;
stz = 100*sum(h>0&s.tstat<0)/246;

subplot(3,3,[3 6]); title(sprintf('Freq: %d Hz',pup_freq))
imagesc(nanmean(src_r_all(:,:,pup_freq,:),4).*h,[-0.02 0.02]);
title(sprintf('Freq: %.3f Hz',pupil_freqoi(pup_freq)))
xlabel('MEG freq. [Hz]')
set(gca,'xtick',[1:4:25],'xticklabel',freqoi(1:4:end))
tp_editplots

subplot(3,3,9); hold on
plot(ltz,'color','r')
plot(stz,'color','b')
axis([1 25 0 100])
set(gca,'xtick',[1:4:25],'xticklabel',freqoi(1:4:end))
set(gca,'ytick',[0 50 100],'yticklabel',[0 50 100])
tp_editplots

colormap(cmap)

print(gcf,'-dpdf',sprintf('~/pp/plots/pp_gla_src_crossfreq_examples_v%d.pdf',v))
