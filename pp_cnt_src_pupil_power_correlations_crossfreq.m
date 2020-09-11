%% pp_cnt_src_pupil_power_correlations_crossfreq
% correlate bandpass filtered (or via wavelets) pupil and MEG signals

clear
restoredefaultpath

SUBJLIST  = [4 5 6 7 8 9 10 11 12 13 15 16 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34];
% pupil lag not accounted for
v = 1;
% pupil lag accounted for
% v = 3;
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
  
  for iblock = 1:2
    % %
    fn = sprintf('pp_cnt_src_pupil_power_correlations_crossfreq_s%d_b%d_v%d',isubj,iblock,v);
    if tp_parallel(fn,outdir,1,0)
      continue
    end
    %
    fprintf('Processing subj%d block%d ...\n',isubj,iblock);
    
    try
      % load meg data
      load(sprintf('~/pp/proc/pp_cnt_sens_cleandat_s%d_m%d_b%d_v%d.mat',isubj,im,iblock,1))
      % load pupil data
      load(sprintf('~/pp/proc/pup/pp_pupil_diameter_cleaned_counting_s%d_m%d_b%d.mat',isubj,im,iblock))
      
    catch me
      sprintf('%s',me.message)
      src_r = nan(246,length(freqoi),length(pupil_freqoi));
      save([outdir fn '.mat'],'src_r')
      continue
    end
    
    data.trial = data.trial(:,1:data.end_of_recording);
    data.trial = data.trial(:,end:-1:1);
    
    if size(pupil,2)>1
      raw_pupil=pupil(:,4);
    else
      raw_pupil=pupil;
    end
    
    
    clear pupil nanidx
    
    
    for ipupfreq = 1 : size(pupil_freqoi,1)
      
      f_sample = 1000;
      
      pup = resample(raw_pupil,400,1000);
      pup = pup(end:-1:1);
      
      if ipupfreq== 1
        len = min([size(pup,1) size(data.trial,2)]);
        if len/400 > 600
          len = 400*600;
        end
      end
      
      
      pup = pup(1:len);
      pup = pup(end:-1:1);
      
      pup_shift = round(400*0.93); % 930s from hoeks and levelt (1992?)
      pup = pup(pup_shift:end);
      
      pup = tp_pupil_bandpass(pup,[pupil_freqoi(ipupfreq,1) pupil_freqoi(ipupfreq,2)],f_sample);
      
      pup(end+1:end+pup_shift-1)=nan;
      
      if ipupfreq== 1
        data.trial = data.trial(:,1:len);
        data.trial = data.trial(:,end:-1:1);
        data.trial(:,isnan(pup))=nan;
      end
      % ----------------
      % THIS NEEDS TO BE CHANGED
      % ----------------
      pupil(:,ipupfreq) = pup; clear pup
    end
    % ----------------
    
    load(['/home/tpfeffer/pp/proc/src/' sprintf('pp_cnt_sa_s%d_m%d_b%d_v%d.mat',isubj,im,iblock,1)],'sa');
    clear nanidx
    
    for ifreq=1:numel(freqoi)
      
      ifreq
      
      para          = [];
      para.freq     = freqoi(ifreq);
      para.fsample  = 400;
      %       para.meth     = 'conv';
      [csd, dataf,opt]=tp_compute_csd_wavelets(dat,para);
      
      %       conv_dataf = conv_dataf(:,round(tfPnts));
      
      outp.tp_sens_pow(:,ifreq) = diag(abs(csd));
      
      % -------------------------------
      % prepare pupil signal
      % -------------------------------
      nseg=floor((size(dat,2)-opt.n_win)/opt.n_shift+1);
      pup = nan(nseg,1);
      for j=1:nseg
        tmp = pupil((j-1)*opt.n_shift+1:(j-1)*opt.n_shift+opt.n_win);
        pup(j) = mean(tmp.*gausswin(opt.n_win,3));
      end
      
      for ipupfreq = 1 : size(pupil_freqoi,1)
        ipupfreq
        for j=1:nseg
          dloc2=dat(:,(j-1)*n_shift+1:(j-1)*n_shift+n_win)';
          if any(isnan(dloc2(:,1)))
            %             warning('NaN detected')
            pup{ipupfreq}{ifreq}(:,j) = nan;
            continue
          else
            tmp = pupil((j-1)*n_shift+1:(j-1)*n_shift+n_win,ipupfreq);
            pup{ipupfreq}{ifreq}(:,j) = mean(tmp.*gausswin(n_win,3));
          end
        end
      end
      
      % Compute cross spectrum & find NaN segments
      idx_valid = find(~isnan(dataf(1,:)));
      nanidx{ifreq} = idx_valid;
      % -----------
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
      
      % compute average
      for igrid = 1 : max(BNA.tissue_5mm(:))
        src_pow{ifreq}(:,igrid) = mean(src(BNA.tissue_5mm == igrid,:));
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
    
    
    clear data pupil pup dat tf src_r nanidx src_pow
  end
  
end


error('!')
%%
SUBJLIST  = [4 5 6 7 8 9 10 11 12 13 15 16 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34];
v = 6;
clear all_corr
for isubj = SUBJLIST
  isubj
  for iblock = 1 : 2
    clear src_r
    %     try
    %
    load(sprintf([outdir 'pp_cnt_src_pupil_power_correlations_crossfreq_s%d_b%d_v%d.mat'],isubj,iblock,v));
    all_corr(:,:,:,isubj,2,iblock) = src_r;
    %
    load(sprintf([outdir 'pp_src_pupil_power_correlations_crossfreq_s%d_b%d_v%d.mat'],isubj,iblock,v));
    all_corr(:,:,:,isubj,1,iblock) = src_r;
    
    if ~exist('src_r','var')
      save(sprintf([outdir 'pp_cnt_src_pupil_power_correlations_crossfreq_s%d_b%d_v%d.mat'],isubj,iblock,v),'src_r');
    end
    %     catch me
    %       all_corr(:,:,:,isubj,1:2,iblock) = nan(246,25,21,2);
    %       continue
    %     end
    
  end
  
end

all_corr= nanmean(all_corr(:,:,:,SUBJLIST,:,:),6);

%%
figure; set(gcf,'color','w','Position',[100 100 1200 800]);


icond = 1;

addpath ~/Documents/MATLAB/cbrewer/cbrewer/
[h,p,~,s]=ttest(zeros(size(squeeze(mean(all_corr(:,:,:,:,icond),1)))),squeeze(mean(all_corr(:,:,:,:,icond),1)),'dim',3);
% h = p<fdr1(p(:),0.05);
h=p<.01;

cmap = cbrewer('div','RdBu',128); cmap = cmap(end:-1:1,:)
src_r = mean(squeeze(mean(all_corr(:,:,:,:,icond),1)),3);

subplot(2,3,1);
imagesc(src_r(:,1:19)',[-0.02 0.02])
set(gca,'xtick',[1 5 9 13 17 21 25],'xticklabel',num2cell(freqoi([1 5 9 13 17 21 25])),'ydir','normal')
set(gca,'ytick',1:2:19,'yticklabel',num2cell(round(mean(pupil_freqoi(1:2:19,:),2).*1000)/1000),'ydir','normal')
tp_editplots; xlabel('Carrier frequency (MEG)'); ylabel('Pupil frequency')
title('Unmasked')
colormap(cmap); axis square
colorbar
subplot(2,3,4);
imagesc((src_r(:,1:19).*h(:,1:19))',[-0.02 0.02])
set(gca,'xtick',[1 5 9 13 17 21 25],'xticklabel',num2cell(freqoi([1 5 9 13 17 21 25])),'ydir','normal')
set(gca,'ytick',1:2:19,'yticklabel',num2cell(round(mean(pupil_freqoi(1:2:19,:),2).*1000)/1000),'ydir','normal')
tp_editplots; xlabel('Carrier frequency (MEG)'); ylabel('Pupil frequency')
title(sprintf('Masked at P < 0.05\n(FDR-corrected)'))
colormap(cmap); axis square
colorbar

icond = 2;

addpath ~/Documents/MATLAB/cbrewer/cbrewer/
[h,p,~,s]=ttest(zeros(size(squeeze(mean(all_corr(:,:,:,:,icond),1)))),squeeze(mean(all_corr(:,:,:,:,icond),1)),'dim',3);
% h = p<fdr1(p(:),0.05);
h=p<.01;

cmap = cbrewer('div','RdBu',128); cmap = cmap(end:-1:1,:)
src_r = mean(squeeze(mean(all_corr(:,:,:,:,icond),1)),3);

subplot(2,3,2);
imagesc(src_r(:,1:19)',[-0.02 0.02])
set(gca,'xtick',[1 5 9 13 17 21 25],'xticklabel',num2cell(freqoi([1 5 9 13 17 21 25])),'ydir','normal')
set(gca,'ytick',1:2:19,'yticklabel',num2cell(round(mean(pupil_freqoi(1:2:19,:),2).*1000)/1000),'ydir','normal')
tp_editplots; xlabel('Carrier frequency (MEG)'); ylabel('Pupil frequency')
title('Unmasked')
colormap(cmap); axis square
colorbar
subplot(2,3,5);
imagesc((src_r(:,1:19).*h(:,1:19))',[-0.02 0.02])
set(gca,'xtick',[1 5 9 13 17 21 25],'xticklabel',num2cell(freqoi([1 5 9 13 17 21 25])),'ydir','normal')
set(gca,'ytick',1:2:19,'yticklabel',num2cell(round(mean(pupil_freqoi(1:2:19,:),2).*1000)/1000),'ydir','normal')
tp_editplots; xlabel('Carrier frequency (MEG)'); ylabel('Pupil frequency')
title(sprintf('Masked at P < 0.05\n(FDR-corrected)'))
colormap(cmap); axis square
colorbar

[h,p,~,s]=ttest(squeeze(mean(all_corr(:,:,:,:,2),1)),squeeze(mean(all_corr(:,:,:,:,1),1)),'dim',3);
% h = p<fdr1(p(:),0.05);
h=p<.01;
cmap = cbrewer('div','RdBu',128); cmap = cmap(end:-1:1,:)
src_r = mean(squeeze(mean(all_corr(:,:,:,:,2),1)),3)-mean(squeeze(mean(all_corr(:,:,:,:,1),1)),3);

subplot(2,3,3);
imagesc(src_r(:,1:19)',[-0.02 0.02])
set(gca,'xtick',[1 5 9 13 17 21 25],'xticklabel',num2cell(freqoi([1 5 9 13 17 21 25])),'ydir','normal')
set(gca,'ytick',1:2:19,'yticklabel',num2cell(round(mean(pupil_freqoi(1:2:19,:),2).*1000)/1000),'ydir','normal')
tp_editplots; xlabel('Carrier frequency (MEG)'); ylabel('Pupil frequency')
title('Unmasked')
colormap(cmap); axis square
colorbar
subplot(2,3,6);
imagesc((src_r(:,1:19).*h(:,1:19))',[-0.02 0.02])
set(gca,'xtick',[1 5 9 13 17 21 25],'xticklabel',num2cell(freqoi([1 5 9 13 17 21 25])),'ydir','normal')
set(gca,'ytick',1:2:19,'yticklabel',num2cell(round(mean(pupil_freqoi(1:2:19,:),2).*1000)/1000),'ydir','normal')
tp_editplots; xlabel('Carrier frequency (MEG)'); ylabel('Pupil frequency')
title(sprintf('Masked at P < 0.05\n(FDR-corrected)'))
colormap(cmap); axis square
colorbar

print(gcf,'-dpdf',sprintf('~/pp/plots/pp_src_pupil_power_correlations_crossfreq_v%d.pdf',v))



