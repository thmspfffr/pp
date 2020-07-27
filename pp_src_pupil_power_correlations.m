%% pp_src_pupil_power_correlations
% correlate bandpass filtered (or via wavelets) pupil and MEG signals

clear
restoredefaultpath

% -------------------------
% VERSION 4
% -------------------------
v = 1;
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
      % load meg data
      load(sprintf('~/pupmod/proc/sens/pupmod_rest_sens_cleandat_s%d_m%d_b%d_v%d.mat',isubj,im,iblock,1))
      % load pupil data
      load(sprintf('~/pp/proc/pup/pp_pupil_diameter_cleaned_s%d_m%d_b%d.mat',isubj,im,iblock))
      
    catch me
      src_r = nan(246,25);
      save([outdir fn '.mat'],'src_r')
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
      
    dat(:,isnan(pupil))=nan(size(dat,1),sum(isnan(pupil)));
    
    load(['/home/tpfeffer/pp/proc/src/' sprintf('pp_sa_s%d_m%d_b%d_v%d.mat',isubj,im,iblock,1)],'sa');
    lf = sa.L_genemaps_aal;
    
    save(sprintf('~/pp/proc/sens/pp_meg_pupil_lfs_s%d_b%d.mat',isubj,iblock),'dat','pupil','lf');

    clear all_csd
    
    for ifreq=1:numel(freqoi)
      ifreq

      [KERNEL,~,opt]=tp_mkwavelet(freqoi(ifreq), .5, 400);
      
      nseg=floor((size(dat,2)-opt.n_win)/opt.n_shift+1);
      
      clear dataf pup csd
      %
      kk = 0;
      for j=1:nseg
        
        dloc2=dat(:,(j-1)*opt.n_shift+1:(j-1)*opt.n_shift+opt.n_win)';
        
        if any(isnan(dloc2(:,1)))
          warning('NaN detected')
          dataf(:,j)=nan(1,size(dloc2,2));
          pup(:,j) = nan;
          continue
        else
          kk = kk + 1;
          dataf(:,j)=dloc2'*KERNEL;
          tmp = pupil((j-1)*opt.n_shift+1:(j-1)*opt.n_shift+opt.n_win);
          pup(:,j) = mean(tmp.*gausswin(opt.n_win,3));
          if kk == 1
            csd=dataf(:,j)*dataf(:,j)';
          else
            csd=csd+dataf(:,j)*dataf(:,j)';
          end
        end
      end
      
      csd = csd/kk;
      
      all_csd(:,:,ifreq) = csd;

      % find indices of non-NAN segments
      idx_valid = find(~isnan(dataf(1,:)));
 
      % correlate pupil with sensor level signal
%       sens_r = corr(pup(:,idx_valid)',abs(dataf(:,idx_valid).^2)');
      
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
      
      % correlate with pupil
      outp.src_r(:,ifreq) = corr(pup(idx_valid)',src(:,idx_valid)');
      
      clear src

    end

    save([outdir fn '.mat'],'outp','all_csd')
    tp_parallel(fn,outdir,0)
    
    clear src_r all_nai  
    
  end
end

error('!')
%%
v = 1;
clear all_corr
i = 0
for isubj = SUBJLIST
  for iblock = 1 : 2
    clear src_r
    try
      load(sprintf([outdir 'pp_src_pupil_power_correlations_s%d_b%d_v%d.mat'],isubj,iblock,v));
      isubj
      all_corr(:,:,isubj,iblock) = outp.src_r;
      all_nai(:,:,isubj,iblock) = outp.src_nai;
    catch me
      all_corr(:,:,isubj,iblock) = nan(8799,25);
      all_nai(:,:,isubj,iblock) = nan(8799,25);
      continue
    end    
  end
end

all_corr= nanmean(all_corr(:,:,SUBJLIST,:),4);
all_nai = nanmean(all_nai(:,:,SUBJLIST,:),4);

for igrid = 1 : max(BNA.tissue_5mm(:))
  all_nai_BNA(igrid,:,:) = mean(all_nai(BNA.tissue_5mm == igrid,:,:));
end

%% PLOT POWER / NAI

figure; set(gcf,'color','w')
h=subplot(2,2,[1 3]); hold on
imagesc(nanmean(all_nai_BNA,3)./max(max(nanmean(all_nai_BNA,3))),[0 1])
set(gca,'xtick',[1 5 9 13 17 21 25],'xticklabel',num2cell(freqoi([1 5 9 13 17 21 25])),'ydir','reverse')
xlabel('Frequency [Hz]'); ylabel('BNA Region')
colormap(plasma); colorbar; axis([1 25 1 246]); tp_editplots

% h=subplot(2,2,[2 4]); hold on
% imagesc(nanmean(all_pow_BNA,3)./max(max(nanmean(all_pow_BNA,3))),[0 1])
% set(gca,'xtick',[1 5 9 13 17 21 25],'xticklabel',num2cell(freqoi([1 5 9 13 17 21 25])),'ydir','reverse')
% xlabel('Frequency [Hz]'); ylabel('BNA Region')
% colormap(plasma); colorbar; axis([1 25 1 246]); tp_editplots

print

%% STATISTICS (simple t-test)
addpath ~/Documents/MATLAB/cbrewer/cbrewer/
cmap = cbrewer('div','RdBu',128); cmap = cmap(end:-1:1,:)

[h,p]=ttest(zeros(size(all_corr)),all_corr,'dim',3)

figure; set(gcf,'color','w')
h=subplot(2,2,[1 3]); hold on
imagesc(nanmean(all_corr,3),[-0.05 0.05])
set(gca,'xtick',[1 5 9 13 17 21 25],'xticklabel',num2cell(freqoi([1 5 9 13 17 21 25])),'ydir','reverse')
xlabel('Frequency [Hz]'); ylabel('BNA Region')
colormap(cmap); colorbar; axis([1 25 1 246]); tp_editplots

h=subplot(2,2,[2 4]); hold on
imagesc(nanmean(all_corr,3).*(p<0.0157),[-0.05 0.05])
set(gca,'xtick',[1 5 9 13 17 21 25],'xticklabel',num2cell(freqoi([1 5 9 13 17 21 25])),'ydir','reverse')
xlabel('Frequency [Hz]'); ylabel('BNA Region')
colormap(cmap); colorbar; axis([1 25 1 246]); tp_editplots

print(gcf,'-dpdf',sprintf('~/pp/plots/pp_src_pupil_power_correlations_v%d.pdf',v))
%% PROJECT TO SURFACE

par = squeeze(nanmean(all_corr(:,freqoi>8&freqoi<18,:),2));
par = nanmean(par,2);
close all

cmap = cbrewer('div','RdBu',256); cmap = cmap(end:-1:1,:)


para = [];
para.clim = [-0.02 0.02];
para.cmap = cmap;
para.grid = BNA.centroids;
para.dd = 0.5;
para.fn = sprintf('~/pp/plots/pp_surface_v%d.png',v);
tp_plot_surface(par,para)






%%

for ifoi = 1 : 25
  
m_f(ifoi) = mean(squeeze(nanmean(all_corr(:,ifoi,:),1)));
s_f(ifoi) = nanstd(squeeze(nanmean(all_corr(:,ifoi,:),1)))/sqrt(size(all_corr,3));
[h(ifoi),p(ifoi)] = ttest(squeeze(nanmean(all_corr(:,ifoi,:),1)));

end



figure; set (gcf,'color','w')
subplot(2,2,1); hold on
shadedErrorBar([],m_f,[s_f],[]);
% plot(f,'linewidth',3);
line([1 25],[0 0],'linestyle',':','color','k')
set(gca,'xtick',[1 5 9 13 17 21 25],'xticklabel',num2cell(freqoi([1 5 9 13 17 21 25])))
xlabel('Frequency [Hz]'); ylabel('Correlation coeff.')
tp_editplots; title('Average over sensors')
axis([0 25 -0.04 0.04])
% subplot(1,2,2); hold on
% plot(t,'linewidth',3);
% line([1 25],[0 0],'linestyle',':','color','k')
% set(gca,'xtick',[1 5 9 13 17 21 25],'xticklabel',num2cell(freqoi([1 5 9 13 17 21 25])))
% xlabel('Frequency [Hz]'); ylabel('T-Value')
% tp_editplots

print(gcf,'-dpdf',sprintf('~/pp/plots/pp_sens_pupil_power_correlations_lineplot_v%d.pdf',v))
