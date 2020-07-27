%%

clear
restoredefaultpath

addpath ~/Documents/MATLAB/fieldtrip-20160919/
addpath ~/pconn/matlab/
load(sprintf('~/pp/proc/pp_atlas_BNA.mat'))
load /home/gnolte/meth/templates/mri.mat

ft_defaults

outdir = '~/pp/proc/src/';
ord    = pconn_randomization;

freqoi=2.^(1:(1/4):7); % 2-128 Hz as per Hipp et al. (2012) Nat Neurosci

SUBJLIST = [4 5 6 7 8 9 10 11 12 13 15 16 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34];

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


%% PLOT POWER (SENSOR LEVEL) - HAMBURG & GLASGOW
for isubj = SUBJLIST
  isubj
  for iblock = 1 : 2
    clear src_r
    try
      load(sprintf([outdir 'pp_src_pupil_power_correlations_s%d_b%d_v%d.mat'],isubj,iblock,v));
      
      for ifreq = 1 : 25
        
        tmp = all_csd(:,:,ifreq); 
        pow_ham(ifreq,isubj,iblock) = mean(diag(abs(tmp)));
      end
      
     catch me
      pow_ham(ifreq,isubj,iblock) = nan;
    end
  end
end

pow_ham = nanmean(pow_ham(:,SUBJLIST,:),3);

for isubj = 1:24
  isubj
  for iblock = 1 : 2
    clear src_r
    try
      load(sprintf([outdir 'pp_gla_src_pupil_power_correlations_s%d_b%d_v%d.mat'],isubj,iblock,v));
      
      for ifreq = 1 : 25
        
        tmp = all_csd(:,:,ifreq); 
        pow_gla(ifreq,isubj,iblock) = mean(diag(abs(tmp)));
      end
      
     catch me
    end
  end
end

pow_gla= nanmean(pow_gla(:,1:24,:),3);
%%

subj_gla = [1 2 3 4 6:24];

figure_w; 
subplot(4,4,[1 2 3 4 5 6 7 8]); hold on
plot(log10(freqoi), mean(zscore(pow_ham,0,1),2))
plot(log10(freqoi), mean(zscore(pow_gla(:,subj_gla),0,1),2))
axis([0 2.5 -1.5 2.5]); xlabel('Frequency [Hz]'); 
tp_editplots
set(gca,'xtick',log10(freqoi(1:4:end)),'xticklabel',num2cell(freqoi(1:4:end)))

subplot(4,4,[9 10 13 14]); hold on
plot(log10(freqoi), log10(pow_ham),'k')
xlabel('Frequency [Hz]'); 
tp_editplots
set(gca,'xtick',log10(freqoi(1:4:end)),'xticklabel',num2cell(freqoi(1:4:end)))

subplot(4,4,[11 12 15 16]); hold on
plot(log10(freqoi), log10(pow_gla(:,subj_gla)),'k')
xlabel('Frequency [Hz]'); 
tp_editplots
set(gca,'xtick',log10(freqoi(1:4:end)),'xticklabel',num2cell(freqoi(1:4:end)))



%% 
clear all_corr all_pow
addpath /home/gnolte/meth/highlevel/
addpath ~/Documents/MATLAB/cbrewer/cbrewer/
cmap = cbrewer('div', 'RdBu', 256,'pchip'); cmap = cmap(end:-1:1,:);
load /home/gnolte/meth/templates/mri.mat

SUBJLIST= 1 : 24; SUBJLIST([5 9]) = [];
v = 3;
clear all_corr
i = 0
for isubj = SUBJLIST
  isubj
  for iblock = 1 : 1
    clear src_r
    try
      load(sprintf([outdir 'pp_gla_src_pupil_power_correlations_s%d_b%d_v%d.mat'],isubj,iblock,v));
      isubj
      corr_gla(:,:,isubj,iblock) = outp.src_r;
      nai_gla(:,:,isubj,iblock) = outp.src_nai;
    catch me
      warning('!!!')
      corr_gla(:,:,isubj,iblock) = nan(8799,25);
      nai_gla(:,:,isubj,iblock) =  nan(8799,25);
      continue
    end
  end
end

corr_gla= nanmean(corr_gla(:,:,SUBJLIST,:),4);
nai_gla= nanmean(nai_gla(:,:,SUBJLIST,:),4);

for igrid = 1 : max(BNA.tissue_5mm(:))
  corr_gla_BNA(igrid,:,:) = tanh(mean(atanh(corr_gla(BNA.tissue_5mm == igrid,:,:))));
  nai_gla_BNA(igrid,:,:) = mean(nai_gla(BNA.tissue_5mm == igrid,:,:));
end

SUBJLIST = [4 5 6 7 8 9 10 11 12 13 15 16 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34];
v = 1
for isubj = SUBJLIST
  isubj
  for iblock = 1 : 2
    clear src_r
    try
      load(sprintf([outdir 'pp_src_pupil_power_correlations_s%d_b%d_v%d.mat'],isubj,iblock,v));
      isubj
      corr_ham(:,:,isubj,iblock) = outp.src_r;
      nai_ham(:,:,isubj,iblock) = outp.src_nai;
    catch me
      warning('!!!')
      corr_ham(:,:,isubj,iblock) = nan(8799,25);
      nai_ham(:,:,isubj,iblock) =  nan(8799,25);
      continue
    end
  end
end

corr_ham = nanmean(corr_ham(:,:,SUBJLIST,:),4);
nai_ham  = nanmean(nai_ham(:,:,SUBJLIST,:),4);

for igrid = 1 : max(BNA.tissue_5mm(:))
  corr_ham_BNA(igrid,:,:) = tanh(mean(atanh(corr_ham(BNA.tissue_5mm == igrid,:,:))));
  nai_ham_BNA(igrid,:,:) = mean(nai_ham(BNA.tissue_5mm == igrid,:,:));
end

corr_all_BNA = cat(3,corr_gla_BNA,corr_ham_BNA);
nai_all_BNA  = cat(3,nai_gla_BNA,nai_ham_BNA);

corr_all = cat(3,corr_gla,corr_ham);
nai_all  = cat(3,nai_gla,nai_ham);


%% PLOT NAI / POWER
figure_w;

subplot(1,2,1);

imagesc(nanmean(nai_ham_BNA,3)./max(max(nanmean(nai_ham_BNA,3))),[0 .3])
colormap(inferno)

subplot(1,2,2);

imagesc(nanmean(nai_gla_BNA,3)./max(max(nanmean(nai_gla_BNA,3))),[0 .3])
colormap(inferno)


% 
% para.ifoi = [10];
% para.orientation = 'axial';
% para.dslice_shown = 0.5;
% 
% showmri_transp(mri,para,[BNA.centroids./10 nanmean(nanmean(nai_all_BNA(:,para.ifoi,:),2),3)])


%% CORR METHODS

clear all_corr all_pow
addpath /home/gnolte/meth/highlevel/
addpath ~/Documents/MATLAB/cbrewer/cbrewer/
cmap = cbrewer('div', 'RdBu', 256,'pchip'); cmap = cmap(end:-1:1,:);
load /home/gnolte/meth/templates/mri.mat

SUBJLIST= 1 : 24; SUBJLIST([5 9]) = [];
v = 3;
clear all_corr
for isubj = SUBJLIST
  isubj
  for iblock = 1 : 1
    clear src_r
    try
    load(sprintf([outdir 'pp_gla_src_pupil_power_correlations_s%d_b%d_v%d.mat'],isubj,iblock,v));
    corr_meth(:,:,isubj) = outp.corr_meth_src;
    catch me
      corr_meth(:,:,isubj) = nan(8799,25);
    end
  end
end

corr_meth  = corr_meth(:,:,SUBJLIST);
%%
figure_w;
subplot(2,7,[1 2 3 8 9 10])
imagesc(nanmean(corr_meth,3),[0.8 1]); 
colormap(plasma); colorbar
set(gca,'xtick',[(1:4:25)],'xticklabel',num2cell([freqoi(1:4:end)]))
xlabel('Frequency [Hz]'); ylabel('Source Region')
tp_editplots

subplot(2,7,[5 6 7])
plot(1:25,nanmean(nanmean(corr_meth,3),1),'k')
set(gca,'xtick',[(1:4:25)],'xticklabel',num2cell([freqoi(1:4:end)]))
axis([0 26 0 1.1]); tp_editplots;
xlabel('Frequency [Hz]'); ylabel('Avg. correlation')

