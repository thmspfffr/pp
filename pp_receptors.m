clear rec_map

load ~/receptormaps/gene_values.mat
load(sprintf('~/pp/proc/pp_atlas_BNA.mat'))
addpath /home/gnolte/meg_toolbox/toolbox_nightly/


orig = [46 64 37];
locs=2*(locs-repmat(orig,[size(locs,1) 1]));
idx=BNA.centroids(:,1)>0;

load /home/gnolte/meth/templates/sa_template.mat;

% g2    = sa_template.grid_cortex3000;
% idx = g2(:,1)>0;
% g2    = BNA.centroids(idx,:);
g2    = BNA.centroids;
g1    = locs;
dd    = 0.75;

for i = 1 : 54
rec_map(:,i) = spatfiltergauss(zscore(gdat(:,i)),g1,dd,g2);
end

rec_map = rec_map(idx,:);
%%
SUBJLIST  = [4 5 6 7 8 9 10 11 12 13 15 16 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34];
v = 2;
clear all_corr
for isubj = SUBJLIST
  
  for iblock = 1 : 2
    clear src_r
    try
      
      load(sprintf([outdir 'pp_cnt_src_pupil_power_correlations_crossfreq_s%d_b%d_v%d.mat'],isubj,iblock,v));
      all_corr(:,:,:,isubj,2,iblock) = src_r(:,:,1:19);
      
      load(sprintf([outdir 'pp_src_pupil_power_correlations_crossfreq_s%d_b%d_v%d.mat'],isubj,iblock,v));
      all_corr(:,:,:,isubj,1,iblock) = src_r(:,:,1:19);

      if ~exist('src_r','var')
        save(sprintf([outdir 'pp_cnt_src_pupil_power_correlations_crossfreq_s%d_b%d_v%d.mat'],isubj,iblock,v),'src_r');
      end
    catch me
      all_corr(:,:,:,isubj,1:2,iblock) = nan(246,25,19,2);
      continue
    end
    
  end

end

all_corr= nanmean(all_corr(:,:,:,SUBJLIST,:,:),6);

% first meg freq, then pupil freq
% 122 x 25 (meg) x 25 (pup
%%


for ifoi = 1 : 25
  ifoi
  for ipupfreq = 1 : 19
    for isubj = 1:28
      tmp = abs(all_corr(:,ifoi,ipupfreq,isubj,1));
%       tmp = spatfiltergauss(tmp,BNA.centroids,dd,g2);
      tmp = corr(repmat(tmp(idx,:),[1 30]),rec_map);
      
      r(:,ifoi,ipupfreq,isubj)  = tmp(1,:);
%       for igdat = 1  : 30
% 
%         r(igdat,ifoi,ipupfreq,isubj) = corr(tmp,rec_map(:,igdat));
%       end
    end
  end
end


%%
clim = [-0.1 0.1]

NA1a = squeeze(mean(r(1:3,:,:,:),1));
NA2a = squeeze(mean(r(4:6,:,:,:),1));
NA1b = squeeze(mean(r(7:9,:,:,:),1));
AChN = squeeze(mean(r([20 21 22 23 24 25 26 27],:,:,:),1));
AChM = squeeze(mean(r([15:19],:,:,:),1));
D1 = squeeze(mean(r([10 14],:,:,:),1));
D2 = squeeze(mean(r([11 12 13],:,:,:),1));

figure; set(gca,'color','w'); hold on
subplot(2,5,1);
[h,p]=ttest(NA1a,zeros(size(NA1a)),'dim',3); h = p<fdr1(p(:),0.05);
imagesc(nanmean(NA1a,3)',clim)
set(gca,'xtick',[1 5 9 13 17 21 25],'xticklabel',num2cell(freqoi([1 5 9 13 17 21 25])),'ydir','normal')
set(gca,'ytick',1:2:19,'yticklabel',num2cell(round(mean(pupil_freqoi(1:2:19,:),2).*1000)/1000),'ydir','normal')
tp_editplots; xlabel('Carrier frequency (MEG)'); ylabel('Pupil frequency')
title('NE - a1')

subplot(2,5,6);
imagesc(nanmean(NA1a,3)'.*h',clim)
set(gca,'xtick',[1 5 9 13 17 21 25],'xticklabel',num2cell(freqoi([1 5 9 13 17 21 25])),'ydir','normal')
set(gca,'ytick',1:2:19,'yticklabel',num2cell(round(mean(pupil_freqoi(1:2:19,:),2).*1000)/1000),'ydir','normal')
tp_editplots; xlabel('Carrier frequency (MEG)'); ylabel('Pupil frequency')

subplot(2,5,2);
[h,p]=ttest(NA2a,zeros(size(NA2a)),'dim',3); h = p<fdr1(p(:),0.05);
imagesc(nanmean(NA2a,3)',clim)
set(gca,'xtick',[1 5 9 13 17 21 25],'xticklabel',num2cell(freqoi([1 5 9 13 17 21 25])),'ydir','normal')
set(gca,'ytick',1:2:19,'yticklabel',num2cell(round(mean(pupil_freqoi(1:2:19,:),2).*1000)/1000),'ydir','normal')
tp_editplots; xlabel('Carrier frequency (MEG)'); ylabel('Pupil frequency')
title('NE - a2')

subplot(2,5,7);
imagesc(nanmean(NA2a,3)'.*h',clim)
set(gca,'xtick',[1 5 9 13 17 21 25],'xticklabel',num2cell(freqoi([1 5 9 13 17 21 25])),'ydir','normal')
set(gca,'ytick',1:2:19,'yticklabel',num2cell(round(mean(pupil_freqoi(1:2:19,:),2).*1000)/1000),'ydir','normal')
tp_editplots; xlabel('Carrier frequency (MEG)'); ylabel('Pupil frequency')

subplot(2,5,3);
[h,p]=ttest(NA1b,zeros(size(NA1b)),'dim',3); h = p<fdr1(p(:),0.05);
imagesc(nanmean(NA1b,3)',clim)
set(gca,'xtick',[1 5 9 13 17 21 25],'xticklabel',num2cell(freqoi([1 5 9 13 17 21 25])),'ydir','normal')
set(gca,'ytick',1:2:19,'yticklabel',num2cell(round(mean(pupil_freqoi(1:2:19,:),2).*1000)/1000),'ydir','normal')
tp_editplots; xlabel('Carrier frequency (MEG)'); ylabel('Pupil frequency')
title('NE - b1')

subplot(2,5,8);
imagesc(nanmean(NA1b,3)'.*h',clim)
set(gca,'xtick',[1 5 9 13 17 21 25],'xticklabel',num2cell(freqoi([1 5 9 13 17 21 25])),'ydir','normal')
set(gca,'ytick',1:2:19,'yticklabel',num2cell(round(mean(pupil_freqoi(1:2:19,:),2).*1000)/1000),'ydir','normal')
tp_editplots; xlabel('Carrier frequency (MEG)'); ylabel('Pupil frequency')

subplot(2,5,4);
[h,p]=ttest(AChN,zeros(size(AChN)),'dim',3); h = p<fdr1(p(:),0.05);
imagesc(nanmean(AChN,3)',clim)
set(gca,'xtick',[1 5 9 13 17 21 25],'xticklabel',num2cell(freqoi([1 5 9 13 17 21 25])))
title('ACh - N')

subplot(2,5,9);
imagesc(nanmean(AChN,3)'.*h',clim)
set(gca,'xtick',[1 5 9 13 17 21 25],'xticklabel',num2cell(freqoi([1 5 9 13 17 21 25])),'ydir','normal')
set(gca,'ytick',1:2:19,'yticklabel',num2cell(round(mean(pupil_freqoi(1:2:19,:),2).*1000)/1000),'ydir','normal')
tp_editplots; xlabel('Carrier frequency (MEG)'); ylabel('Pupil frequency')

subplot(2,5,5);
[h,p]=ttest(AChM,zeros(size(AChM)),'dim',3); h = p<fdr1(p(:),0.05);
imagesc(nanmean(AChM,3)',clim)
set(gca,'xtick',[1 5 9 13 17 21 25],'xticklabel',num2cell(freqoi([1 5 9 13 17 21 25])),'ydir','normal')
set(gca,'ytick',1:2:19,'yticklabel',num2cell(round(mean(pupil_freqoi(1:2:19,:),2).*1000)/1000),'ydir','normal')
tp_editplots; xlabel('Carrier frequency (MEG)'); ylabel('Pupil frequency')
title('ACh - M')

subplot(2,5,10);
imagesc(nanmean(AChM,3)'.*h',clim)
set(gca,'xtick',[1 5 9 13 17 21 25],'xticklabel',num2cell(freqoi([1 5 9 13 17 21 25])),'ydir','normal')
set(gca,'ytick',1:2:19,'yticklabel',num2cell(round(mean(pupil_freqoi(1:2:19,:),2).*1000)/1000),'ydir','normal')
tp_editplots; xlabel('Carrier frequency (MEG)'); ylabel('Pupil frequency')
% 
% subplot(2,5,4);
% [h,p]=ttest(D1,zeros(size(D1)),'dim',3); h = p<fdr1(p(:),0.05);
% imagesc(nanmean(D1,3)',[-0.05 0.05])
% set(gca,'xtick',[1 5 9 13 17 21 25],'xticklabel',num2cell(freqoi([1 5 9 13 17 21 25])),'ydir','normal')
% set(gca,'ytick',1:2:19,'yticklabel',num2cell(round(mean(pupil_freqoi(1:2:19,:),2).*1000)/1000),'ydir','normal')
% tp_editplots; xlabel('Carrier frequency (MEG)'); ylabel('Pupil frequency')
% title('D1')
% 
% subplot(2,5,9);
% imagesc(nanmean(D1,3)'.*h',[-0.05 0.05])
% set(gca,'xtick',[1 5 9 13 17 21 25],'xticklabel',num2cell(freqoi([1 5 9 13 17 21 25])),'ydir','normal')
% set(gca,'ytick',1:2:19,'yticklabel',num2cell(round(mean(pupil_freqoi(1:2:19,:),2).*1000)/1000),'ydir','normal')
% tp_editplots; xlabel('Carrier frequency (MEG)'); ylabel('Pupil frequency')
% 
% subplot(2,5,5);
% [h,p]=ttest(D2,zeros(size(D2)),'dim',3); h = p<fdr1(p(:),0.05);
% imagesc(nanmean(D2,3)',[-0.05 0.05])
% set(gca,'xtick',[1 5 9 13 17 21 25],'xticklabel',num2cell(freqoi([1 5 9 13 17 21 25])),'ydir','normal')
% set(gca,'ytick',1:2:19,'yticklabel',num2cell(round(mean(pupil_freqoi(1:2:19,:),2).*1000)/1000),'ydir','normal')
% tp_editplots; xlabel('Carrier frequency (MEG)'); ylabel('Pupil frequency')
% title('D2')
% 
% subplot(2,5,10);
% imagesc(nanmean(D2,3)'.*h',[-0.05 0.05])
% set(gca,'xtick',[1 5 9 13 17 21 25],'xticklabel',num2cell(freqoi([1 5 9 13 17 21 25])),'ydir','normal')
% set(gca,'ytick',1:2:19,'yticklabel',num2cell(round(mean(pupil_freqoi(1:2:19,:),2).*1000)/1000),'ydir','normal')
% tp_editplots; xlabel('Carrier frequency (MEG)'); ylabel('Pupil frequency')
% 

%%
i = 15
figure; set(gca,'color','w'); hold on

null = zeros(size(r));
[~,~,~,s]=ttest(squeeze(r(i,:,:,:)),squeeze(null(i,:,:,:)),'dim',3);
subplot(1,2,1)
imagesc(s.tstat'.*h',[-2 2]); colormap(cmap); 
set(gca,'xtick',[1 5 9 13 17 21 25],'xticklabel',num2cell(freqoi([1 5 9 13 17 21 25])),'ydir','normal')
tp_editplots; xlabel('Carrier frequency (MEG)'); ylabel('Pupil frequency')
set(gca,'ytick',1:2:19,'yticklabel',num2cell(round(mean(pupil_freqoi(1:2:19,:),2).*1000)/1000),'ydir','normal')
