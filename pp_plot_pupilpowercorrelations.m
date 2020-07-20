%%
clear all_corr
v = 6;
for isubj = SUBJLIST
  for iblock = 1 : 2
            
    load(sprintf([outdir 'pp_src_pupil_power_correlations_s%d_b%d_v%d.mat'],isubj,iblock,v));
    all_corr(:,:,isubj,1,iblock) = src_r;
    load(sprintf([outdir 'pp_cnt_src_pupil_power_correlations_s%d_b%d_v%d.mat'],isubj,iblock,v));
    all_corr(:,:,isubj,2,iblock) = src_r;

  end
end

all_corr= nanmean(all_corr(:,:,SUBJLIST,:,:),5);
%
%% STATISTICS (simple t-test)
freqoi=2.^(1:(1/4):7); % 1-128 Hz as per Hipp et al. (2012) Nat Neurosci
addpath ~/Documents/MATLAB/cbrewer/cbrewer/
cmap = cbrewer('div','RdBu',128); cmap = cmap(end:-1:1,:)

[h,p]=ttest(zeros(size(all_corr)),all_corr,'dim',3)
h=p(:,:,:,1)<fdr1(reshape(p(:,:,:,1),[246*25 1]),0.05,0);

figure; set(gcf,'color','w')
k=subplot(2,2,[1 3]); hold on
imagesc(nanmean(all_corr(:,:,:,1),3),[-0.05 0.05])
set(gca,'xtick',[1 5 9 13 17 21 25],'xticklabel',num2cell(freqoi([1 5 9 13 17 21 25])),'ydir','reverse')
xlabel('Frequency [Hz]'); ylabel('BNA Region')
colormap(cmap); colorbar; axis([1 25 1 246]); tp_editplots

k=subplot(2,2,[2 4]); hold on
imagesc(nanmean(all_corr(:,:,:,1),3).*h,[-0.05 0.05])
set(gca,'xtick',[1 5 9 13 17 21 25],'xticklabel',num2cell(freqoi([1 5 9 13 17 21 25])),'ydir','reverse')
xlabel('Frequency [Hz]'); ylabel('BNA Region')
colormap(cmap); colorbar; axis([1 25 1 246]); tp_editplots

print(gcf,'-dpdf',sprintf('~/pp/plots/pp_src_pupil_power_correlations_v%d.pdf',v))

figure; set(gcf,'color','w')
k=subplot(2,2,[1 3]); hold on
h=p(:,:,:,2)<fdr1(reshape(p(:,:,:,2),[246*25 1]),0.05,0);
imagesc(nanmean(all_corr(:,:,:,2),3),[-0.05 0.05])
set(gca,'xtick',[1 5 9 13 17 21 25],'xticklabel',num2cell(freqoi([1 5 9 13 17 21 25])),'ydir','reverse')
xlabel('Frequency [Hz]'); ylabel('BNA Region')
colormap(cmap); colorbar; axis([1 25 1 246]); tp_editplots

k=subplot(2,2,[2 4]); hold on
imagesc(nanmean(all_corr(:,:,:,2),3).*h,[-0.05 0.05])
set(gca,'xtick',[1 5 9 13 17 21 25],'xticklabel',num2cell(freqoi([1 5 9 13 17 21 25])),'ydir','reverse')
xlabel('Frequency [Hz]'); ylabel('BNA Region')
colormap(cmap); colorbar; axis([1 25 1 246]); tp_editplots

print(gcf,'-dpdf',sprintf('~/pp/plots/pp_cnt_src_pupil_power_correlations_v%d.pdf',v))

%%

figure; set(gcf,'color','w')
[r,p]=corr(nanmean(all_corr(:,:,:,1),3),nanmean(all_corr(:,:,:,2),3))
imagesc(r.*(p<fdr1(p(:),0.05,0)),[-0.8 0.8])
set(gca,'ytick',[1 5 9 13 17 21 25],'yticklabel',num2cell(freqoi([1 5 9 13 17 21 25])),'ydir','reverse')
set(gca,'xtick',[1 5 9 13 17 21 25],'xticklabel',num2cell(freqoi([1 5 9 13 17 21 25])),'ydir','reverse')
xlabel('Frequency [Hz]'); ylabel('Frequency [Hz]'); axis square
tp_editplots; colormap(cmap); colorbar; set(gca,'ydir','normal')

figure; set(gcf,'color','w')
[r,p]=corr(nanmean(all_corr(:,:,:,1),3));
h=p<fdr1(p(:),0.05,0);
imagesc(r.*h,[-0.75 0.75])
set(gca,'ytick',[1 5 9 13 17 21 25],'yticklabel',num2cell(freqoi([1 5 9 13 17 21 25])),'ydir','reverse')
set(gca,'xtick',[1 5 9 13 17 21 25],'xticklabel',num2cell(freqoi([1 5 9 13 17 21 25])),'ydir','reverse')
xlabel('Frequency [Hz]'); ylabel('Frequency [Hz]'); axis square
tp_editplots; colormap(cmap); colorbar; set(gca,'ydir','normal')

figure; set(gcf,'color','w')
[r,p]=corr(nanmean(all_corr(:,:,:,2),3));
h=p<fdr1(p(:),0.05,0);
imagesc(r.*h,[-0.75 0.75])
set(gca,'ytick',[1 5 9 13 17 21 25],'yticklabel',num2cell(freqoi([1 5 9 13 17 21 25])),'ydir','reverse')
set(gca,'xtick',[1 5 9 13 17 21 25],'xticklabel',num2cell(freqoi([1 5 9 13 17 21 25])),'ydir','reverse')
xlabel('Frequency [Hz]'); ylabel('Frequency [Hz]'); axis square
tp_editplots; colormap(cmap); colorbar; set(gca,'ydir','normal')


    

%% PROJECT TO SURFACE
% 
% par = squeeze(nanmean(all_corr(:,freqoi>8&freqoi<18,:),2));
% par = nanmean(par,2);
% close all
% 
% cmap = cbrewer('div','RdBu',128); cmap = cmap(end:-1:1,:)
% 
% 
% para = [];
% para.clim = [-0.05 0.05];
% para.cmap = cmap;
% para.grid = BNA.centroids;
% para.dd = 1.5;
% para.fn = sprintf('~/pp/plots/pp_surface_v%d.png',v);
% tp_plot_surface(par,para)


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
