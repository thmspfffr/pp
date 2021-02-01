


%% REGIONS OF INTEREST 

M1 = [-42 -26 54; 38 -32 48];
V1 = [-20 -86 18; 16 -80 26];
A1 = [-54 -22 10; 52 -24 12];

for i= 1 : size(BNA.grid_5mm,1)
  
  d_A1_l(i) = sqrt((BNA.grid_5mm(i,1)-A1(1,1))^2 + (BNA.grid_5mm(i,2)-A1(1,2))^2 + (BNA.grid_5mm(i,3)-A1(1,3))^2);
  d_A1_r(i) = sqrt((BNA.grid_5mm(i,1)-A1(2,1))^2 + (BNA.grid_5mm(i,2)-A1(2,2))^2 + (BNA.grid_5mm(i,3)-A1(2,3))^2);
  d_V1_l(i) = sqrt((BNA.grid_5mm(i,1)-V1(1,1))^2 + (BNA.grid_5mm(i,2)-V1(1,2))^2 + (BNA.grid_5mm(i,3)-V1(1,3))^2);
  d_V1_r(i) = sqrt((BNA.grid_5mm(i,1)-V1(2,1))^2 + (BNA.grid_5mm(i,2)-V1(2,2))^2 + (BNA.grid_5mm(i,3)-V1(2,3))^2);
  d_M1_l(i) = sqrt((BNA.grid_5mm(i,1)-M1(1,1))^2 + (BNA.grid_5mm(i,2)-M1(1,2))^2 + (BNA.grid_5mm(i,3)-M1(1,3))^2);
  d_M1_r(i) = sqrt((BNA.grid_5mm(i,1)-M1(2,1))^2 + (BNA.grid_5mm(i,2)-M1(2,2))^2 + (BNA.grid_5mm(i,3)-M1(2,3))^2);
  
end

[~,A1_l]=min(d_A1_l);
[~,A1_r]=min(d_A1_r);
[~,V1_l]=min(d_V1_l);
[~,V1_r]=min(d_V1_r);
[~,M1_l]=min(d_M1_l);
[~,M1_r]=min(d_M1_r);
dlpfc_l=BNA.tissue_5mm==3;
dlpfc_r=BNA.tissue_5mm==4;


%% Figure X 
% -------------------------------------
% Plot Fig. X, showing pupil-MEG correlations across source space
% Panel A: Pooled across all data sets and locations
% Panel B: For all 246 regions of the BNA atlas
% Panel C: Separately for 4 regions of interest
% -------------------------------------

[plt_gla,plt_hh,plt_mue,plt_hh_cnt,plt_all]=pp_load_results(2);

[~,idx_sorted] = sort(BNA.centroids(:,2),'descend');

figure_w; set(gcf,'Position',[500 500 750 350])

subplot(8,4,[1 5]); hold on;

% pool across datasets and compute cluster-based permutation statistics
tmp = squeeze(mean(cat(3,plt_gla.corr_src,plt_hh.corr_src,plt_mue.corr_src),1));
[c,p]=permutest(tmp,zeros(size(squeeze(tmp))),1,0.01,1000,2); h=[c{p<0.05}];
% Mean and SEM (across subjectds)
par = mean(tmp,2); 
par_std = std(tmp,[],2)/sqrt(size(tmp,2));

% plot pupil-MEG correlation, averaged across source locations
shadedErrorBar(log10(freqoi),par,par_std,'k')
plot(log10(freqoi(h)),par(h),'k.','markersize',8)
axis([.3 2.11 -0.05 0.05])
line([.3 2.11], [0 0],'color',[.7 .7 .7],'linestyle',':')
tp_editplots; ylabel('Mean corr.','fontsize',6)
set(gca,'xtick',log10(freqoi(1:4:25)),'xticklabel',''); 
text(log10(3),0.04,'o P<0.05 (corrected)','fontsize',6); title('Pooled')
set(gca,'fontsize',6)

% visual cortex
subplot(8,4,[2 6]); hold on;

par = mean(plt_all.corr_src([V1_l V1_r],:,:),1);
[c,p]=permutest(squeeze(par),zeros(size(squeeze(par))),1,0.01,1000,2); h=[c{p<0.05}];
shadedErrorBar(log10(freqoi),nanmean(par,3),std(par,[],3)/sqrt(size(par,3)),'k')
if sum(h)>0
  plot(log10(freqoi(h)),nanmean(par(:,h,:),3),'k.','markersize',8)
end
line([.3 2.11], [0 0],'color',[.7 .7 .7],'linestyle',':')
tp_editplots; ylabel('Mean corr.','fontsize',6)
set(gca,'xtick',log10(freqoi(1:4:25)),'xticklabel',{''})
axis([.3 2.11 -0.05 0.05]); 
text(log10(3),0.04,'Visual cortex','fontsize',6)
set(gca,'fontsize',6)

% auditory cortex
subplot(8,4,[10 14]); hold on;

par = mean(plt_all.corr_src([A1_l A1_r],:,:),1);
[c,p]=permutest(squeeze(par),zeros(size(squeeze(par))),1,0.01,1000,2); h=[c{p<0.05}];
shadedErrorBar(log10(freqoi),nanmean(par,3),std(par,[],3)/sqrt(size(par,3)),'k')
if sum(h)>0
  plot(log10(freqoi(h)),nanmean(par(:,h,:),3),'k.','markersize',8)
end
line([.3 2.11], [0 0],'color',[.7 .7 .7],'linestyle',':')
tp_editplots; ylabel('Mean corr.','fontsize',6)
set(gca,'xtick',log10(freqoi(1:4:25)),'xticklabel',{''})
axis([.3 2.11 -0.05 0.05]); %title(sprintf('A1: r = %.3f',r_A1_lr))
text(log10(3),0.04,'Auditory cortex','fontsize',6)
set(gca,'fontsize',6)

% somatosensory cortex
subplot(8,4,[18 22]); hold on;

par = mean(plt_all.corr_src([M1_l M1_r],:,:),1);
[c,p]=permutest(squeeze(par),zeros(size(squeeze(par))),1,0.01,1000,2); h=[c{p<0.05}];
shadedErrorBar(log10(freqoi),nanmean(par,3),std(par,[],3)/sqrt(size(par,3)),'k')
if sum(h)>0
  plot(log10(freqoi(h)),nanmean(par(:,h,:),3),'k.','markersize',8)
end
line([.3 2.11], [0 0],'color',[.7 .7 .7],'linestyle',':')
tp_editplots; ylabel('Mean corr.','fontsize',6)
set(gca,'xtick',log10(freqoi(1:4:25)),'xticklabel',{''})
axis([.3 2.11 -0.05 0.05]); %title(sprintf('A1: r = %.3f',r_A1_lr))
text(log10(3),0.04,'Somatosensory cortex','fontsize',6)
set(gca,'fontsize',6)

% dorsolateral PFC
subplot(8,4,[26 30]); hold on;

par = mean(plt_all.corr_src(find([dlpfc_l+dlpfc_r]),:,:),1);
[c,p]=permutest(squeeze(par),zeros(size(squeeze(par))),1,0.01,1000,2); h=[c{p<0.05}];
shadedErrorBar(log10(freqoi),nanmean(par,3),std(par,[],3)/sqrt(size(par,3)),'k')
if sum(h)>0
  plot(log10(freqoi(h)),nanmean(par(:,h,:),3),'k.','markersize',8)
end
line([.3 2.11], [0 0],'color',[.7 .7 .7],'linestyle',':')
tp_editplots; xlabel('Frequency [Hz]','fontsize',6);ylabel('Mean corr.','fontsize',6)
set(gca,'xtick',log10(freqoi(1:4:25)),'xticklabel',round(freqoi(1:4:25)))
axis([.3 2.11 -0.05 0.05]); %title(sprintf('A1: r = %.3f',r_A1_lr))
text(log10(3),0.04,'Dorsolateral PFC','fontsize',6)
set(gca,'fontsize',6)

% imagesc plot, correlation across all BNA regions
subplot(8,4,[9 13 17 21 25 29]); 

pooled= cat(3,plt_hh.corr_src_BNA(idx_sorted,:,:),plt_gla.corr_src_BNA(idx_sorted,:,:),plt_mue.corr_src_BNA(idx_sorted,:,:));

[~,p]=ttest(pooled,zeros(size(pooled)),'dim',3); h = p<fdr1(p(:),0.05,1);
imagesc(nanmean(pooled,3).*h,[-0.05 0.05])
tp_editplots; colormap(cmap)
set(gca,'xtick',1:4:25,'xticklabel',round(freqoi(1:4:25)),'fontsize',6)
set(gca,'ytick',[0],'yticklabel',{''})
xlabel('Frequency [Hz]'); ylabel(sprintf('Brain regions\nPosterior <---- Anterior'),'fontsize',6)
set(gca,'box','off')
Pos = get(gca,'Position');
cb = colorbar; set(gca,'Position',Pos)
cb.Position=[0.29 0.1114 0.005 0.6029];
cb.Ticks=[-0.05 0 0.05];
set(gca,'fontsize',6)

% ---------------------------------
% LOAD DATA V=1 (no lag) AND PLOT TEMPORAL DERIVATIVE
% ---------------------------------

[plt_gla,plt_hh,plt_mue,plt_hh_cnt,plt_all]=pp_load_results(1);

[~,idx_sorted] = sort(BNA.centroids(:,2),'descend');

subplot(8,4,[3 7]); hold on;

tmp = squeeze(mean(cat(3,plt_gla.corr_src_df,plt_hh.corr_src_df,plt_mue.corr_src_df),1));
[c,p]=permutest(tmp,zeros(size(squeeze(tmp))),1,0.01,1000,2); h=[c{p<0.05}];
par = mean(tmp,2); 
par_std = std(tmp,[],2)/sqrt(size(tmp,2));

shadedErrorBar(log10(freqoi),par,par_std,'k')
plot(log10(freqoi(h)),par(h),'k.','markersize',8)
axis([.3 2.11 -0.05 0.05])
line([.3 2.11], [0 0],'color',[.7 .7 .7],'linestyle',':')
tp_editplots; ylabel('Mean corr.','fontsize',6)
set(gca,'xtick',log10(freqoi(1:4:25)),'xticklabel','')

text(log10(3),0.04,'o P<0.05 (corrected)','fontsize',6); title('Pooled')

subplot(8,4,[4 8]); hold on;

par = mean(plt_all.corr_src_df([V1_l V1_r],:,:),1);
[c,p]=permutest(squeeze(par),zeros(size(squeeze(par))),1,0.01,1000,2); h=[c{p<0.05}];
shadedErrorBar(log10(freqoi),nanmean(par,3),std(par,[],3)/sqrt(size(par,3)),'k')
if sum(h)>0
  plot(log10(freqoi(h)),nanmean(par(:,h,:),3),'k.','markersize',8)
end
line([.3 2.11], [0 0],'color',[.7 .7 .7],'linestyle',':')
tp_editplots; ylabel('Mean corr.','fontsize',6)
set(gca,'xtick',log10(freqoi(1:4:25)),'xticklabel',{''})
axis([.3 2.11 -0.05 0.05]); %title(sprintf('A1: r = %.3f',r_A1_lr))

text(log10(3),0.04,'Visual cortex','fontsize',6)
subplot(8,4,[12 16]); hold on;

par = mean(plt_all.corr_src_df([A1_l A1_r],:,:),1);
[c,p]=permutest(squeeze(par),zeros(size(squeeze(par))),1,0.01,1000,2); h=[c{p<0.05}];
shadedErrorBar(log10(freqoi),nanmean(par,3),std(par,[],3)/sqrt(size(par,3)),'k')
if sum(h)>0
  plot(log10(freqoi(h)),nanmean(par(:,h,:),3),'k.','markersize',8)
end
line([.3 2.11], [0 0],'color',[.7 .7 .7],'linestyle',':')
tp_editplots; ylabel('Mean corr.','fontsize',6)
set(gca,'xtick',log10(freqoi(1:4:25)),'xticklabel',{''})
axis([.3 2.11 -0.05 0.05]); %title(sprintf('A1: r = %.3f',r_A1_lr))

text(log10(3),0.04,'Auditory cortex','fontsize',6)

subplot(8,4,[20 24]); hold on;

par = mean(plt_all.corr_src_df([M1_l M1_r],:,:),1);
[c,p]=permutest(squeeze(par),zeros(size(squeeze(par))),1,0.01,1000,2); h=[c{p<0.05}];
shadedErrorBar(log10(freqoi),nanmean(par,3),std(par,[],3)/sqrt(size(par,3)),'k')
if sum(h)>0
  plot(log10(freqoi(h)),nanmean(par(:,h,:),3),'k.','markersize',8)
end
line([.3 2.11], [0 0],'color',[.7 .7 .7],'linestyle',':')
tp_editplots; ylabel('Mean corr.','fontsize',6)
set(gca,'xtick',log10(freqoi(1:4:25)),'xticklabel',{''})
axis([.3 2.11 -0.05 0.05]); %title(sprintf('A1: r = %.3f',r_A1_lr))

text(log10(3),0.04,'Somatosensory cortex','fontsize',6)

subplot(8,4,[28 32]); hold on;

par = mean(plt_all.corr_src_df(find([dlpfc_l+dlpfc_r]),:,:),1);
[c,p]=permutest(squeeze(par),zeros(size(squeeze(par))),1,0.01,1000,2); h=[c{p<0.05}];
shadedErrorBar(log10(freqoi),nanmean(par,3),std(par,[],3)/sqrt(size(par,3)),'k')
if sum(h)>0
  plot(log10(freqoi(h)),nanmean(par(:,h,:),3),'k.','markersize',8)
end
line([.3 2.11], [0 0],'color',[.7 .7 .7],'linestyle',':')
tp_editplots; xlabel('Frequency [Hz]','fontsize',6);ylabel('Mean corr.','fontsize',6)
set(gca,'xtick',log10(freqoi(1:4:25)),'xticklabel',round(freqoi(1:4:25)),'fontsize',6)
axis([.3 2.11 -0.05 0.05]); %title(sprintf('A1: r = %.3f',r_A1_lr))

text(log10(3),0.04,'Dorsolateral PFC','fontsize',6)

subplot(8,4,[11 15 19 23 27 31]); 

pooled= cat(3,plt_hh.corr_src_df_BNA(idx_sorted,:,:),plt_gla.corr_src_df_BNA(idx_sorted,:,:),plt_mue.corr_src_df_BNA(idx_sorted,:,:));

[~,p]=ttest(pooled,zeros(size(pooled)),'dim',3); h = p<fdr1(p(:),0.05,1);
imagesc(nanmean(pooled,3).*h,[-0.05 0.05])
tp_editplots; colormap(cmap)
set(gca,'xtick',1:4:25,'xticklabel',round(freqoi(1:4:25)),'fontsize',6)
set(gca,'ytick',[0],'yticklabel',{''})
xlabel('Frequency [Hz]','fontsize',6); ylabel(sprintf('Brain regions\nPosterior <---- Anterior'),'fontsize',6)
set(gca,'box','off')
Pos = get(gca,'Position');
cb = colorbar; set(gca,'Position',Pos)
cb.Position=[0.7 0.1114 0.005 0.6029];
cb.Ticks=[-0.05 0 0.05];
cb.Box='off';

print(gcf,'-dpdf',sprintf('~/pp/plots/figure3.pdf'))
