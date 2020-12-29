


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

[~,idx_sorted] = sort(BNA.centroids(:,2),'descend');

figure_w; set(gcf,'Position',[500 500 400 400])
subplot(8,2,[1 3]); hold on;

tmp = squeeze(mean(cat(3,plt_gla.corr_src,plt_hh.corr_src,plt_mue.corr_src),1));
[c,p]=permutest(tmp,zeros(size(squeeze(tmp))),1,0.01,1000,2); h=[c{p<0.05}];

shadedErrorBar(log10(freqoi),par_all,std_all,'k')
plot(log10(freqoi(h)),par_all(h),'k.','markersize',8)
axis([.3 2.11 -0.05 0.05])
line([.3 2.11], [0 0],'color',[.7 .7 .7],'linestyle',':')
tp_editplots; ylabel('Mean corr.')
set(gca,'xtick',log10(freqoi(1:4:25)),'xticklabel','')

text(log10(3),0.04,'o P<0.05 (corrected)','fontsize',7); title('Pooled')

subplot(8,2,[2 4]); hold on;

par = mean(plt_all.corr_src([V1_l V1_r],:,:),1);
[c,p]=permutest(squeeze(par),zeros(size(squeeze(par))),1,0.01,1000,2); h=[c{p<0.05}];
shadedErrorBar(log10(freqoi),nanmean(par,3),std(par,[],3)/sqrt(size(par,3)),'k')
if sum(h)>0
  plot(log10(freqoi(h)),nanmean(par(:,h,:),3),'k.','markersize',8)
end
line([.3 2.11], [0 0],'color',[.7 .7 .7],'linestyle',':')
tp_editplots; ylabel('Mean corr.')
set(gca,'xtick',log10(freqoi(1:4:25)),'xticklabel',{''})
axis([.3 2.11 -0.05 0.05]); %title(sprintf('A1: r = %.3f',r_A1_lr))

text(log10(3),0.04,'Visual cortex','fontsize',7)
subplot(8,2,[6 8]); hold on;

par = mean(plt_all.corr_src([A1_l A1_r],:,:),1);
[c,p]=permutest(squeeze(par),zeros(size(squeeze(par))),1,0.01,1000,2); h=[c{p<0.05}];
shadedErrorBar(log10(freqoi),nanmean(par,3),std(par,[],3)/sqrt(size(par,3)),'k')
if sum(h)>0
  plot(log10(freqoi(h)),nanmean(par(:,h,:),3),'k.','markersize',8)
end
line([.3 2.11], [0 0],'color',[.7 .7 .7],'linestyle',':')
tp_editplots; ylabel('Mean corr.')
set(gca,'xtick',log10(freqoi(1:4:25)),'xticklabel',{''})
axis([.3 2.11 -0.05 0.05]); %title(sprintf('A1: r = %.3f',r_A1_lr))

text(log10(3),0.04,'Auditory cortex','fontsize',7)

subplot(8,2,[10 12]); hold on;

par = mean(plt_all.corr_src([M1_l M1_r],:,:),1);
[c,p]=permutest(squeeze(par),zeros(size(squeeze(par))),1,0.01,1000,2); h=[c{p<0.05}];
shadedErrorBar(log10(freqoi),nanmean(par,3),std(par,[],3)/sqrt(size(par,3)),'k')
if sum(h)>0
  plot(log10(freqoi(h)),nanmean(par(:,h,:),3),'k.','markersize',8)
end
line([.3 2.11], [0 0],'color',[.7 .7 .7],'linestyle',':')
tp_editplots; ylabel('Mean corr.')
set(gca,'xtick',log10(freqoi(1:4:25)),'xticklabel',{''})
axis([.3 2.11 -0.05 0.05]); %title(sprintf('A1: r = %.3f',r_A1_lr))

text(log10(3),0.04,'Somatosensory cortex','fontsize',7)

subplot(8,2,[14 16]); hold on;

par = mean(plt_all.corr_src(find([dlpfc_l+dlpfc_r]),:,:),1);
[c,p]=permutest(squeeze(par),zeros(size(squeeze(par))),1,0.01,1000,2); h=[c{p<0.05}];
shadedErrorBar(log10(freqoi),nanmean(par,3),std(par,[],3)/sqrt(size(par,3)),'k')
if sum(h)>0
  plot(log10(freqoi(h)),nanmean(par(:,h,:),3),'k.','markersize',8)
end
line([.3 2.11], [0 0],'color',[.7 .7 .7],'linestyle',':')
tp_editplots; xlabel('Frequency [Hz]');ylabel('Mean corr.')
set(gca,'xtick',log10(freqoi(1:4:25)),'xticklabel',round(freqoi(1:4:25)))
axis([.3 2.11 -0.05 0.05]); %title(sprintf('A1: r = %.3f',r_A1_lr))

text(log10(3),0.04,'Dorsolateral PFC','fontsize',7)

subplot(8,2,[5 7 9 11 13 15]); 

pooled= cat(3,plt_hh.corr_src_BNA(idx_sorted,:,:),plt_gla.corr_src_BNA(idx_sorted,:,:),plt_mue.corr_src_BNA(idx_sorted,:,:));

[~,p]=ttest(pooled,zeros(size(pooled)),'dim',3); h = p<fdr1(p(:),0.05,1);
imagesc(nanmean(pooled,3).*h,[-0.05 0.05])
tp_editplots; colormap(cmap)
set(gca,'xtick',1:4:25,'xticklabel',round(freqoi(1:4:25)))
set(gca,'ytick',[0],'yticklabel',{''})
xlabel('Frequency [Hz]'); ylabel(sprintf('Brain regions\nPosterior <---- Anterior'))
set(gca,'box','off')
