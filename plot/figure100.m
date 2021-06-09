
fooof22 = pp_load_fooof_results(2);
addpath ~/Documents/MATLAB/cbrewer/cbrewer/
%%
load ~/pp/proc/pp_atlas_BNA.mat
[~,idx_sorted] = sort(BNA.centroids(:,2),'descend');
colors = cbrewer('qual', 'Set3', 10,'pchip');
colors = colors(4:6,:);



%% SLOPES PER SEGMENT

pooled = squeeze(cat(5,nanmean(fooof22.slopes_seg_hh,6),fooof22.slopes_seg_gla,fooof22.slopes_seg_mue));

figure_w;

subplot(3,4,1); hold on

plot(1:20,squeeze(mean(mean(pooled,3),1)),'.','markersize',20)
lsline

tp_editplots; xlabel('Pupil bin'); ylabel('Spectral exponent')
axis([0 21 0.91 0.96])

print(gcf,'-dpdf',sprintf('~/pp/plots/pp_slope_vs_pupil.pdf',v))
% 
% figure_w
% % subplot(3,4,2); 
% hold on
% 
% isubj = 10;
% 
% tmp  = squeeze(nanmean(fooof22.pxx_seg_gla(5:100,160,:,:,:),5));
% 
% plot(log10(fooof22.fxx(5:100)),log10(squeeze(mean(tmp(1:end,1,isubj),3))))
% plot(log10(fooof22.fxx(5:100)),log10(squeeze(mean(tmp(1:end,20,isubj),3))))
%%

close

V1_r = find(BNA.tissue_5mm==205);
V1_l = find(BNA.tissue_5mm==206);
A1_r = find(BNA.tissue_5mm==71);
A1_l = find(BNA.tissue_5mm==72);
M1_r = find(BNA.tissue_5mm==159);
M1_l = find(BNA.tissue_5mm==160);
ACC_l= find(BNA.tissue_5mm==179);
ACC_r= find(BNA.tissue_5mm==180);


pooled = squeeze(cat(2,fooof22.slope_hh,fooof22.slope_gla,fooof22.slope_mue));


figure_w
subplot(5,3,[1 2 4 5]); hold on

r = (rand(81,1)-0.5)/4;

plot(ones(81,1)+r,squeeze(mean(pooled(:,:),1)),'k.','markersize',12)
% plot(1.4,squeeze(mean(mean(pooled(:,:),1),2)),'k.','markersize',20)
[~,p]=ttest(squeeze(mean(pooled(:,:),1)))

plot(2*ones(81,1)+r,squeeze(mean(pooled([205, 206],:),1)),'.','markersize',12)
% plot(2.4,squeeze(mean(mean(pooled([205, 206],:),1),2)),'.','markersize',20)

[~,p1]=ttest(squeeze(mean(pooled([205, 206],:),1)))

plot(3*ones(81,1)+r,squeeze(mean(pooled([71, 72],:),1)),'.','markersize',12)
[~,p2]=ttest(squeeze(mean(pooled([71, 72],:),1)))

plot(4*ones(81,1)+r,squeeze(mean(pooled([159, 160],:),1)),'.','markersize',12)
[~,p3]=ttest(squeeze(mean(pooled([159, 160],:),1)))

plot(5*ones(81,1)+r,squeeze(mean(pooled([179, 180],:),1)),'.','markersize',12)
[~,p4]=ttest(squeeze(mean(pooled([179, 180],:),1)))

axis([0.5 5.5 -0.80 0.80])

line([0.8 1.2],[mean(mean(pooled(:,:),1),2) mean(mean(pooled(:,:),1),2)],'linewidth',2)
line([1.8 2.2],[mean(mean(pooled([205, 206],:),1),2) mean(mean(pooled([205, 206],:),1),2)],'linewidth',2)
line([2.8 3.2],[mean(mean(pooled([71, 72],:),1),2) mean(mean(pooled([71, 72],:),1),2)],'linewidth',2)
line([3.8 4.2],[mean(mean(pooled([159, 160],:),1),2) mean(mean(pooled([159, 160],:),1),2)],'linewidth',2)
line([4.8 5.2],[mean(mean(pooled([179, 180],:),1),2) mean(mean(pooled([179, 180],:),1),2)],'linewidth',2)

line([0.5 4.5],[0 0],'linewidth',0.5,'color',[0.7 0.7 0.7],'linestyle','--')

tp_editplots

set(gca,'xtick',1:4,'xticklabels',{'Visual cortex';'Auditory cortex';'Somatosensory cortex';'Anterior cigulate cortex'})
set(gca,'ytick',[-0.8 -0.4 0 0.4 0.8],'yticklabels',[-0.8 -0.4 0 0.4 0.8])
ylabel('Correlation')
% mean(mean(pooled([71, 72],:),1),2)
% mean(mean(pooled([159, 160],:),1),2)
% mean(mean(pooled([179, 180],:),1),2)
 
% 
% for isubj = 1 : size(pooled,3)
%   
%   
%   
%   
% end



%%
figure_w;

subplot(4,4,1)
plot(nanmean(fooof22.slope_gla(idx_sorted,:),2),'color',colors(1,:))
axis([1 246 -0.1 0.02]); tp_editplots; box off
ylabel('Correlation','fontsize',6); xlabel(sprintf('Anterior -> Posterior\n Atlas region'),'fontsize',6)
line([1 246], [0 0],'color',[.7 .7 .7],'linestyle',':')
title('Glasgow')

subplot(4,4,5)
plot(nanmean(fooof22.slope_hh(idx_sorted,:),2),'color',colors(2,:))
axis([1 246 -0.1 0.02]); tp_editplots; box off
xlabel(sprintf('Anterior -> Posterior\n Atlas region'),'fontsize',6)
line([1 246], [0 0],'color',[.7 .7 .7],'linestyle',':')
title('Hamburg')

subplot(4,4,9)
plot(nanmean(fooof22.slope_mue(idx_sorted,:),2),'color',colors(3,:))
axis([1 246 -0.1 0.02]); tp_editplots; box off
xlabel(sprintf('Anterior -> Posterior\n Atlas region'),'fontsize',6)
line([1 246], [0 0],'color',[.7 .7 .7],'linestyle',':')
title('Muenster')

subplot(4,4,13); hold on

pooled = cat(2,fooof22.slope_gla,fooof22.slope_hh,fooof22.slope_mue);
[h,p] = ttest(pooled,zeros(size(pooled)),'dim',2); h = p<fdr1(p(:),0.01,0);
par = nanmean(pooled(idx_sorted,:),2);
plot(par,'k')
plot(find(h(idx_sorted)),par(find(h(idx_sorted))),'r.','markersize',6)
xlabel(sprintf('Anterior -> Posterior\n Atlas region'),'fontsize',6)
axis([1 246 -0.1 0.02]); tp_editplots; box off
line([1 246], [0 0],'color',[.7 .7 .7],'linestyle',':')
title('Pooled')


subplot(4,4,3)
plot(nanmean(fooof11.slope_df_gla(idx_sorted,:),2),'color',colors(1,:))
axis([1 246 -0.07 0.07]); tp_editplots; box off
ylabel('Correlation','fontsize',6); xlabel(sprintf('Anterior -> Posterior\n Atlas region'),'fontsize',6)
line([1 246], [0 0],'color',[.7 .7 .7],'linestyle',':')
title('Glasgow')

subplot(4,4,7)
plot(nanmean(fooof11.slope_df_hh(idx_sorted,:),2),'color',colors(2,:))
axis([1 246 -0.07 0.07]); tp_editplots; box off
xlabel(sprintf('Anterior -> Posterior\n Atlas region'),'fontsize',6)
line([1 246], [0 0],'color',[.7 .7 .7],'linestyle',':')
title('Hamburg')

subplot(4,4,11)
plot(nanmean(fooof11.slope_df_mue(idx_sorted,:),2),'color',colors(3,:))
axis([1 246 -0.07 0.07]); tp_editplots; box off
xlabel(sprintf('Anterior -> Posterior\n Atlas region'),'fontsize',6)
line([1 246], [0 0],'color',[.7 .7 .7],'linestyle',':')
title('Muenster')

subplot(4,4,15); hold on

pooled = cat(2,fooof11.slope_df_gla,fooof11.slope_df_hh,fooof11.slope_df_mue);
[h,p] = ttest(pooled,zeros(size(pooled)),'dim',2); h = p<fdr1(p(:),0.01,0);
par = nanmean(pooled(idx_sorted,:),2);
plot(par,'k')
plot(find(h(idx_sorted)),par(find(h(idx_sorted))),'r.','markersize',6)
xlabel(sprintf('Anterior -> Posterior\n Atlas region'),'fontsize',6)
axis([1 246 -0.07 0.07]); tp_editplots; box off
line([1 246], [0 0],'color',[.7 .7 .7],'linestyle',':')
title('Pooled')

print(gcf,'-dpdf',sprintf('~/pp/plots/figure100.pdf'))


%%

addpath /home/gnolte/meth/highlevel/
load /home/gnolte/meth/templates/sa_template;

load /home/gnolte/meth/templates/mri.mat
figure_w
is_slope = 1; is_dt = 0;


if is_slope && ~is_dt
  pooled = cat(2,fooof22.slope_gla,fooof22.slope_hh,fooof22.slope_mue);
elseif ~is_slope && ~is_dt
  pooled = cat(2,fooof22.offset_gla,fooof22.offset_hh,fooof22.offset_mue);
elseif is_slope && is_dt
  pooled = cat(2,fooof11.slope_df_gla,fooof11.slope_df_hh,fooof11.slope_df_mue);
elseif ~is_slope && is_dt
  pooled = cat(2,fooof11.offset_df_gla,fooof11.offset_df_hh,fooof11.offset_df_mue);
end

cmap = redblue;
par = zeros(8799,81);

for i = 1 : size(pooled,1)
  par(BNA.tissue_5mm==i,:)=repmat(pooled(i,:),[sum(BNA.tissue_5mm==i) 1]); 
end

[h,p]=ttest(par,zeros(size(par)),'dim',2); h = p < fdr1(p(:),0.05,0);

par = nanmean(par,2).*h;
% project onto fine grid
par=spatfiltergauss(par,BNA.grid_5mm./10,.5,sa_template.grid_fine);


clim = [-max([abs([min(par(:)) max(par(:))])]) max([abs([min(par(:)) max(par(:))])])];
para = [];
para.colorlimits = clim
para.colormaps{1} = cmap;
para.orientation = 'axial';

para.dslice_shown = 0.95;
para.colorbar= 0;

tp_showmri_transp(mri,para,[sa_template.grid_fine par]); f=get(gcf)

move_x = [0:0.08:0.32];
move_y = [0.08:-0.02:0];

for i=4:-1:1
  for j = 0:3
    if j ~= 4
      f.Children(i*4-j).Position(1)=f.Children(i*4-j).Position(1)-move_x(j+1);
    end
    f.Children(i*4-j).Position(2)=f.Children(i*4-j).Position(2)+move_y(i);
  end
end

text(1,1,sprintf('[%.3f %.3f]',clim(1),clim(2)))

set(gcf,'renderer','painters')
print(gcf,'-dpdf',sprintf('~/pp/plots/pp_fooof_sourcemap_dt%d_slp%d_v%d.tiff',is_dt,is_slope,v))


%% EXAMPLE FOOOF
close
SUBJLIST = [4 5 6 7 8 9 10 11 12 13 15 16 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34];
isubj = 16;
iblock = 1;
v = 1;

load(sprintf('~/pp/proc/src/pp_gla_collected_fooof_s%d_b%d_v%d.mat',isubj,iblock,v))
load(sprintf('~/pp/proc/src/pp_gla_src_powerspectra_s%d_b%d_v%d.mat',isubj,iblock,v))
 
pup = pup(~isnan(pxx(1,1,:)));
pxx=pxx(:,:,~isnan(pxx(1,1,:)));

[r,i]=min(corr(squeeze(aper(2,:,:))',pup(~isnan(pup))'));

pxx = pxx(:,:,~isnan(pup));
% 
[idx1] = pup>prctile(pup,75);
[idx2] = pup<prctile(pup,25);

% idx1 = 3;
% idx2 = 100;

ap1 = aper(:,i,idx1);
aper1 = (-ap1(2,:)'.*log10(3:0.5:40)+repmat(ap1(1),[75 1])')';
aper1 = mean(aper1,2);

gg1 = squeeze(mean(gg(i,:,idx1),3));

ap2 = aper(:,i,idx2);
aper2 = (-ap2(2,:)'.*log10(3:0.5:40)+repmat(ap2(1),[75 1])')';
aper2 = mean(aper2,2);

figure_w; hold on
% plot(3:0.5:40,gg1,'k'); hold on;
% plot(3:0.5:40,gg2);
% plot(log10(3:0.5:40),aper1,'k:')
plot(log10(3:0.5:40),mean(mean(pxx(3:77,:,idx1),2),3),'k:')
% plot(log10(3:0.5:40),aper2,'r-')
plot(log10(3:0.5:40),mean(mean(pxx(3:77,:,idx2),2),3),'r-')

%% PLOT POWER AFTER TAKING OUT SLOPE
colors = cbrewer('qual', 'Set3', 10,'pchip');
colors = colors(4:6,:);

figure_w;

subplot(4,3,1); hold on
% GLASGOW
% plot(log10(3:0.5:40),nanmean(nanmean(g(:,:,200),1),3))
freqoi = 3:0.5:40;

tmp = squeeze(nanmean(fooof22.psfit_gla,1));
par = nanmean(tmp,2);
par_std = std(tmp,[],2)/sqrt(size(tmp,2));

shadedErrorBar(log10(freqoi),par,par_std,{'color',colors(1,:)})
line([log10(3) log10(40)], [0 0],'color',[.7 .7 .7],'linestyle',':')
axis([log10(3) log10(40) -0.07 0.07]); 
set(gca,'xtick',log10([4 8 16 32]),'xticklabel',[4 8 16 32]);
xlabel('Frequency [Hz]'); tp_editplots


subplot(4,3,2); hold on
% HAMBURG
% plot(log10(3:0.5:40),nanmean(nanmean(g(:,:,200),1),3))
freqoi = 3:0.5:40;

tmp = squeeze(nanmean(fooof22.psfit_hh,1));
par = nanmean(tmp,2);
par_std = std(tmp,[],2)/sqrt(size(tmp,2));

shadedErrorBar(log10(freqoi),par,par_std,{'color',colors(2,:)})
line([log10(3) log10(40)], [0 0],'color',[.7 .7 .7],'linestyle',':')
axis([log10(3) log10(40) -0.07 0.07]); 
set(gca,'xtick',log10([4 8 16 32]),'xticklabel',[4 8 16 32]);
xlabel('Frequency [Hz]'); tp_editplots

subplot(4,3,3); hold on
% HAMBURG
% plot(log10(3:0.5:40),nanmean(nanmean(g(:,:,200),1),3))
freqoi = 3:0.5:40;

tmp = squeeze(nanmean(fooof22.psfit_mue,1));
par = nanmean(tmp,2);
par_std = std(tmp,[],2)/sqrt(size(tmp,2));

shadedErrorBar(log10(freqoi),par,par_std,{'color',colors(3,:)})
line([log10(3) log10(40)], [0 0],'color',[.7 .7 .7],'linestyle',':')
axis([log10(3) log10(40) -0.07 0.07]); 
set(gca,'xtick',log10([4 8 16 32]),'xticklabel',[4 8 16 32]);
xlabel('Frequency [Hz]'); tp_editplots

subplot(4,3,5); hold on
% HAMBURG
% plot(log10(3:0.5:40),nanmean(nanmean(g(:,:,200),1),3))
freqoi = 3:0.5:40;

tmp = cat(2,squeeze(nanmean(fooof22.psfit_gla,1)),squeeze(nanmean(fooof22.psfit_hh,1)),squeeze(nanmean(fooof22.psfit_mue,1)));
par = nanmean(tmp,2);
par_std = std(tmp,[],2)/sqrt(size(tmp,2));

[c,p]=permutest(tmp,zeros(size(squeeze(tmp))),1,0.01,1000,2); h=[c{p<0.05}];

shadedErrorBar(log10(freqoi),par,par_std,'k')
plot(log10(freqoi(h)),par(h),'k.','markersize',8)

line([log10(3) log10(40)], [0 0],'color',[.7 .7 .7],'linestyle',':')
axis([log10(3) log10(40) -0.07 0.07]); 
set(gca,'xtick',log10([4 8 16 32]),'xticklabel',[4 8 16 32]);
xlabel('Frequency [Hz]'); tp_editplots

subplot(4,3,[8 11]); 
[~,idx_sorted] = sort(BNA.centroids(:,2),'descend');
tmp = cat(3,fooof22.psfit_gla,fooof22.psfit_hh,fooof22.psfit_mue);
[h,p] = ttest(tmp,zeros(size(tmp)),'dim',3); h = p<fdr1(p(:),0.05,1);
imagesc(freqoi,1:246,nanmean(tmp(idx_sorted,:,:),3),[-0.05 0.05]); tp_editplots; colormap(cmap)
% set(gca,'xtick',log10([4 8 16 32]),'xticklabel',[4 8 16 32]);
set(gca,'xtick',[0:5:40],'xticklabel',0:5:40,'fontsize',6)

print(gcf,'-dpdf',sprintf('~/pp/plots/pp_fig100b_v%d.pdf',v))

%% TAKE fitted 1/f distribution out of empirical spectra
freqoi=2.^(1:(1/4):7); ffx = 2:0.5:128;
colors = cbrewer('qual', 'Set3', 10,'pchip'); colors = colors(4:6,:);

for ifreq = 1 : length(freqoi)
  [wavelet,opt]=tp_mkwavelet(freqoi(ifreq),0.5,400,0.8);
  idx = find(ffx>=opt.freq_range(1) & ffx<=opt.freq_range(2));
  ps_hh(:,ifreq,:) = squeeze(nanmean(fooof22.ps_hh_corrected(:,idx,:),2));
  ps_gla(:,ifreq,:) = squeeze(nanmean(fooof22.ps_gla_corrected(:,idx,:),2));
  ps_mue(:,ifreq,:) = squeeze(nanmean(fooof22.ps_mue_corrected(:,idx,:),2));
  ps_hh_df(:,ifreq,:) = squeeze(nanmean(fooof11.ps_hh_df_corrected(:,idx,:),2));
  ps_gla_df(:,ifreq,:) = squeeze(nanmean(fooof11.ps_gla_df_corrected(:,idx,:),2));
  ps_mue_df(:,ifreq,:) = squeeze(nanmean(fooof11.ps_mue_df_corrected(:,idx,:),2));
end

figure_w;

subplot(2,4,1)
imagesc(nanmean(ps_gla(idx_sorted,:,:),3),[-0.05 0.05])
set(gca,'xtick',1:4:25,'xticklabel',round(freqoi(1:4:25))); tp_editplots
colormap(cmap); xlabel('Frequency [Hz]')

subplot(2,4,2)
imagesc(nanmean(ps_hh(idx_sorted,:,:),3),[-0.05 0.05])
set(gca,'xtick',1:4:25,'xticklabel',round(freqoi(1:4:25))); tp_editplots
colormap(cmap); xlabel('Frequency [Hz]')

subplot(2,4,3)
imagesc(nanmean(ps_mue(idx_sorted,:,:),3),[-0.05 0.05])
set(gca,'xtick',1:4:25,'xticklabel',round(freqoi(1:4:25))); tp_editplots
colormap(cmap); xlabel('Frequency [Hz]')

subplot(2,4,4)
tmp = cat(3,ps_gla,ps_hh,ps_mue);
[~,p]=ttest(tmp,zeros(size(tmp)),'dim',3); h = p<fdr1(p(:),0.05,1);
imagesc(nanmean(tmp(idx_sorted,:,:),3).*h,[-0.05 0.05])
set(gca,'xtick',1:4:25,'xticklabel',round(freqoi(1:4:25))); tp_editplots
colormap(cmap); xlabel('Frequency [Hz]')

subplot(2,4,5)
imagesc(nanmean(ps_gla_df(idx_sorted,:,:),3),[-0.05 0.05])
set(gca,'xtick',1:4:25,'xticklabel',round(freqoi(1:4:25))); tp_editplots
colormap(cmap); xlabel('Frequency [Hz]')

subplot(2,4,6)
imagesc(nanmean(ps_hh_df(idx_sorted,:,:),3),[-0.05 0.05])
set(gca,'xtick',1:4:25,'xticklabel',round(freqoi(1:4:25))); tp_editplots
colormap(cmap); xlabel('Frequency [Hz]')

subplot(2,4,7)
imagesc(nanmean(ps_mue_df(idx_sorted,:,:),3),[-0.05 0.05])
set(gca,'xtick',1:4:25,'xticklabel',round(freqoi(1:4:25))); tp_editplots
colormap(cmap); xlabel('Frequency [Hz]')

subplot(2,4,8)
tmp = cat(3,ps_gla_df,ps_hh_df,ps_mue_df);
[~,p]=ttest(tmp,zeros(size(tmp)),'dim',3); h = p<fdr1(p(:),0.05,1);
imagesc(nanmean(tmp(idx_sorted,:,:),3).*h,[-0.05 0.05])
set(gca,'xtick',1:4:25,'xticklabel',round(freqoi(1:4:25))); tp_editplots
colormap(cmap); xlabel('Frequency [Hz]')

print(gcf,'-dpdf',sprintf('~/pp/plots/pp_figure100c.pdf',v))



figure_w;
subplot(4,4,[1]); hold on
par = squeeze(nanmean(ps_gla(idx_sorted,:,:),1));
par_std = std(par,[],2)/sqrt(size(par,2));
[c,p]=permutest(par,zeros(size(squeeze(par))),1,0.01,1000,2); h=[c{p<0.05}];
shadedErrorBar(log10(freqoi),nanmean(par,2),par_std,{'color',colors(1,:)})
plot(log10(freqoi(h)),nanmean(par(h,:),2),'.','markersize',8,'color',colors(1,:))
set(gca,'xtick',log10(freqoi(1:4:25)),'xticklabel',round(freqoi(1:4:25))); tp_editplots
colormap(cmap); xlabel('Frequency [Hz]')
axis([log10(freqoi(1)) log10(40) -0.04 0.04])
line([log10(freqoi(1)) log10(40)], [0 0],'color',[.7 .7 .7],'linestyle',':')

subplot(4,4,[2]); hold on
par = squeeze(nanmean(ps_hh(idx_sorted,:,:),1));
par_std = std(par,[],2)/sqrt(size(par,2));
[c,p]=permutest(par,zeros(size(squeeze(par))),1,0.01,1000,2); h=[c{p<0.05}];
shadedErrorBar(log10(freqoi),nanmean(par,2),par_std,{'color',colors(2,:)})
plot(log10(freqoi(h)),nanmean(par(h,:),2),'.','markersize',8,'color',colors(2,:))
set(gca,'xtick',log10(freqoi(1:4:25)),'xticklabel',round(freqoi(1:4:25))); tp_editplots
colormap(cmap); xlabel('Frequency [Hz]')
axis([log10(freqoi(1)) log10(40) -0.04 0.04]);
line([log10(freqoi(1)) log10(40)], [0 0],'color',[.7 .7 .7],'linestyle',':')

subplot(4,4,[3]); hold on
par = squeeze(nanmean(ps_mue(idx_sorted,:,:),1));
par_std = std(par,[],2)/sqrt(size(par,2));
[c,p]=permutest(par,zeros(size(squeeze(par))),1,0.01,1000,2); h=[c{p<0.05}];
shadedErrorBar(log10(freqoi),nanmean(par,2),par_std,{'color',colors(3,:)})
plot(log10(freqoi(h)),nanmean(par(h,:),2),'.','markersize',8,'color',colors(3,:))
set(gca,'xtick',log10(freqoi(1:4:25)),'xticklabel',round(freqoi(1:4:25))); tp_editplots
colormap(cmap); xlabel('Frequency [Hz]')
axis([log10(freqoi(1)) log10(40) -0.04 0.04])
line([log10(freqoi(1)) log10(40)], [0 0],'color',[.7 .7 .7],'linestyle',':')

subplot(4,4,[4]); hold on
par = squeeze(nanmean(cat(3,ps_gla,ps_hh,ps_mue),1));
par_std = std(par,[],2)/sqrt(size(par,2));
[c,p]=permutest(par,zeros(size(squeeze(par))),1,0.01,1000,2); h=[c{p<0.05}];
shadedErrorBar(log10(freqoi),nanmean(par,2),par_std,{'color','k'})
plot(log10(freqoi(h)),nanmean(par(h,:),2),'.','markersize',8,'color','k')
set(gca,'xtick',log10(freqoi(1:4:25)),'xticklabel',round(freqoi(1:4:25))); tp_editplots
colormap(cmap); xlabel('Frequency [Hz]')
axis([log10(freqoi(1)) log10(40) -0.04 0.04])
line([log10(freqoi(1)) log10(40)], [0 0],'color',[.7 .7 .7],'linestyle',':')

subplot(4,4,[5]); hold on
par = squeeze(nanmean(ps_gla_df(idx_sorted,:,:),1));
par_std = std(par,[],2)/sqrt(size(par,2));
[c,p]=permutest(par,zeros(size(squeeze(par))),1,0.01,1000,2); h=[c{p<0.05}];
shadedErrorBar(log10(freqoi),nanmean(par,2),par_std,{'color',colors(1,:)})
plot(log10(freqoi(h)),nanmean(par(h,:),2),'.','markersize',8,'color',colors(1,:))
set(gca,'xtick',log10(freqoi(1:4:25)),'xticklabel',round(freqoi(1:4:25))); tp_editplots
colormap(cmap); xlabel('Frequency [Hz]')
axis([log10(freqoi(1)) log10(40) -0.04 0.04])
line([log10(freqoi(1)) log10(40)], [0 0],'color',[.7 .7 .7],'linestyle',':')

subplot(4,4,[6]); hold on
par = squeeze(nanmean(ps_hh_df(idx_sorted,:,:),1));
par_std = std(par,[],2)/sqrt(size(par,2));
[c,p]=permutest(par,zeros(size(squeeze(par))),1,0.05,1000,2); h=[c{p<0.05}];
shadedErrorBar(log10(freqoi),nanmean(par,2),par_std,{'color',colors(2,:)})
plot(log10(freqoi(h)),nanmean(par(h,:),2),'.','markersize',8,'color',colors(2,:))
set(gca,'xtick',log10(freqoi(1:4:25)),'xticklabel',round(freqoi(1:4:25))); tp_editplots
colormap(cmap); xlabel('Frequency [Hz]')
axis([log10(freqoi(1)) log10(40) -0.04 0.04]);
line([log10(freqoi(1)) log10(40)], [0 0],'color',[.7 .7 .7],'linestyle',':')

subplot(4,4,[7]); hold on
par = squeeze(nanmean(ps_mue_df(idx_sorted,:,:),1));
par_std = std(par,[],2)/sqrt(size(par,2));
[c,p]=permutest(par,zeros(size(squeeze(par))),1,0.05,1000,2); h=[c{p<0.05}];
shadedErrorBar(log10(freqoi),nanmean(par,2),par_std,{'color',colors(3,:)})
plot(log10(freqoi(h)),nanmean(par(h,:),2),'.','markersize',8,'color',colors(3,:))
set(gca,'xtick',log10(freqoi(1:4:25)),'xticklabel',round(freqoi(1:4:25))); tp_editplots
colormap(cmap); xlabel('Frequency [Hz]')
axis([log10(freqoi(1)) log10(40) -0.04 0.04])
line([log10(freqoi(1)) log10(40)], [0 0],'color',[.7 .7 .7],'linestyle',':')

subplot(4,4,[8]); hold on
par = squeeze(nanmean(cat(3,ps_gla_df,ps_hh_df,ps_mue),1));
par_std = std(par,[],2)/sqrt(size(par,2));
[c,p]=permutest(par,zeros(size(squeeze(par))),1,0.01,1000,2); h=[c{p<0.05}];
shadedErrorBar(log10(freqoi),nanmean(par,2),par_std,{'color','k'})
plot(log10(freqoi(h)),nanmean(par(h,:),2),'.','markersize',8,'color','k')
set(gca,'xtick',log10(freqoi(1:4:25)),'xticklabel',round(freqoi(1:4:25))); tp_editplots
colormap(cmap); xlabel('Frequency [Hz]')
axis([log10(freqoi(1)) log10(40) -0.04 0.04])
line([log10(freqoi(1)) log10(40)], [0 0],'color',[.7 .7 .7],'linestyle',':')

print(gcf,'-dpdf',sprintf('~/pp/plots/pp_figure100d.pdf',v))

%%
addpath /home/gnolte/meth/highlevel/
load /home/gnolte/meth/templates/mri.mat
figure_w
is_slope = 1; is_dt = 0;

par = zeros(8799,80);

ifoi = 22;

pooled = squeeze(cat(3,ps_gla(:,ifoi,:),ps_hh(:,ifoi,:),ps_mue(:,ifoi,:)));


for i = 1 : size(pooled,1)
  par(BNA.tissue_5mm==i,:)=repmat(pooled(i,:),[sum(BNA.tissue_5mm==i) 1]); 
end

[~,p]=ttest(par,zeros(size(par)),'dim',2); h = p < fdr1(p(:),0.05,0);

par = nanmean(par,2).*h;
% project onto fine grid
par=spatfiltergauss(par,BNA.grid_5mm./10,.5,sa_template.grid_fine);


clim = [-max([abs([min(par(:)) max(par(:))])]) max([abs([min(par(:)) max(par(:))])])];
para = [];
para.colorlimits = clim
para.colormaps{1} = cmap;
para.orientation = 'axial';

para.dslice_shown = 0.95;
para.colorbar= 0;

tp_showmri_transp(mri,para,[sa_template.grid_fine par]); f=get(gcf)

move_x = [0:0.08:0.32];
move_y = [0.08:-0.02:0];

for i=4:-1:1
  for j = 0:3
    if j ~= 4
      f.Children(i*4-j).Position(1)=f.Children(i*4-j).Position(1)-move_x(j+1);
    end
    f.Children(i*4-j).Position(2)=f.Children(i*4-j).Position(2)+move_y(i);
  end
end

text(1,1,sprintf('[%.3f %.3f]',clim(1),clim(2)))

% set(gcf,'renderer','painters')
% print(gcf,'-dpdf',sprintf('~/pp/plots/pp_fooof_sourcemap_dt%d_slp%d_v%d.tiff',is_dt,is_slope,v))

%% PLOT iverted U
is_dt=0; 
if is_dt == 0
    pooled=double(cat(4,double(nanmean(fooof22.pxx_seg_hh,5)),double(fooof22.pxx_seg_gla),double(fooof22.pxx_seg_mue)));
%     pooled_slp=squeeze(cat(3,mean(fooof22.slopes_seg_hh,4),fooof22.slopes_seg_gla,fooof22.slopes_seg_mue));
else
  pooled=cat(3,mean(fooof11.pxx_seg_dt_hh,4),fooof11.pxx_seg_dt_gla,fooof11.pxx_seg_dt_mue);
end
tmp = mean(pooled,4);

idx_low=fooof22.fxx>2 & fooof22.fxx<4;
idx_alpha=fooof22.fxx>8 & fooof22.fxx<16;
idx_high=fooof22.fxx>64 & fooof22.fxx<128;

[~,idx_sorted] = sort(BNA.centroids(:,2),'descend');

p_low=  zeros(size(pooled,1),size(pooled,2),3,size(pooled,4));
for isubj =1  : size(pooled,4)
  isubj
  for iff = 1  : 253
%     iff
  for i = 1 : 246
    [p_low(iff,i,:,isubj),S_low(i,:)]=polyfit(1:20,mean(tmp(iff,i,:),1),2);
  %   [p_alpha(i,:),S_alpha(i,:)]=polyfit(1:20,mean(tmp(idx_alpha,i,:),1),2);
  %   [p_high(i,:),S_high(i,:)]=polyfit(1:20,mean(tmp(idx_high,i,:),1),2);

  end
  end
end

p_low=p_low(:,idx_sorted,:,:);
figure_w

subplot(1,2,1);
imagesc(fooof22.fxx,1:246,mean(p_low(:,:,1,:),4))

[h,p]=ttest(p_low(:,:,1,:),zeros(size(p_low(:,:,1,:))),'dim',4)
subplot(1,2,2);
imagesc(fooof22.fxx,1:246,mean(p_low(:,:,1,:),4).*h)



figure_w
subplot(3,3,1); hold on
plot(mean(tmp(idx_low,:),1),'o','markersize',7,'markeredgecolor','w','markerfacecolor','k')
plot(1:20,p_low(1).*(1:20).^2 + p_low(2).*(1:20) + p_low(3))
tp_editplots; xlabel('Pupil bin'); ylabel('Power [a.u.]')

subplot(3,3,2); hold on
plot(mean(tmp(idx_alpha,:),1),'o','markersize',7,'markeredgecolor','w','markerfacecolor','k')
plot(1:20,p_alpha(1).*(1:20).^2 + p_alpha(2).*(1:20) + p_alpha(3))
tp_editplots; xlabel('Pupil bin'); ylabel('Power [a.u.]')


subplot(3,3,3); hold on
plot(mean(tmp(idx_high,:),1),'o','markersize',7,'markeredgecolor','w','markerfacecolor','k')
plot(1:20,p_high(1).*(1:20).^2 + p_high(2).*(1:20) + p_high(3))
tp_editplots; xlabel('Pupil bin'); ylabel('Power [a.u.]')

[p_slp,S_slp]=polyfit(1:20,mean(pooled_slp,2),2);

subplot(3,3,4); hold on
plot(mean(pooled_slp,2),'o','markersize',7,'markeredgecolor','w','markerfacecolor','k')
plot(1:20,p_slp(1).*(1:20).^2 + p_slp(2).*(1:20) + p_slp(3))
tp_editplots; xlabel('Pupil bin'); ylabel('Power [a.u.]')

print(gcf,'-dpdf',sprintf('~/pp/plots/pp_invU_isdt%d_v%d.pdf',is_dt,v))

%%

tmp = nanmean(fooof22.pxx_seg_hh,4);
for i = 1 : 246
  i
  for j = 1 : 253
    
  [p(:,i,j),S_low]=polyfit(1:20,squeeze(tmp(j,i,:)),2);
end
  end