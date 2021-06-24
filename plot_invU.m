

% Load data computed in powerspectra.m
fooof22 = pp_load_fooof_results(2);
fooof11 = pp_load_fooof_results(1);

tmp22 = cat(4,nanmean(fooof22.pxx_seg_hh,5),fooof22.pxx_seg_gla,fooof22.pxx_seg_mue);
tmp11 = cat(4,nanmean(fooof11.pxx_seg_hh,5),fooof11.pxx_seg_gla,fooof11.pxx_seg_mue);


fxx = fooof22.fxx;
fooof22.fxx = fxx;

% set notch filter range to nan
idx = (fooof22.fxx>=47 & fooof22.fxx<=53) |  (fooof22.fxx>=97 & fooof22.fxx<=103);
tmp22(idx,:,:,:) = nan;
tmp11(idx,:,:,:) = nan;

pooled22 = zeros(size(tmp22),'double');
pooled11 = zeros(size(tmp11),'double');

% interpolate nans (speed up at some point)
for ireg = 1 : size(pooled22,2)
  ireg
  for iseg = 1 : size(pooled22,3)
    for isubj = 1 : size(pooled22,4)
      
      pooled22(:,ireg,iseg,isubj) = fillmissing(tmp22(:,ireg,iseg,isubj),'spline');
      pooled11(:,ireg,iseg,isubj) = fillmissing(tmp11(:,ireg,iseg,isubj),'spline');
      
    end
  end  
end

%%

clear mean_dat22 mean_dat11 pooled_f11

pooled_f22 = zeros(25,size(pooled22,2),size(pooled22,3),size(pooled22,4),'double');
% pooled_f11 = zeros(25,size(pooled11,2),size(pooled11,3),size(pooled11,4),'double');

f = 2.^(1:0.25:7);
for ii = 1 : 25
  [wavelet, outp]= tp_mkwavelet(f(ii),0.5,400,0);
  freq_range(ii,:) = outp.freq_range;
  pooled_f22(ii,:,:,:) = squeeze(mean(pooled22(fxx>=freq_range(ii,1) & fxx<=freq_range(ii,2),:,:,:),1));
  pooled_f11(ii,:,:,:) = squeeze(mean(pooled11(fxx>freq_range(ii,1) & fxx<freq_range(ii,2),:,:,:),1));
end

x = 1 :25;

p_low22 = zeros(size(pooled_f22,1),size(pooled_f22,2),3,size(pooled_f22,4),'double');
p_low11 = zeros(size(pooled_f22,1),size(pooled_f22,2),3,size(pooled_f22,4),'double');

for isubj = 1  : size(pooled_f22,4)
  
%   m_subj(isubj) = mean(mean(mean1/(pooled_f22(:,:,:,isubj),1),2),3);
  isubj
  rand_idx = randperm(20);
  for iff = 1  : 25
    for i = 1 : 246
      
      dat = squeeze(pooled_f22(iff,i,:,isubj))';
      
      % POLYFIT WITH DEG 2 and 
      if any(isnan(dat))
        sum(isnan(dat))
        dat = zscore(fillmissing(dat,'spline')); 
      else
        dat = zscore(dat);  
      end
    
      p_low22(iff,i,:,isubj)=polyfit(x,dat,2);
      mean_dat22(iff,i,:,isubj) = dat;

      % PERMUTE PUPIL SEGMENTS
      dat = squeeze(pooled_f11(iff,i,:,isubj))';
      % POLYFIT WITH DEG 2 and 
      if any(isnan(dat))
        % interpolation is only needed if pupil segments are not ordered by
        % quantiles
        dat = zscore(fillmissing(dat,'spline')); 
      else
        dat = zscore(dat);  
      end
      
      p_low11(iff,i,:,isubj)=polyfit(x,dat,2);
      
      mean_dat11(iff,i,:,isubj) = dat;

      clear dat

    end
  end
end

%%
close

v=2;
conjunction = 0;
for ifoi = 14
 
%   load(sprintf('~/pp/proc/src/pp_sourcemaps_mask_f%d_v%d.mat',ifoi,v))
  
  % [h,p] = ttest(par_slp_all,zeros(size(par_slp_all)),'dim',2); h = p<fdr1(p(:),0.05,1);
%  ifoi=11
% h = p_max(:,ifoi)<0.05;
h = ones(246,1);
idx = find(h);
% 
% tmp_par = log10(R_sq2)./log10(R_sq1) - log10(R_sq2_rnd)./log10(R_sq1_rnd);
% par_src = zeros(8799,1);
% for i = 1:sum(h)
%   par_src(BNA.tissue_5mm == idx(i),:) = repmat(mean(tmp_par(ifoi,idx(i),:),3)',[sum(BNA.tissue_5mm == idx(i)), 1]);
% end
% 
if v==2
  ppp = p_low22(:,:,1,:);%./mean(pooled_f(:,:,:,:),3);
else
  ppp = p_low11(:,:,1,:);%./mean(pooled_f(:,:,:,:),3);
end
% 
% load(sprintf('~/pvals_%d.mat',ifoi),'h')
% mm = h; clear h 
[h,p]=ttest(ppp(ifoi,:,:,:),zeros(size(ppp(ifoi,:,:,:))),'dim',4);
h = p<fdr1(p(:),0.1,0);

masked = 1;
% prc_max = mean(abs(10.5-idx_max),3);
if masked == 1
  prc_max = mean(ppp(ifoi,:,:,:),4).*h;
else
  prc_max = mean(ppp(ifoi,:,:,:),4);
end


par = zeros(8799,1);
hh=zeros(8799,1);
% mask=zeros(8799,1);
for i = 1:246
  par(BNA.tissue_5mm == idx(i),:) = repmat(prc_max(:,idx(i)),[sum(BNA.tissue_5mm == idx(i)), 1]);
  hh(BNA.tissue_5mm == idx(i),:) = repmat(h(idx(i)),[sum(BNA.tissue_5mm == idx(i)), 1]);
%   mask(BNA.tissue_5mm == idx(i),:) = repmat(mm(idx(i)),[sum(BNA.tissue_5mm == idx(i)), 1]);
end

if conjunction==1
  par = zeros(8799,1);
%   par(mask&(~hh))=30; iii = 5;
% par(mask&hh)=30; iii = 1;
  par((~mask)&hh)=30; iii = 3;

end
sum(par==1)
if conjunction==1
  cmap = cbrewer('qual', 'Set1', 5,'pchip');
  cmap = cmap(iii,:);
%   cmap = cmap;
else
  cmap = redblue;
end
para = [];
para.colorlimits = [-max([abs([min(par(:)) max(par(:))])]) max([abs([min(par(:)) max(par(:))])])];
% para.colorlimits = [min(par(:)) max(par(:))-5];
% para.colorlimits = [5.1 8]
tmp=plasma;
para.colormaps{1} = cmap;
para.orientation = 'axial';
para.dslice_shown = 0.95;
para.colorbar= 0;

tp_showmri_transp(mri,para,[BNA.grid_5mm./10 par]); f=get(gcf)
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

text(1,1,sprintf('[%.8f %.3f]\n [%.8f Hz]',para.colorlimits(1),para.colorlimits(2),ifoi))

set(gcf,'renderer','painters')
if ~conjunction
print(gcf,'-dpdf',sprintf('~/pp/plots/pp_quad_masked%d_spatialmap_f%d_v%d.tiff',masked,ifoi,1))
else
  print(gcf,'-dpdf',sprintf('~/pp/plots/pp_quad_spatialmap_conj%d_f%d_v%d.tiff',iii,ifoi,1))
end
end

%%

ifoi = 22;

[h,p]=ttest(squeeze(p_low22(ifoi,:,1,:)),zeros(size(squeeze(p_low22(ifoi,:,1,:)))),'dim',2);
h=p<fdr1(p,0.1,0);
% h = p<0.01;
par = squeeze(mean(p_low22(ifoi,:,1,:),4));
idx = 1:246;
par1 = zeros(8799,1);
for i = 1:246
  par1(BNA.tissue_5mm == idx(i)) = repmat(par(:,idx(i)).*h(idx(i)),[sum(BNA.tissue_5mm == idx(i)), 1]);
end

% [h,p]=ttest(pooled,zeros(size(pooled)),'dim',3);
% par = mean(pooled(:,ifoi,:),3).*(p(:,ifoi)<0.0005);

cmap      = redblue;
para      = [];
para.clim = [-max([abs([min(par1(:)) max(par1(:))])]) max([abs([min(par1(:)) max(par1(:))])])];
para.cmap = cmap;
para.grid = BNA.grid_5mm/10;
para.dd   = 0.5;
para.fn   = sprintf('~/test_ACC.png');
tp_plot_surface(par1,para)

%%
% 
figure_w;

ppp = p_low22(:,:,3,:);%./mean(pooled_f(:,:,:,:),3);

[h,p,~,s]=ttest(ppp,zeros(size(ppp)),'dim',4);
h=p<0.05;
imagesc(s.tstat(:,:)'.*h',[-3 3])
colormap(redblue)
%%

v = 2;
freqs = 2.^(1:0.25:7);

figure_w

ppp = p_low22(:,:,1,1:end);

load ~/pp/proc/pp_atlas_BNA.mat

if v == 2
is_dt = 0; 
is_task = 0;
else
is_dt = 1; 
is_task = 0;
end


figure_w


subplot(4,3,1); hold on; box off

if v==2
  ppp = p_low22(:,:,1,:);
  sig_V1 = squeeze(mean(p_low22(:,[205 206],1,:),2));
  sig_A1 = squeeze(mean(p_low22(:,[71 72],1,:),2));
  sig_M1 = squeeze(mean(p_low22(:,[159 160],1,:),2));
  sig_ACC = squeeze(mean(p_low22(:,[179 180],1,:),2));
elseif v==1
  ppp = p_low11(:,:,1,:);
  sig_V1 = squeeze(mean(p_low11(:,[205 206],1,:),2));
  sig_A1 = squeeze(mean(p_low11(:,[71 72],1,:),2));
  sig_M1 = squeeze(mean(p_low11(:,[159 160],1,:),2));
  sig_ACC = squeeze(mean(p_low11(:,[179 180],1,:),2));
end

par = squeeze(mean(ppp,2));
[h,p]=ttest(par,zeros(size(par)),'dim',2)
p_fdr = fdr1(p(:),0.1,0);

shadedErrorBar(log10(freqs),mean(par,2),std(par,[],2)/sqrt(size(par,2)))
plot(log10(freqs(find(p<0.05))),mean(par(find(p<0.05),:),2),'k.','markersize',8)
plot(log10(freqs(find(p<0.01))),mean(par(find(p<0.01),:),2),'r.','markersize',8)
plot(log10(freqs(find(p<p_fdr))),mean(par(find(p<p_fdr),:),2),'b.','markersize',8)

tp_editplots
line([log10(freqs(1)) log10(freqs(end))], [0 0])
axis([.3 2.11 -0.013 0.013])
set(gca,'xtick',log10(freqs(1:4:end)),'xticklabels',freqs(1:4:end))


subplot(4,3,2); hold on; box off
[h,p]=ttest(sig_V1,zeros(size(sig_V1)),'dim',2);
p_fdr = fdr1(p(:),0.1,0);

shadedErrorBar(log10(freqs),mean(sig_V1,2),std(par,[],2)/sqrt(size(sig_V1,2)))
plot(log10(freqs(find(p<0.05))),mean(sig_V1(find(p<0.05),:),2),'k.','markersize',8)
plot(log10(freqs(find(p<0.01))),mean(sig_V1(find(p<0.01),:),2),'r.','markersize',8)
plot(log10(freqs(find(p<p_fdr))),mean(sig_V1(find(p<p_fdr),:),2),'r.','markersize',8)

tp_editplots
line([log10(freqs(1)) log10(freqs(end))], [0 0])
set(gca,'xtick',log10(freqs(1:4:end)),'xticklabels',freqs(1:4:end))
axis([.3 2.11 -0.013 0.013])

subplot(4,3,5); hold on; box off
[h,p]=ttest(sig_A1,zeros(size(sig_A1)),'dim',2);
p_fdr = fdr1(p(:),0.1,0);

shadedErrorBar(log10(freqs),mean(sig_A1,2),std(par,[],2)/sqrt(size(sig_A1,2)))
plot(log10(freqs(find(p<0.05))),mean(sig_A1(find(p<0.05),:),2),'k.','markersize',8)
plot(log10(freqs(find(p<0.01))),mean(sig_A1(find(p<0.01),:),2),'r.','markersize',8)
plot(log10(freqs(find(p<p_fdr))),mean(sig_A1(find(p<p_fdr),:),2),'b.','markersize',8)

tp_editplots
line([log10(freqs(1)) log10(freqs(end))], [0 0])
set(gca,'xtick',log10(freqs(1:4:end)),'xticklabels',freqs(1:4:end))
axis([.3 2.11 -0.013 0.013])

subplot(4,3,8); hold on; box off
[h,p]=ttest(sig_M1,zeros(size(sig_M1)),'dim',2);
p_fdr = fdr1(p(:),0.1,0);

shadedErrorBar(log10(freqs),mean(sig_M1,2),std(par,[],2)/sqrt(size(sig_M1,2)))
plot(log10(freqs(find(p<0.05))),mean(sig_M1(find(p<0.05),:),2),'k.','markersize',8)
plot(log10(freqs(find(p<0.01))),mean(sig_M1(find(p<0.01),:),2),'r.','markersize',8)
plot(log10(freqs(find(p<p_fdr))),mean(sig_M1(find(p<p_fdr),:),2),'b.','markersize',8)

tp_editplots
line([log10(freqs(1)) log10(freqs(end))], [0 0])
set(gca,'xtick',log10(freqs(1:4:end)),'xticklabels',freqs(1:4:end))
axis([.3 2.11 -0.013 0.013])

subplot(4,3,11);hold on; box off
[h,p]=ttest(sig_ACC,zeros(size(sig_ACC)),'dim',2);
p_fdr = fdr1(p(:),0.1,0);

shadedErrorBar(log10(freqs),mean(sig_ACC,2),std(par,[],2)/sqrt(size(sig_ACC,2)))
plot(log10(freqs(find(p<0.05))),mean(sig_ACC(find(p<0.05),:),2),'k.','markersize',8)
plot(log10(freqs(find(p<0.01))),mean(sig_ACC(find(p<0.01),:),2),'r.','markersize',8)
plot(log10(freqs(find(p<p_fdr))),mean(sig_ACC(find(p<p_fdr),:),2),'b.','markersize',8)

tp_editplots
line([log10(freqs(1)) log10(freqs(end))], [0 0])
set(gca,'xtick',log10(freqs(1:4:end)),'xticklabels',freqs(1:4:end))
axis([.3 2.11 -0.013 0.013])
% close al
print(gcf,'-dpdf',sprintf('~/pp/plots/pp_src_ROIs_invertedU_dt%d_v%d.pdf',0,v))


%%
clear m22  poly_par22 m11 poly_par11
% is_dt=0; 
% pooled11=double(cat(4,double(nanmean(fooof11.pxx_seg_dt_hh,5)),double(fooof11.pxx_seg_dt_gla),double(fooof11.pxx_seg_dt_mue)));

% tmp22 = squeeze(mean(mean(pooled22,4),2));
% tmp11 = squeeze(mean(mean(pooled11,4),2));

for i=1:81
  i
  for ii=1:246
    fxx = 2:0.5:128;
    idx_low=freqs>=2 & freqs<=4;
    idx_alpha=freqs>=8 & freqs<=16;
    idx_high=freqs>=64 & freqs<=128;
    
    m22(1,ii,i,:)=mean(mean_dat22(idx_low,ii,:,i),1);
    m22(2,ii,i,:)=mean(mean_dat22(idx_alpha,ii,:,i),1);
    m22(3,ii,i,:)=mean(mean_dat22(idx_high,ii,:,i),1);

    m11(1,ii,i,:)=mean(mean_dat11(idx_low,ii,:,i),1);
    m11(2,ii,i,:)=mean(mean_dat11(idx_alpha,ii,:,i),1);
    m11(3,ii,i,:)=mean(mean_dat11(idx_high,ii,:,i),1);
    
  end
end
%
poly_par22(1,:) = polyfit(1:size(m22,4),squeeze(mean(mean(m22(1,:,:,:),2),3)),2);
poly_par22(2,:) = polyfit(1:size(m22,4),squeeze(mean(mean(m22(2,:,:,:),2),3)),2);
poly_par22(3,:) = polyfit(1:size(m22,4),squeeze(mean(mean(m22(3,:,:,:),2),3)),2);
  
poly_par11(1,:) = polyfit(1:size(m11,4),squeeze(mean(mean(m11(1,:,:,:),2),3)),2);
poly_par11(2,:) = polyfit(1:size(m11,4),squeeze(mean(mean(m11(2,:,:,:),2),3)),2);
poly_par11(3,:) = polyfit(1:size(m11,4),squeeze(mean(mean(m11(3,:,:,:),2),3)),2);
 
figure_w
subplot(3,3,1); hold on
plot(squeeze(mean(mean(m22(1,:,:,:),2),3)),'o','markersize',7,'markeredgecolor','w','markerfacecolor','k')
plot(1:size(m22,4),poly_par22(1,1).*(1:size(m22,4)).^2 + poly_par22(1,2).*(1:size(m22,4)) + poly_par22(1,3))
tp_editplots; xlabel('Pupil bin'); ylabel('Power [norm.]')
xlim([-1 size(m22,4)+1])
% ylim([0.9 1.1])

subplot(3,3,2); hold on
plot(squeeze(mean(mean(m22(2,:,:,:),2),3)),'o','markersize',7,'markeredgecolor','w','markerfacecolor','k')
plot(1:size(m22,4),poly_par22(2,1).*(1:size(m22,4)).^2 + poly_par22(2,2).*(1:size(m22,4)) + poly_par22(2,3))
tp_editplots; xlabel('Pupil bin'); ylabel('Power [norm.]')
xlim([-1 size(m22,4)+1])
% ylim([0.9 1.1])

subplot(3,3,3); hold on
plot(squeeze(mean(mean(m22(3,:,:,:),2),3)),'o','markersize',7,'markeredgecolor','w','markerfacecolor','k')
plot(1:size(m22,4),poly_par22(3,1).*(1:size(m22,4)).^2 + poly_par22(3,2).*(1:size(m22,4)) + poly_par22(3,3))
tp_editplots; xlabel('Pupil bin'); ylabel('Power [norm.]')
xlim([-1 size(m22,4)+1])
% ylim([0.95 1.03])

% 
subplot(3,3,4); hold on
plot(squeeze(mean(mean(m11(1,:,:,:),2),3)),'o','markersize',7,'markeredgecolor','w','markerfacecolor','k')
plot(1:25,poly_par11(1,1).*(1:25).^2 + poly_par11(1,2).*(1:25) + poly_par11(1,3))
tp_editplots; xlabel('Pupil bin'); ylabel('Power [norm.]')
xlim([-1 size(m22,4)+1])
% ylim([0.9 1.1])

subplot(3,3,5); hold on
plot(squeeze(mean(mean(m11(2,:,:,:),2),3)),'o','markersize',7,'markeredgecolor','w','markerfacecolor','k')
plot(1:25,poly_par11(2,1).*(1:25).^2 + poly_par11(2,2).*(1:25) + poly_par11(2,3))
tp_editplots; xlabel('Pupil bin'); ylabel('Power [norm.]')
xlim([-1 size(m22,4)+1])
% ylim([0.9 1.1])

subplot(3,3,6); hold on
plot(squeeze(mean(mean(m11(3,:,:,:),2),3)),'o','markersize',7,'markeredgecolor','w','markerfacecolor','k')
plot(1:25,poly_par11(3,1).*(1:25).^2 + poly_par11(3,2).*(1:25) + poly_par11(3,3))
tp_editplots; xlabel('Pupil bin'); ylabel('Power [norm.]')
xlim([-1 size(m22,4)+1])
% ylim([0.9 1.1])

print(gcf,'-dpdf',sprintf('~/pp/plots/pp_invU_avg_v%d.pdf',v))