%% PLOT RESULTS
addpath ~/Documents/MATLAB/Colormaps/Colormaps' (5)'/Colormaps/
addpath ~/Documents/MATLAB/cbrewer/cbrewer/
SUBJLIST  = [4 5 6 7 8 9 10 11 12 13 15 16 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34];

cmap = cbrewer('div','RdBu',128); cmap = cmap(end:-1:1,:)
%% PUPIL POWER CORREALTIONS IN SOURCE SPACE (BNA ATLAS)
freqoi=2.^(1:(1/4):7); % 1-128 Hz as per Hipp et al. (2012) Nat Neurosci

v = 1;
clear all_corr
for isubj = SUBJLIST
  for iblock = 1 : 2
    clear src_r
    load(sprintf([outdir 'pp_src_pupil_power_correlations_s%d_b%d_v%d.mat'],isubj,iblock,v));
    all_corr(:,:,isubj,1,iblock) = src_r;
    load(sprintf([outdir 'pp_cnt_src_pupil_power_correlations_s%d_b%d_v%d.mat'],isubj,iblock,v));
    all_corr(:,:,isubj,2,iblock) = src_r;
  end
end

all_corr= nanmean(all_corr(:,:,SUBJLIST,:,:),5);

subj = 1:28;

figure; set(gcf,'color','w','Position',[100 100 800 800]); 

subplot(2,3,1)
imagesc(squeeze(nanmean(all_corr(:,:,subj,1),3)),[-0.05 0.05]); 
set(gca,'xtick',[1 5 9 13 17 21 25],'xticklabel',num2cell(freqoi([1 5 9 13 17 21 25])))
xlabel('Carrier frequency [Hz]'); ylabel('BNA region'); title('Fixation'); tp_editplots

subplot(2,3,2)
imagesc(squeeze(nanmean(all_corr(:,:,subj,2),3)),[-0.05 0.05]); 
set(gca,'xtick',[1 5 9 13 17 21 25],'xticklabel',num2cell(freqoi([1 5 9 13 17 21 25])))
xlabel('Carrier frequency [Hz]'); ylabel('BNA region'); title('Counting'); tp_editplots

subplot(2,3,3)
imagesc((nanmean(all_corr(:,:,subj,2),3)-nanmean(all_corr(:,:,subj,1),3)),[-0.05 0.05]);
set(gca,'xtick',[1 5 9 13 17 21 25],'xticklabel',num2cell(freqoi([1 5 9 13 17 21 25])))
xlabel('Carrier frequency [Hz]'); ylabel('BNA region'); title('Count. - Fix.'); tp_editplots

subplot(2,3,4)
[~,p]=ttest(zeros(size(all_corr(:,:,:,1))),all_corr(:,:,:,1),'dim',3); h = p < fdr1(p(:),0.05,1);
imagesc(squeeze(nanmean(all_corr(:,:,:,1),3)).*h,[-0.05 0.05]); 
set(gca,'xtick',[1 5 9 13 17 21 25],'xticklabel',num2cell(freqoi([1 5 9 13 17 21 25])))
xlabel('Carrier frequency [Hz]'); ylabel('BNA region'); title('Fixation'); tp_editplots

subplot(2,3,5)
[~,p]=ttest(zeros(size(all_corr(:,:,:,1))),all_corr(:,:,:,2),'dim',3); h = p < fdr1(p(:),0.05,1);
imagesc(squeeze(nanmean(all_corr(:,:,:,2),3)).*h,[-0.05 0.05]); 
set(gca,'xtick',[1 5 9 13 17 21 25],'xticklabel',num2cell(freqoi([1 5 9 13 17 21 25])))
xlabel('Carrier frequency [Hz]'); ylabel('BNA region'); title('Counting'); tp_editplots

subplot(2,3,6)
[~,p]=ttest(all_corr(:,:,:,2),all_corr(:,:,:,1),'dim',3); h = p < fdr1(p(:),0.05,1);
imagesc((nanmean(all_corr(:,:,:,2),3)-nanmean(all_corr(:,:,:,1),3)).*h,[-0.05 0.05]);
set(gca,'xtick',[1 5 9 13 17 21 25],'xticklabel',num2cell(freqoi([1 5 9 13 17 21 25])))
xlabel('Carrier frequency [Hz]'); ylabel('BNA region'); title('Count. - Fix.'); tp_editplots


colormap(cmap)
%% CHECK BELOW. MADE CHANGES HERE.

para.time = [3000 7000];
pup = pp_loadpupil(para);


%%
SUBJLIST1  = [4 5 6 7 8 9 10 11 12 13 15 16 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34];
para.str_behav = 'count';
behav = pconn_read_behavioral_data(SUBJLIST1,para);
behav_cnt = behav;
% behav_cnt(27,:) = [];

% delta_cnt = 100*(behav_cnt(:,2)-behav_cnt(:,1))./behav_cnt(:,1);


para.str_behav = 'numb_switches';
behav = pconn_read_behavioral_data(SUBJLIST1,para);
behav_bttn = behav;
behav_bttn = permute(behav_bttn,[2 1 3]);
% behav_bttn(27,:) = []

% behav=   (behav_cnt+behav_bttn)./2;
% behav(isnan(behav_bttn))=behav_cnt(isnan(behav_bttn));
% 
% delta_behav = 100*(behav(:,2)-behav(:,1))./behav(:,1);

%%
ord = pconn_randomization;
for isubj = SUBJLIST
  for m = 1 : 3
    im = find(ord(isubj,:)==m);
    for iblock = 1 : 2
      try
      load(sprintf(['~/pconn_bttn/proc/' 'pconn_bttn_trig_samples_s%d_m%d_b%d_v%d.mat'],isubj,im,iblock,1));
      
      num(isubj,m,iblock) =  length(trig);
%       [tmp, lags]=autocorr(diff(trig));
%       acf(:,isubj,m,iblock)= tmp(1:10); clear tmp
      catch me
        num(isubj,m,iblock) =  nan;
%         num_std(isubj,m,iblock)
% = nan;
%         acf(:,isubj,m,iblock)=nan(1,10);
      end
      
    end
  end
end

num = nanmean(num(SUBJLIST,:,:),3);
% acf = nanmean(acf(:,SUBJLIST,:,:),4);

%% PUPIL CONNECTIVITY
outdir = '~/pp/proc/conn/';

idx = 1:46;
idx = ~ismember(idx,[21 22 23 44 45 46]);
clear pc powcorr
ord = pconn_randomization;

for isubj= SUBJLIST
  isubj
  for m = 1 : 3
    im = find(ord(isubj,:)==m);
    for iblock = 1 : 2
      for ifoi = 1 : 13
        
        try
          load(sprintf([outdir 'pp_src_powcorr_s%d_m%d_b%d_f%d_v%d.mat'],isubj,im,iblock,ifoi,20))
          pc(:,:,isubj,m,1,ifoi,iblock) = powcorr;
        catch me
          pc(:,:,isubj,m,1,ifoi,iblock) = nan(46,46);
        end
        
        try
          load(sprintf([outdir 'pp_task_src_powcorr_s%d_m%d_b%d_f%d_v%d.mat'],isubj,im,iblock,ifoi,20))
          pc(:,:,isubj,m,2,ifoi,iblock) = powcorr;
        catch me
          pc(:,:,isubj,m,2,ifoi,iblock) = nan(46,46);
        end
        
        
        clear powcorr
      end
      
    end
  end
end

pc_mean   = squeeze(nanmean(pc(idx,idx,SUBJLIST,:,:,:,:),7));
pc_blocks = squeeze(pc(idx,idx,SUBJLIST,:,:,:,:));

%%
% tmp = tp_create_grid('vtpm');
clim = [-0.025 0.025];
reg = tmp.tissuelabel_4mm(idx);
% [h,p] = ttest(pc(:,:,:,2),pc(:,:,:,1),'dim',3,'alpha',0.01);
% h = p<fdr1(p(:),0.025);

d = nanmean(pc_mean(:,:,:,2)-pc_mean(:,:,:,1),3);
tp_plot_matrix(d,'clim',clim)
set(gca,'yTick',1:2:40,'yTickLabels',reg(1:2:40),'ticklabelinterpreter','none')
set(gca,'xTick',1:2:40,'xTickLabels',reg(1:2:40),'ticklabelinterpreter','none')
xtickangle(90); axis square
% colormap(cmap)

%%
figure; set(gcf,'color','w');
% behav_cnt = nanmean(behav_cnt,3)';
mask = logical(tril(ones(40,40),-1));

d_behav = behav_cnt(:,2)-behav_cnt(:,1);
d_fc = pc_mean(:,:,:,2)-pc_mean(:,:,:,1);

for i = 1 :40
  for j = 1 : 40
    if i == j; r(i,j)=nan; continue; end
    
    [r(i,j) p(i,j)] = corr(squeeze(d_fc(i,j,:)),d_behav);
    
  end
end

subplot(2,2,1); imagesc(r,[-0.3 0.3]); axis square; title('Correlation w behav')
set(gca,'xTick',1:2:40,'xTickLabels',reg(1:2:40),'ticklabelinterpreter','none');xtickangle(90)
set(gca,'yTick',1:2:40,'yTickLabels',reg(1:2:40),'ticklabelinterpreter','none')
tp_editplots
subplot(2,2,2); imagesc(r.*(p<0.05),[-0.3 0.3]); axis square;
set(gca,'xTick',1:2:40,'xTickLabels',reg(1:2:40),'ticklabelinterpreter','none');xtickangle(90)
set(gca,'yTick',1:2:40,'yTickLabels',reg(1:2:40),'ticklabelinterpreter','none')
colormap(cmap); tp_editplots

%%
% pup = nanmean(pup(:,:,:,2),3);
% figure; set(gcf,'color','w');

for i = 1 :40
  for j = 1 : 40
    if i == j; r(i,j)=nan; continue; end
    
    [r(i,j) p(i,j)] = corr(squeeze(d_fc(i,j,:)),pup(:,2)-pup(:,1));
    
  end
end

subplot(2,2,3); imagesc(r,[-0.3 0.3]); axis square; title('Correlation w pupil')
set(gca,'xTick',1:2:40,'xTickLabels',reg(1:2:40),'ticklabelinterpreter','none');xtickangle(90)
set(gca,'yTick',1:2:40,'yTickLabels',reg(1:2:40),'ticklabelinterpreter','none')
tp_editplots
subplot(2,2,4); imagesc(r.*(p<0.05),[-0.3 0.3]); axis square;
set(gca,'xTick',1:2:40,'xTickLabels',reg(1:2:40),'ticklabelinterpreter','none');xtickangle(90)
set(gca,'yTick',1:2:40,'yTickLabels',reg(1:2:40),'ticklabelinterpreter','none')
colormap(cmap); tp_editplots

%% CONSISTENCY FC MATRICES

for isubj = 1 : 28
  tmp = squeeze(pc_blocks(:,:,isubj,2,:))-squeeze(pc_blocks(:,:,isubj,1,:));
  tmp1 = tmp(:,:,1); tmp1 = tmp1(mask);
  tmp2 = tmp(:,:,2); tmp2 = tmp2(mask);
  [r(isubj) p(isubj)] = corr(tmp1,tmp2);
end

%% PLOT VTPM FC - BEHAVIOR CORRELATIONS FOR BOTH BLOCKS

ipharm = 2;

d_pc1 = pc_blocks(:,:,:,ipharm,1)-pc_blocks(:,:,:,1,1);
d_pc2 = pc_blocks(:,:,:,ipharm,2)-pc_blocks(:,:,:,1,2);

d_beh1 = squeeze(behav_cnt(ipharm,:,1)- behav_cnt(1,:,1))';
d_beh2 = squeeze(behav_cnt(ipharm,:,2)- behav_cnt(1,:,2))';
% d_beh1 = squeeze(behav_bttn(2,:,1)- behav_bttn(1,:,1))';
% d_beh2 = squeeze(behav_bttn(2,:,2)- behav_bttn(1,:,2))';

for i = 1 : 40
  for j = 1 : 40
    if i==j; r_pc_behav(i,j,:)=nan(2,1); continue; end
     idx_nan1 = any([isnan(squeeze(d_pc1(i,j,:))) isnan(d_beh1)],2);
     [r_pc_behav(i,j,1) p_pc_behav(i,j,1)] = corr(squeeze(d_pc1(i,j,~idx_nan1)),d_beh1(~idx_nan1));
     idx_nan2 = any([isnan(squeeze(d_pc2(i,j,:))) isnan(d_beh2)],2);
     [r_pc_behav(i,j,2) p_pc_behav(i,j,2)] = corr(squeeze(d_pc2(i,j,~idx_nan2)),d_beh2(~idx_nan2));
  end
end
  
figure; set(gcf,'color','w');
subplot(2,2,1); imagesc(r_pc_behav(:,:,1),[-0.3 0.3]); axis square; title('Correlation w behav: B1')
set(gca,'xTick',1:2:40,'xTickLabels',reg(1:2:40),'ticklabelinterpreter','none');xtickangle(90)
set(gca,'yTick',1:2:40,'yTickLabels',reg(1:2:40),'ticklabelinterpreter','none')
tp_editplots; colorbar
subplot(2,2,2); imagesc(r_pc_behav(:,:,1).*(p_pc_behav(:,:,1)<0.05),[-0.3 0.3]); axis square;
set(gca,'xTick',1:2:40,'xTickLabels',reg(1:2:40),'ticklabelinterpreter','none');xtickangle(90)
set(gca,'yTick',1:2:40,'yTickLabels',reg(1:2:40),'ticklabelinterpreter','none')
colormap(cmap); tp_editplots

subplot(2,2,3); imagesc(r_pc_behav(:,:,2),[-0.3 0.3]); axis square; title('Correlation w behav: B2')
set(gca,'xTick',1:2:40,'xTickLabels',reg(1:2:40),'ticklabelinterpreter','none');xtickangle(90)
set(gca,'yTick',1:2:40,'yTickLabels',reg(1:2:40),'ticklabelinterpreter','none')
tp_editplots; colorbar
subplot(2,2,4); imagesc(r_pc_behav(:,:,2).*(p_pc_behav(:,:,2)<0.05),[-0.3 0.3]); axis square;
set(gca,'xTick',1:2:40,'xTickLabels',reg(1:2:40),'ticklabelinterpreter','none');xtickangle(90)
set(gca,'yTick',1:2:40,'yTickLabels',reg(1:2:40),'ticklabelinterpreter','none')
colormap(cmap); tp_editplots

print(gcf,'-dpdf',sprintf('~/pp/plots/pp_fc_behav.pdf'))

figure; set(gcf,'color','w');
subplot(2,2,1);
scatter(reshape(r_pc_behav(:,:,1),[40*40 1]),reshape(r_pc_behav(:,:,2),[40*40 1]),'markerfacecolor','k','markeredgecolor','w')
axis([-0.5 0.5 -0.5 0.5])
line([-0.5 0.5],[0 0],'color',[0.8 0.8 0.8],'linestyle',':')
line([0 0],[-0.5 0.5],'color',[0.8 0.8 0.8],'linestyle',':')
line([-0.5 0.5],[-0.5 0.5],'color',[0 0 0],'linestyle','-')
tp_editplots; xlabel('Correlation (Block #1)'); ylabel('Correlation (Block #2)')

%% PLOT VTPM PUPIL - BEHAVIOR CORRELATIONS FOR BOTH BLOCKS

para.time = [3000 7000];
if any(size(pup)~=[28 3 2 3])
  pup = pp_loadpupil(para);
end

d_pup1 = pup(:,ipharm,1,2)-pup(:,1,1,2); 
d_pup2 = pup(:,ipharm,2,2)-pup(:,1,2,2); 

for i = 1 : 40
  for j = 1 : 40
    if i==j; r_pc_pup(i,j,:)=nan(2,1); p_pc_pup(i,j,:)=nan(2,1); continue; end
     idx_nan1 = any([isnan(squeeze(d_pc1(i,j,:))) isnan(d_pup1)],2);
     [r_pc_pup(i,j,1), p_pc_pup(i,j,1)] = corr(squeeze(d_pc1(i,j,~idx_nan1)),d_pup1(~idx_nan1));
     idx_nan2 = any([isnan(squeeze(d_pc2(i,j,:))) isnan(d_pup2)],2);
     [r_pc_pup(i,j,2), p_pc_pup(i,j,2)] = corr(squeeze(d_pc2(i,j,~idx_nan2)),d_pup2(~idx_nan2));
  end
end
  
figure; set(gcf,'color','w'); set(gcf,'Position',[100 100 1200 600])

subplot(2,4,1); imagesc(nanmean(d_pc1,3),[-0.02 0.02]); axis square; title('\DeltaFC: B1')
set(gca,'xTick',1:2:40,'xTickLabels',reg(1:2:40),'ticklabelinterpreter','none');xtickangle(90)
set(gca,'yTick',1:2:40,'yTickLabels',reg(1:2:40),'ticklabelinterpreter','none')
tp_editplots; 
subplot(2,4,5); imagesc(nanmean(d_pc1,3).*(ttest(d_pc1,zeros(size(d_pc1)),'dim',3)),[-0.02 0.02]); axis square;
set(gca,'xTick',1:2:40,'xTickLabels',reg(1:2:40),'ticklabelinterpreter','none');xtickangle(90)
set(gca,'yTick',1:2:40,'yTickLabels',reg(1:2:40),'ticklabelinterpreter','none')
colormap(cmap); tp_editplots

subplot(2,4,2); imagesc(nanmean(d_pc2,3),[-0.02 0.02]); axis square; title('\DeltaFC: B2')
set(gca,'xTick',1:2:40,'xTickLabels',reg(1:2:40),'ticklabelinterpreter','none');xtickangle(90)
set(gca,'yTick',1:2:40,'yTickLabels',reg(1:2:40),'ticklabelinterpreter','none')
tp_editplots; 
subplot(2,4,6); imagesc(nanmean(d_pc2,3).*(ttest(d_pc2,zeros(size(d_pc2)),'dim',3)),[-0.02 0.02]); axis square;
set(gca,'xTick',1:2:40,'xTickLabels',reg(1:2:40),'ticklabelinterpreter','none');xtickangle(90)
set(gca,'yTick',1:2:40,'yTickLabels',reg(1:2:40),'ticklabelinterpreter','none')
colormap(cmap); tp_editplots

subplot(2,4,3); imagesc(r_pc_pup(:,:,1),[-0.3 0.3]); axis square; title('Correlation w pupil: B1')
set(gca,'xTick',1:2:40,'xTickLabels',reg(1:2:40),'ticklabelinterpreter','none');xtickangle(90)
set(gca,'yTick',1:2:40,'yTickLabels',reg(1:2:40),'ticklabelinterpreter','none')
tp_editplots; 
subplot(2,4,7); imagesc(r_pc_pup(:,:,1).*(p_pc_pup(:,:,1)<0.05),[-0.3 0.3]); axis square;
set(gca,'xTick',1:2:40,'xTickLabels',reg(1:2:40),'ticklabelinterpreter','none');xtickangle(90)
set(gca,'yTick',1:2:40,'yTickLabels',reg(1:2:40),'ticklabelinterpreter','none')
colormap(cmap); tp_editplots

subplot(2,4,4); imagesc(r_pc_pup(:,:,2),[-0.3 0.3]); axis square; title('Correlation w pupil: B2')
set(gca,'xTick',1:2:40,'xTickLabels',reg(1:2:40),'ticklabelinterpreter','none');xtickangle(90)
set(gca,'yTick',1:2:40,'yTickLabels',reg(1:2:40),'ticklabelinterpreter','none')
tp_editplots; 
subplot(2,4,8); imagesc(r_pc_pup(:,:,2).*(p_pc_pup(:,:,2)<0.05),[-0.3 0.3]); axis square;
set(gca,'xTick',1:2:40,'xTickLabels',reg(1:2:40),'ticklabelinterpreter','none');xtickangle(90)
set(gca,'yTick',1:2:40,'yTickLabels',reg(1:2:40),'ticklabelinterpreter','none')
colormap(cmap); tp_editplots

print(gcf,'-dpdf',sprintf('~/pp/plots/pp_fc_pupil.pdf'))

figure; set(gcf,'color','w');
subplot(2,2,1);
scatter(reshape(r_pc_pup(:,:,1),[40*40 1]),reshape(r_pc_pup(:,:,2),[40*40 1]),'markerfacecolor','k','markeredgecolor','w')
lsline
axis([-0.5 0.5 -0.5 0.5])
line([-0.5 0.5],[0 0],'color',[0.8 0.8 0.8],'linestyle',':')
line([0 0],[-0.5 0.5],'color',[0.8 0.8 0.8],'linestyle',':')
line([-0.5 0.5],[-0.5 0.5],'color',[0 0 0],'linestyle','-')
tp_editplots; xlabel('Correlation (Block #1)'); ylabel('Correlation (Block #2)')

par1 = reshape(r_pc_pup(:,:,1),[40*40 1]);
par2 = reshape(r_pc_pup(:,:,2),[40*40 1]);
nan_idx = isnan(par1);
[r,p] = corr(par1(~isnan(par1)),par2(~isnan(par2)),'type','spearman');

%%  

ipharm = 1;

for i = 1 : 40
  for j = 1 : 40
    if i==j; r(i,j)=nan; continue; end
    
    par1 = squeeze(pc_mean(i,j,:,ipharm));
    par2 = squeeze(nanmean(pup(:,ipharm,1,2),3)); 
    nan_idx = any([isnan(par1) isnan(par2)],2);
    r(i,j) = corr(par1(~nan_idx),par2(~nan_idx));
  end
end

imagesc(r,[-0.5 0.5]);

ipharm = 1;

for i = 1 : 40
  for j = 1 : 40
    if i==j; r(i,j)=nan; continue; end
    
    par1 = squeeze(pc_mean(i,j,:,ipharm));
    par2 = squeeze(nanmean(pup(:,ipharm,:,2),3));
    
    nan_idx = any([isnan(par1) isnan(par2)],2);
    
    r(i,j) = corr(par1(~nan_idx),par2(~nan_idx));
  end
end

imagesc(r,[-0.5 0.5])

%% PUPIL BLOCKS

[~,i] = sort(squeeze(pup(:,1,:,2)),2);

for isubj = 1 : 28
  
  new_pc(:,:,isubj,:,[1 2]) = pc_blocks(:,:,isubj,:,i(isubj,:));
  
end

ipharm = 3;
figure; set(gcf,'color','w');
subplot(1,2,1);
[h,p]= ttest(new_pc(:,:,:,ipharm,2),new_pc(:,:,:,ipharm,1),'dim',3);
imagesc(nanmean(new_pc(:,:,:,ipharm,2),3)-nanmean(new_pc(:,:,:,ipharm,1),3),[-0.02 0.02])
set(gca,'xTick',1:1:40,'xTickLabels',reg(1:1:40),'ticklabelinterpreter','none');xtickangle(90)
set(gca,'yTick',1:1:40,'yTickLabels',reg(1:1:40),'ticklabelinterpreter','none')
colormap(cmap); tp_editplots; axis square

subplot(1,2,2);
imagesc((nanmean(new_pc(:,:,:,ipharm,2),3)-nanmean(new_pc(:,:,:,ipharm,1),3)).*h,[-0.02 0.02])
set(gca,'xTick',1:1:40,'xTickLabels',reg(1:1:40),'ticklabelinterpreter','none');xtickangle(90)
set(gca,'yTick',1:1:40,'yTickLabels',reg(1:1:40),'ticklabelinterpreter','none')
colormap(cmap); tp_editplots; axis square
%%
cmap = cbrewer('div','RdBu',128); cmap = cmap(end:-1:1,:)
% fc = pupmod_loadpowcorr(12,0);

figure; set(gcf,'color','w');

subplot(2,2,1); title('Atx: Block 1')
[h,p]=ttest(fc(:,:,:,2,2,6,1),fc(:,:,:,1,2,6,1),'dim',3);
imagesc(nanmean(fc(:,:,:,2,2,6,1),3)-nanmean(fc(:,:,:,1,2,6,1),3),[-0.02 0.02]); axis square

subplot(2,2,2); title('Atx: Block 2')
[h,p]=ttest(fc(:,:,:,2,2,6,1),fc(:,:,:,1,2,6,1),'dim',3);
imagesc(nanmean(fc(:,:,:,2,2,6,2),3)-nanmean(fc(:,:,:,1,2,6,2),3),[-0.02 0.02]); axis square

subplot(2,2,3); title('Dpz: Block 1')
[h,p]=ttest(fc(:,:,:,3,2,6,1),fc(:,:,:,1,2,6,1),'dim',3);
imagesc(nanmean(fc(:,:,:,3,2,6,1),3)-nanmean(fc(:,:,:,1,2,6,1),3),[-0.02 0.02]); axis square

subplot(2,2,4); title('Dpz: Block 2')
[h,p]=ttest(fc(:,:,:,3,2,6,1),fc(:,:,:,1,2,6,1),'dim',3);
imagesc(nanmean(fc(:,:,:,3,2,6,2),3)-nanmean(fc(:,:,:,1,2,6,2),3),[-0.02 0.02]); axis square
colormap(cmap)


%% BEHAVIOR EFFECT for each BLOCK

par_cnt1 = squeeze(behav_cnt(2,:,1)-behav_cnt(1,:,1));
par_cnt2 = squeeze(behav_cnt(2,:,2)-behav_cnt(1,:,2));

par_bttn1 = squeeze(behav_bttn(2,:,1)-behav_bttn(1,:,1));
par_bttn2 = squeeze(behav_bttn(2,:,2)-behav_bttn(1,:,2));

figure; set(gcf,'color','w');
subplot(2,2,1);

m1 = nanmean(par_cnt1);
s1 = nanstd(par_cnt1)./length(par_cnt1);
m2 = nanmean(par_cnt2);
s2 = nanstd(par_cnt2)./length(par_cnt2);

bar([m1 m2]); axis([0 3 -15 15])
set(gca,'xtick',[1 2],'xticklabels',{'Block1';'Block2'}); xtickangle(90);
tp_editplots; ylabel('Change in numb. of switches');
title('Atomoxetine: Counting')

subplot(2,2,2);

m1 = nanmean(par_bttn1);
s1 = nanstd(par_bttn1)./length(par_bttn1);
m2 = nanmean(par_bttn2);
s2 = nanstd(par_bttn2)./length(par_bttn2);

bar([m1 m2]); axis([0 3 -15 15])
set(gca,'xtick',[1 2],'xticklabels',{'Block1';'Block2'}); xtickangle(90);
tp_editplots; ylabel('Change in numb. of switches');
title('Atomoxetine: Pressing')

par_cnt1 = squeeze(behav_cnt(3,:,1)-behav_cnt(1,:,1));
par_cnt2 = squeeze(behav_cnt(3,:,2)-behav_cnt(1,:,2));

par_bttn1 = squeeze(behav_bttn(3,:,1)-behav_bttn(1,:,1));
par_bttn2 = squeeze(behav_bttn(3,:,2)-behav_bttn(1,:,2));

subplot(2,2,3);

m1 = nanmean(par_cnt1);
s1 = nanstd(par_cnt1)./length(par_cnt1);
m2 = nanmean(par_cnt2);
s2 = nanstd(par_cnt2)./length(par_cnt2);

bar([m1 m2]); axis([0 3 -15 15])
set(gca,'xtick',[1 2],'xticklabels',{'Block1';'Block2'}); xtickangle(90);
tp_editplots; ylabel('Change in numb. of switches');
title('Donepezil: Counting')

subplot(2,2,4);

m1 = nanmean(par_bttn1);
s1 = nanstd(par_bttn1)./length(par_bttn1);
m2 = nanmean(par_bttn2);
s2 = nanstd(par_bttn2)./length(par_bttn2);

bar([m1 m2]); axis([0 3 -15 15])
set(gca,'xtick',[1 2],'xticklabels',{'Block1';'Block2'}); xtickangle(90);
tp_editplots; ylabel('Change in numb. of switches');
title('Donepezil: Pressing')


%% PUPIL EFFECT FOR EACH BLOCK
if ~exist('pup','var')
  para = [];
  para.time = [3000 7000];
  pup = pp_loadpupil(para);
end

figure; set(gcf,'color','w');
subplot(2,2,1);

par_cnt1  = squeeze(pup(:,2,1,2)-pup(:,1,1,2));
par_cnt2  = squeeze(pup(:,2,2,2)-pup(:,1,2,2));
par_bttn1 = squeeze(pup(:,2,1,3)-pup(:,1,1,3));
par_bttn2 = squeeze(pup(:,2,2,3)-pup(:,1,2,3));

m1 = nanmean(par_cnt1);
s1 = nanstd(par_cnt1)./length(par_cnt1);
m2 = nanmean(par_cnt2);
s2 = nanstd(par_cnt2)./length(par_cnt2);

bar([m1 m2]); axis([0 3 -500 1200])
set(gca,'xtick',[1 2],'xticklabels',{'Block1';'Block2'}); xtickangle(90);
tp_editplots; ylabel('Change in pupil diameter'); title('Atomoxetine: Counting')

subplot(2,2,2);

m1 = nanmean(par_bttn1);
s1 = nanstd(par_bttn1)./length(par_bttn1);
m2 = nanmean(par_bttn2);
s2 = nanstd(par_bttn2)./length(par_bttn2);

bar([m1 m2]); axis([0 3 -500 1200])
set(gca,'xtick',[1 2],'xticklabels',{'Block1';'Block2'}); xtickangle(90);
tp_editplots; ylabel('Change in pupil diameter'); title('Atomoxetine: Button')


par_cnt1  = squeeze(pup(:,3,1,2)-pup(:,1,1,2));
par_cnt2  = squeeze(pup(:,3,2,2)-pup(:,1,2,2));
par_bttn1 = squeeze(pup(:,3,1,3)-pup(:,1,1,3));
par_bttn2 = squeeze(pup(:,3,2,3)-pup(:,1,2,3));

subplot(2,2,3);

m1 = nanmean(par_cnt1);
s1 = nanstd(par_cnt1)./length(par_cnt1);
m2 = nanmean(par_cnt2);
s2 = nanstd(par_cnt2)./length(par_cnt2);

bar([m1 m2]); axis([0 3 -500 1200])
set(gca,'xtick',[1 2],'xticklabels',{'Block1';'Block2'}); xtickangle(90);
tp_editplots; ylabel('Change in pupil diameter'); title('Donepezil: Counting')

subplot(2,2,4);

m1 = nanmean(par_bttn1);
s1 = nanstd(par_bttn1)./length(par_bttn1);
m2 = nanmean(par_bttn2);
s2 = nanstd(par_bttn2)./length(par_bttn2);

bar([m1 m2]); axis([0 3 -500 1200])
set(gca,'xtick',[1 2],'xticklabels',{'Block1';'Block2'}); xtickangle(90);
tp_editplots; ylabel('Change in pupil diameter'); title('Donepezil: Button')

%% FC OVER BLOCKS
% --------------------------------
ifoi = 2; icond = 1; ipharm = 2; clim = [-0.02 0.02];
% --------------------------------

figure; set(gcf,'color','w','Position',[100 0 500 950]);

subplot(4,2,1); 
imagesc(nanmean(pc_blocks(:,:,:,ipharm,icond,ifoi,1),3)-nanmean(pc_blocks(:,:,:,1,icond,ifoi,1),3),clim)
set(gca,'xTick',1:1:40,'xTickLabels',reg(1:1:40),'ticklabelinterpreter','none');xtickangle(90)
set(gca,'yTick',1:1:40,'yTickLabels',reg(1:1:40),'ticklabelinterpreter','none');xtickangle(90)
axis square; tp_editplots; title('Atx: Task, Block1')

subplot(4,2,2); 
imagesc(nanmean(pc_blocks(:,:,:,ipharm,icond,ifoi,2),3)-nanmean(pc_blocks(:,:,:,1,icond,ifoi,2),3),clim)
set(gca,'xTick',1:1:40,'xTickLabels',reg(1:1:40),'ticklabelinterpreter','none');xtickangle(90)
set(gca,'yTick',1:1:40,'yTickLabels',reg(1:1:40),'ticklabelinterpreter','none');xtickangle(90)
axis square; tp_editplots; title('Atx: Task, Block2')

% --------------------------------
ipharm = 3;
% --------------------------------

subplot(4,2,3); 
imagesc(nanmean(pc_blocks(:,:,:,ipharm,icond,ifoi,1),3)-nanmean(pc_blocks(:,:,:,1,icond,ifoi,1),3),clim)
set(gca,'xTick',1:1:40,'xTickLabels',reg(1:1:40),'ticklabelinterpreter','none');xtickangle(90)
set(gca,'yTick',1:1:40,'yTickLabels',reg(1:1:40),'ticklabelinterpreter','none');xtickangle(90)
axis square; tp_editplots; title('Dpz: Rest, Block1')

subplot(4,2,4); 
imagesc(nanmean(pc_blocks(:,:,:,ipharm,icond,ifoi,2),3)-nanmean(pc_blocks(:,:,:,1,icond,ifoi,2),3),clim)
set(gca,'xTick',1:1:40,'xTickLabels',reg(1:1:40),'ticklabelinterpreter','none');xtickangle(90)
set(gca,'yTick',1:1:40,'yTickLabels',reg(1:1:40),'ticklabelinterpreter','none');xtickangle(90)
axis square; tp_editplots; title('Dpz: Rest, Block2')

colormap(cmap)

% --------------------------------
ipharm = 2;
% --------------------------------

subplot(4,2,5); 
[h,p]=ttest(pc_blocks(:,:,:,ipharm,icond,ifoi,1),pc_blocks(:,:,:,1,icond,ifoi,1),'dim',3);
imagesc((nanmean(pc_blocks(:,:,:,ipharm,icond,ifoi,1),3)-nanmean(pc_blocks(:,:,:,1,icond,ifoi,1),3)).*h,clim)
set(gca,'xTick',1:1:40,'xTickLabels',reg(1:1:40),'ticklabelinterpreter','none');xtickangle(90)
set(gca,'yTick',1:1:40,'yTickLabels',reg(1:1:40),'ticklabelinterpreter','none');xtickangle(90)
axis square; tp_editplots; title('Atx: Task, Block1')

subplot(4,2,6); 
[h,p]=ttest(pc_blocks(:,:,:,ipharm,icond,ifoi,2),pc_blocks(:,:,:,1,icond,ifoi,2),'dim',3);
imagesc((nanmean(pc_blocks(:,:,:,ipharm,icond,ifoi,2),3)-nanmean(pc_blocks(:,:,:,1,icond,ifoi,2),3)).*h,clim)
set(gca,'xTick',1:1:40,'xTickLabels',reg(1:1:40),'ticklabelinterpreter','none');xtickangle(90)
set(gca,'yTick',1:1:40,'yTickLabels',reg(1:1:40),'ticklabelinterpreter','none');xtickangle(90)
axis square; tp_editplots; title('Atx: Task, Block2')

% --------------------------------
ipharm = 3;
% --------------------------------

subplot(4,2,7); 
[h,p]=ttest(pc_blocks(:,:,:,ipharm,icond,ifoi,1),pc_blocks(:,:,:,1,icond,ifoi,1),'dim',3);
imagesc((nanmean(pc_blocks(:,:,:,ipharm,icond,ifoi,1),3)-nanmean(pc_blocks(:,:,:,1,icond,ifoi,1),3)).*h,clim)
set(gca,'xTick',1:1:40,'xTickLabels',reg(1:1:40),'ticklabelinterpreter','none');xtickangle(90)
set(gca,'yTick',1:1:40,'yTickLabels',reg(1:1:40),'ticklabelinterpreter','none');xtickangle(90)
axis square; tp_editplots; title('Dpz: Rest, Block1')

subplot(4,2,8); 
[h,p]=ttest(pc_blocks(:,:,:,ipharm,icond,ifoi,2),pc_blocks(:,:,:,1,icond,ifoi,2),'dim',3);
imagesc((nanmean(pc_blocks(:,:,:,ipharm,icond,ifoi,2),3)-nanmean(pc_blocks(:,:,:,1,icond,ifoi,2),3)).*h,clim)
set(gca,'xTick',1:1:40,'xTickLabels',reg(1:1:40),'ticklabelinterpreter','none');xtickangle(90)
set(gca,'yTick',1:1:40,'yTickLabels',reg(1:1:40),'ticklabelinterpreter','none');xtickangle(90)
axis square; tp_editplots; title('Dpz: Rest, Block2')

colormap(cmap)
print(gcf,'-dpdf',sprintf('~/pp/plots/pp_src_powcorr_fcmat_blocks_f%d_cond%d.pdf',ifoi,icond))


%%
