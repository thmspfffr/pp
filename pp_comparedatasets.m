ord = pconn_randomization;

old_fc = zeros(400,400,34,3,2,13,2,'single');
new_fc = zeros(400,400,34,3,2,13,2,'single');

clear r_fc
for isubj = SUBJLIST
  isubj
  for iblock = 1 : 2
    for m =  1 : 3
      im = find(ord(isubj,:)==m);
      for ifoi = 1 : 13
        
        try
          load(sprintf('~/pupmod/proc/conn/pupmod_src_powcorr_s%d_m%d_b%d_f%d_v12.mat',isubj,im,iblock,ifoi))
          old_fc(:,:,isubj,m,1,ifoi,iblock) = single(powcorr); clear powcorr
        catch me
          old_fc(:,:,isubj,m,1,ifoi,iblock) = single(nan(400,400));
        end
        
        try
          load(sprintf('~/pp/proc/conn/pp_src_powcorr_s%d_m%d_b%d_f%d_v12.mat',isubj,im,iblock,ifoi))
          new_fc(:,:,isubj,m,1,ifoi,iblock) = powcorr; clear powcorr
        catch me
          new_fc(:,:,isubj,m,1,ifoi,iblock) = nan(400,400);
        end
        
        try
          load(sprintf('~/pupmod/proc/conn/pupmod_task_src_powcorr_s%d_m%d_b%d_f%d_v12.mat',isubj,im,iblock,ifoi))
          old_fc(:,:,isubj,m,2,ifoi,iblock) = single(powcorr); clear powcorr
        catch me
          old_fc(:,:,isubj,m,2,ifoi,iblock) = single(nan(400,400));
        end
        
        try
          load(sprintf('~/pp/proc/conn/pp_task_src_powcorr_s%d_m%d_b%d_f%d_v12.mat',isubj,im,iblock,ifoi))
          new_fc(:,:,isubj,m,2,ifoi,iblock) = powcorr; clear powcorr
        catch me
          new_fc(:,:,isubj,m,2,ifoi,iblock) = nan(400,400);
        end
        
      end
    end
  end
end

%% ATOMOXETINE BLOCKS

% allpowcorr = fc; clear fc

figure; set(gcf,'color','w')
SUBJ = 1:28;
subplot(2,2,1);
iblock = 1; icond = 1; ipharm = 2; ifoi = 6;
[h,~,~,s]=ttest(allpowcorr(:,:,SUBJ,ipharm,icond,ifoi,iblock),allpowcorr(:,:,SUBJ,1,icond,ifoi,iblock),'dim',3);
imagesc(h.*sign(s.tstat),[-1 1])
% imagesc(nanmean(allpowcorr(:,:,SUBJ,ipharm,icond,ifoi,iblock),3)-nanmean(allpowcorr(:,:,SUBJ,1,icond,ifoi,iblock),3),[-0.02 0.02])
axis square
subplot(2,2,2);
iblock = 2;
[h,~,~,s]=ttest(allpowcorr(:,:,SUBJ,ipharm,icond,ifoi,iblock),allpowcorr(:,:,SUBJ,1,icond,ifoi,iblock),'dim',3);
imagesc(h.*sign(s.tstat),[-1 1])
% imagesc(nanmean(allpowcorr(:,:,SUBJ,ipharm,icond,ifoi,iblock),3)-nanmean(allpowcorr(:,:,SUBJ,1,icond,ifoi,iblock),3),[-0.02 0.02])
axis square

subplot(2,2,3);
iblock = 1; icond = 2;
[h,~,~,s]=ttest(allpowcorr(:,:,SUBJ,ipharm,icond,ifoi,iblock),allpowcorr(:,:,SUBJ,1,icond,ifoi,iblock),'dim',3);
imagesc(h.*sign(s.tstat),[-1 1])
% imagesc(nanmean(allpowcorr(:,:,SUBJ,ipharm,icond,ifoi,iblock),3)-nanmean(allpowcorr(:,:,SUBJ,1,icond,ifoi,iblock),3),[-0.02 0.02])
axis square
subplot(2,2,4);
iblock = 2;
[h,~,~,s]=ttest(allpowcorr(:,:,SUBJ,ipharm,icond,ifoi,iblock),allpowcorr(:,:,SUBJ,1,icond,ifoi,iblock),'dim',3);
imagesc(h.*sign(s.tstat),[-1 1])
% imagesc(nanmean(allpowcorr(:,:,SUBJ,ipharm,icond,ifoi,iblock),3)-nanmean(allpowcorr(:,:,SUBJ,1,icond,ifoi,iblock),3),[-0.02 0.02])
axis square


%% DONEPEZIL BLOCKS
figure; set(gcf,'color','w')

SUBJ = 1:28;
subplot(2,2,1);
iblock = 1; icond = 1; ipharm = 3; ifoi = 7;
[h,~,~,s]=ttest(allpowcorr(:,:,SUBJ,ipharm,icond,ifoi,iblock),allpowcorr(:,:,SUBJ,1,icond,ifoi,iblock),'dim',3);
imagesc(h.*sign(s.tstat),[-1 1])
% imagesc(nanmean(allpowcorr(:,:,SUBJ,ipharm,icond,ifoi,iblock),3)-nanmean(allpowcorr(:,:,SUBJ,1,icond,ifoi,iblock),3),[-0.02 0.02])
axis square
subplot(2,2,2);
iblock = 2;
[h,~,~,s]=ttest(allpowcorr(:,:,SUBJ,ipharm,icond,ifoi,iblock),allpowcorr(:,:,SUBJ,1,icond,ifoi,iblock),'dim',3);
imagesc(h.*sign(s.tstat),[-1 1])
% imagesc(nanmean(allpowcorr(:,:,SUBJ,ipharm,icond,ifoi,iblock),3)-nanmean(allpowcorr(:,:,SUBJ,1,icond,ifoi,iblock),3),[-0.02 0.02])
axis square

subplot(2,2,3);
iblock = 1; icond = 2;
[h,~,~,s]=ttest(allpowcorr(:,:,SUBJ,ipharm,icond,ifoi,iblock),allpowcorr(:,:,SUBJ,1,icond,ifoi,iblock),'dim',3);
imagesc(h.*sign(s.tstat),[-1 1])
% imagesc(nanmean(allpowcorr(:,:,SUBJ,ipharm,icond,ifoi,iblock),3)-nanmean(allpowcorr(:,:,SUBJ,1,icond,ifoi,iblock),3),[-0.02 0.02])
axis square
subplot(2,2,4);
iblock = 2;
[h,~,~,s]=ttest(allpowcorr(:,:,SUBJ,ipharm,icond,ifoi,iblock),allpowcorr(:,:,SUBJ,1,icond,ifoi,iblock),'dim',3);
imagesc(h.*sign(s.tstat),[-1 1])
% imagesc(nanmean(allpowcorr(:,:,SUBJ,ipharm,icond,ifoi,iblock),3)-nanmean(allpowcorr(:,:,SUBJ,1,icond,ifoi,iblock),3),[-0.02 0.02])
axis square

subplot(2,2,4);
iblock = 2;
[h,~,~,s]=ttest(nanmean(allpowcorr(:,:,SUBJ,ipharm,icond,ifoi,:),7),nanmean(allpowcorr(:,:,SUBJ,1,icond,ifoi,:),7),'dim',3);
imagesc(h.*sign(s.tstat),[-1 1])
% imagesc(nanmean(allpowcorr(:,:,SUBJ,ipharm,icond,ifoi,iblock),3)-nanmean(allpowcorr(:,:,SUBJ,1,icond,ifoi,iblock),3),[-0.02 0.02])
axis square

%% FRACTION OF ALTERED CORRELATIONS

mask = logical(tril(ones(400,400),-1));
for ipharm = 2 : 3
  for ifoi = 1 : 13
    ifoi
    for icond = 1 : 2
      
      [h,~,~,s]=ttest(nanmean(old_fc(:,:,:,ipharm,icond,ifoi,:),7),nanmean(old_fc(:,:,:,1,icond,ifoi,:),7),'dim',3);
      
      avg_neg_old(ifoi,icond,ipharm-1) = sum(sum((h(mask)>0 & s.tstat(mask)<0)))/sum(mask(:));
      avg_pos_old(ifoi,icond,ipharm-1) = sum(sum((h(mask)>0 & s.tstat(mask)>0)))/sum(mask(:));
      
      [h,~,~,s]=ttest(nanmean(new_fc(:,:,:,ipharm,icond,ifoi,:),7),nanmean(new_fc(:,:,:,1,icond,ifoi,:),7),'dim',3);
      
      avg_neg_new(ifoi,icond,ipharm-1) = sum(sum((h(mask)>0 & s.tstat(mask)<0)))/sum(mask(:));
      avg_pos_new(ifoi,icond,ipharm-1) = sum(sum((h(mask)>0 & s.tstat(mask)>0)))/sum(mask(:));
      
    end
  end
end

%%
for i = 1 : 2
  
  if i == 1
    avg_pos = avg_pos_old;
    avg_neg = avg_neg_old;
  else
    avg_pos = avg_pos_new;
    avg_neg = avg_neg_new;
  end
  
  
  figure; set(gcf,'color','w');
  
  foi_range       = unique(round(2.^[1:.5:7]));
  
  ipharm = 1;
  
  subplot(3,2,1);
  plot(avg_pos(:,1,ipharm),'r'); hold on
  plot(avg_neg(:,1,ipharm),'b'); hold on
  axis([0 14 -0.05 0.5]); tp_editplots; %text(1,0.45,'Block #1','FontSize',7)
  ylabel(sprintf('Fraction of singificanly\nincreasedcorrelations'))
  xlabel('Carrier frequency'); box off
  set(gca,'xTick',1:2:13,'xTickLabels',num2cell([foi_range(1:2:13)]))
  
  subplot(3,2,3);
  plot(avg_pos(:,2,ipharm),'r'); hold on
  plot(avg_neg(:,2,ipharm),'b'); hold on
  axis([0 14 -0.05 0.5]); tp_editplots; %text(1,0.45,'Block #1','FontSize',7)
  ylabel(sprintf('Fraction of singificanly\nincreasedcorrelations'))
  xlabel('Carrier frequency'); box off
  set(gca,'xTick',1:2:13,'xTickLabels',num2cell([foi_range(1:2:13)]))
  
  subplot(3,2,5);
  plot(avg_pos(:,1,ipharm)-avg_pos(:,2,ipharm),'r'); hold on
  plot(avg_neg(:,1,ipharm)-avg_neg(:,2,ipharm),'b'); hold on
  axis([0 14 -0.5 0.5]); tp_editplots; %text(1,0.45,'Block #1','FontSize',7)
  ylabel(sprintf('Fraction of singificanly\nincreasedcorrelations'))
  xlabel('Carrier frequency'); box off
  set(gca,'xTick',1:2:13,'xTickLabels',num2cell([foi_range(1:2:13)]))
  
  ipharm = 2;
  
  subplot(3,2,2);
  plot(avg_pos(:,1,ipharm),'r'); hold on
  plot(avg_neg(:,1,ipharm),'b'); hold on
  axis([0 14 -0.05 0.5]); tp_editplots; %text(1,0.45,'Block #1','FontSize',7)
  ylabel(sprintf('Fraction of singificanly\nincreasedcorrelations'))
  xlabel('Carrier frequency'); box off
  set(gca,'xTick',1:2:13,'xTickLabels',num2cell([foi_range(1:2:13)]))
  
  subplot(3,2,4);
  plot(avg_pos(:,2,ipharm),'r'); hold on
  plot(avg_neg(:,2,ipharm),'b'); hold on
  axis([0 14 -0.05 0.5]); tp_editplots; %text(1,0.45,'Block #1','FontSize',7)
  ylabel(sprintf('Fraction of singificanly\nincreasedcorrelations'))
  xlabel('Carrier frequency'); box off
  set(gca,'xTick',1:2:13,'xTickLabels',num2cell([foi_range(1:2:13)]))
  
  subplot(3,2,6);
  plot(avg_pos(:,1,ipharm)-avg_pos(:,2,ipharm),'r'); hold on
  plot(avg_neg(:,1,ipharm)-avg_neg(:,2,ipharm),'b'); hold on
  axis([0 14 -0.5 0.5]); tp_editplots; %text(1,0.45,'Block #1','FontSize',7)
  ylabel(sprintf('Fraction of singificanly\nincreasedcorrelations'))
  xlabel('Carrier frequency'); box off
  set(gca,'xTick',1:2:13,'xTickLabels',num2cell([foi_range(1:2:13)]))
  
end

%%
foi_range       = unique(round(2.^[1:.5:7]));
figure; set(gcf,'color','w');

ipharm = 1;

subplot(2,3,1);
plot(pos(:,1,1,ipharm),'r:'); hold on
plot(pos(:,2,1,ipharm),'r-')
axis([0 14 -0.05 0.5]); tp_editplots; text(1,0.45,'Block #1','FontSize',7)
ylabel(sprintf('Fraction of singificanly\nincreasedcorrelations'))
xlabel('Carrier frequency'); box off
set(gca,'xTick',1:2:13,'xTickLabels',num2cell([foi_range(1:2:13)]))

subplot(2,3,2)
plot(pos(:,1,2,ipharm),'r:'); hold on
plot(pos(:,2,2,ipharm),'r-');
axis([0 14 -0.05 0.5]); tp_editplots; text(1,0.45,'Block #2','FontSize',7)
ylabel(sprintf('Fraction of singificanly\nincreasedcorrelations'))
xlabel('Carrier frequency'); box off
set(gca,'xTick',1:2:13,'xTickLabels',num2cell([foi_range(1:2:13)]))

subplot(2,3,3)
plot(pos(:,2,1,ipharm)-pos(:,1,ipharm),'color',[0.7 0.7 0.7]); hold on
plot(pos(:,2,2,ipharm)-pos(:,2,ipharm),'color',[0 0 0]); hold on
axis([0 14 -0.3 0.30001]); tp_editplots;
ylabel(sprintf('Fraction of singificanly\nincreasedcorrelations'))
xlabel('Carrier frequency'); box off
set(gca,'xTick',1:2:13,'xTickLabels',num2cell([foi_range(1:2:13)]))
line([0 14],[0 0],'color',[0.5 0.5 0.5],'linestyle','--')

subplot(2,3,4)
plot(neg(:,1,1,ipharm),'b:'); hold on
plot(neg(:,2,1,ipharm),'b-');
axis([0 14 -0.05 0.5]); tp_editplots; text(1,0.45,'Block #1','FontSize',7)
ylabel(sprintf('Fraction of singificanly\ndecreasedcorrelations'))
xlabel('Carrier frequency'); box off
set(gca,'xTick',1:2:13,'xTickLabels',num2cell([foi_range(1:2:13)]))

subplot(2,3,5)
plot(neg(:,1,2,ipharm),'b:'); hold on
plot(neg(:,2,2,ipharm),'b-');
axis([0 14 -0.05 0.5]); tp_editplots; text(1,0.45,'Block #2','FontSize',7)
ylabel(sprintf('Fraction of singificanly\ndecreasedcorrelations'))
xlabel('Carrier frequency'); box off
set(gca,'xTick',1:2:13,'xTickLabels',num2cell([foi_range(1:2:13)]))

subplot(2,3,6)
plot(neg(:,2,1,ipharm)-neg(:,1,1,ipharm),'color',[0 0 0]); hold on
plot(neg(:,2,2,ipharm)-neg(:,1,2,ipharm),'color',[0.7 0.7 0.7]); hold on
axis([0 14 -0.3  0.30001]); tp_editplots; text(1,0.2,'Block #2','FontSize',7)
ylabel(sprintf('Fraction of singificanly\ndecreasedcorrelations'))
xlabel('Carrier frequency'); box off
set(gca,'xTick',1:2:13,'xTickLabels',num2cell([foi_range(1:2:13)]))
line([0 14],[0 0],'color',[0.5 0.5 0.5],'linestyle','--')

ipharm = 2;

figure; set(gcf,'color','w');

subplot(2,3,1);
plot(pos(:,1,1,ipharm),'r:'); hold on
plot(pos(:,2,1,ipharm),'r-')
axis([0 14 -0.05 0.5]); tp_editplots; text(1,0.45,'Block #1','FontSize',7)
ylabel(sprintf('Fraction of singificanly\nincreasedcorrelations'))
xlabel('Carrier frequency'); box off
set(gca,'xTick',1:2:13,'xTickLabels',num2cell([foi_range(1:2:13)]))

subplot(2,3,2)
plot(pos(:,1,2,ipharm),'r:'); hold on
plot(pos(:,2,2,ipharm),'r-');
axis([0 14 -0.05 0.5]); tp_editplots; text(1,0.45,'Block #2','FontSize',7)
ylabel(sprintf('Fraction of singificanly\nincreasedcorrelations'))
xlabel('Carrier frequency'); box off
set(gca,'xTick',1:2:13,'xTickLabels',num2cell([foi_range(1:2:13)]))

subplot(2,3,3)
plot(pos(:,2,1,ipharm)-pos(:,1,1,ipharm),'color',[0.7 0.7 0.7]); hold on
plot(pos(:,2,2,ipharm)-pos(:,1,2,ipharm),'color',[0 0 0]); hold on
axis([0 14 -0.3 0.30001]); tp_editplots;
ylabel(sprintf('Fraction of singificanly\nincreasedcorrelations'))
xlabel('Carrier frequency'); box off
set(gca,'xTick',1:2:13,'xTickLabels',num2cell([foi_range(1:2:13)]))
line([0 14],[0 0],'color',[0.5 0.5 0.5],'linestyle','--')

subplot(2,3,4)
plot(neg(:,1,1,ipharm),'b:'); hold on
plot(neg(:,2,1,ipharm),'b-');
axis([0 14 -0.05 0.5]); tp_editplots; text(1,0.45,'Block #1','FontSize',7)
ylabel(sprintf('Fraction of singificanly\ndecreasedcorrelations'))
xlabel('Carrier frequency'); box off
set(gca,'xTick',1:2:13,'xTickLabels',num2cell([foi_range(1:2:13)]))

subplot(2,3,5)
plot(neg(:,1,2,ipharm),'b:'); hold on
plot(neg(:,2,2,ipharm),'b-');
axis([0 14 -0.05 0.5]); tp_editplots; text(1,0.45,'Block #2','FontSize',7)
ylabel(sprintf('Fraction of singificanly\ndecreasedcorrelations'))
xlabel('Carrier frequency'); box off
set(gca,'xTick',1:2:13,'xTickLabels',num2cell([foi_range(1:2:13)]))

subplot(2,3,6)
plot(neg(:,2,1,ipharm)-neg(:,1,1,ipharm),'color',[0 0 0]); hold on
plot(neg(:,2,2,ipharm)-neg(:,1,2,ipharm),'color',[0.7 0.7 0.7]); hold on
axis([0 14 -0.3  0.30001]); tp_editplots; text(1,0.2,'Block #2','FontSize',7)
ylabel(sprintf('Fraction of singificanly\ndecreasedcorrelations'))
xlabel('Carrier frequency'); box off
set(gca,'xTick',1:2:13,'xTickLabels',num2cell([foi_range(1:2:13)]))
line([0 14],[0 0],'color',[0.5 0.5 0.5],'linestyle','--')

%% COMPARE
mask = logical(tril(ones(400,400),-1));
for isubj = 1 : 28
  for iblock = 1 : 2
    for m =  1 : 3
      a= new_fc(:,:,isubj,m,2,6,iblock);
      b= old_fc(:,:,isubj,m,2,6,iblock);
      
      r(isubj,m,iblock,1) = corr(a(mask),b(mask))
      %   a= allpowcorr(:,:,isubj,m,2,6,iblock);
      %   b= fc(:,:,isubj,m,2,6,iblock);
      %   r(isubj,m,iblock,2) = corr(a(mask),b(mask))
    end
  end
end

%% CORRELATION BETWEEN BLOCKS
clear r_fc
for isubj = SUBJLIST
  for iblock = 1 : 2
    for m =  1 : 3
      
      load(sprintf('~/pupmod/proc/conn/pupmod_src_powcorr_s%d_m%d_b%d_f6_v12.mat',isubj,m,iblock))
      %       a= allpowcorr(:,:,isubj,m,1,6,iblock);
      a=powcorr; clear powcorr
      %       b= allpowcorr(:,:,isubj,m,2,7,2);
      %       r(isubj,m,1) = corr(a(mask),b(mask));
      load(sprintf('~/pp/proc/conn/pp_src_powcorr_s%d_m%d_b%d_f6_v12.mat',isubj,m,iblock))
      b = powcorr; clear powcorr
      %       b= fc(:,:,isubj,m,1,6,iblock);
      %       b= fc(:,:,isubj,m,2,7,2);
      rfc(isubj,m,iblock) = corr(a(mask),b(mask));
      
    end
  end
end
%%
iblock = 1; ipharm = 3; icond= 2; ifoi = 7


[h,~,~,s]=ttest(allpowcorr(:,:,:,ipharm,icond,ifoi,iblock),allpowcorr(:,:,:,1,icond,ifoi,iblock),'dim',3);
% sum(h(mask).*(s.tstat(mask)>0))/sum(mask(:))
sum(h(mask).*(s.tstat(mask)<0))/sum(mask(:))

% [h,~,~,s]=ttest(fc(:,:,:,ipharm,icond,ifoi,iblock),fc(:,:,:,1,icond,ifoi,iblock),'dim',3);
% sum(h(mask).*(s.tstat(mask)>0))/sum(mask(:))

[h,~,~,s]=ttest(fc(:,:,:,ipharm,icond,ifoi,iblock),fc(:,:,:,1,icond,ifoi,iblock),'dim',3);
% sum(h(mask).*(s.tstat(mask)>0))/sum(mask(:))
sum(h(mask).*(s.tstat(mask)<0))/sum(mask(:))

% [h,~,~,s]=ttest(fc(:,:,:,ipharm,icond,ifoi,iblock),fc(:,:,:,1,icond,ifoi,iblock),'dim',3);
% sum(h(mask).*(s.tstat(mask)>0))/sum(mask(:))

%
% [h,~,~,s]=ttest(nanmean(allpowcorr(:,:,:,3,1,7,:),7),nanmean(allpowcorr(:,:,:,1,1,7,:),7),'dim',3);
% sum(h(mask).*(s.tstat(mask)<0))/sum(mask(:))
%
% [h,~,~,s]=ttest(nanmean(fc(:,:,:,3,1,7,:),7),nanmean(fc(:,:,:,1,1,7,:),7),'dim',3);
% sum(h(mask).*(s.tstat(mask)<0))/sum(mask(:))

%%


%% COUNT LENGTH OF DATA
clear dif_res dif_cnt
for isubj = SUBJLIST
  isubj
  for m = 1 : 3
    im = find(ord(isubj,:)==m);
    for iblock = 1 : 2
      %     try
      load(sprintf('~/pp/proc/pp_cnt_sens_cleandat_s%d_m%d_b%d_v%d.mat',isubj,im,iblock,1))
      len_new = sum(~isnan(data.trial(1,data.start_of_recording:data.end_of_recording)));
      
      f_new = mean(abs(fft(data.trial(:,~isnan(data.trial(1,:))))).^2)';
      
      load(sprintf('/home/tpfeffer/pconn_cnt/proc/preproc/pconn_cnt_postpostproc_s%d_m%d_b%d_v%d.mat',isubj,im,iblock,1))
      len_old = size(data.trial{1}(1,:),2);
      
      %     f_old = mean(abs(fft(data.trial{1})).^2)';
      dif_res(isubj,m,iblock) = len_old - len_new;
      
      %       r_cnt_pow(isubj,m,iblock) = corr(f_new(1:100000),f_old(1:100000));
      
      %     catch me
      %       r_cnt_pow(isubj,m,iblock) = nan;
      %     end
      %     load(sprintf('~/pp/proc/pp_cnt_sens_cleandat_s%d_m%d_b%d_v%d.mat',isubj,m,iblock,1))
      %     len_new = sum(~isnan(data.trial(1,data.start_of_recording:data.end_of_recording)));
      %
      %     load(sprintf('/home/tpfeffer/pconn_cnt/proc/preproc/pconn_cnt_postpostproc_s%d_m%d_b%d_v%d.mat',isubj,m,iblock,1))
      %     len_old = size(data.trial{1}(1,:),2);
      %
      %     len_cnt(isubj,m,iblock) = size(data.trial{1}(1,:),2);
      %     dif_cnt(isubj,m,iblock) = len_old - len_new;
      
    end
  end
end


%% LOAD DATA COMPUTE SRC TIME COURSES
clc
for isubj = SUBJLIST
  fprintf('\n')
  for m =1:3
    for ifoi = 1
      
      
      for iblock = 1:2
        %         clc
        fprintf('Processing s%d m%d b%d...\n', isubj,m,iblock)
        
        %         fprintf('Loading MEG data ...\n');
        
        try
          load(sprintf('~/pp/proc/pp_cnt_sens_cleandat_s%d_m%d_b%d_v%d.mat',isubj,m,iblock,1))
        catch me
          %           if v == 20
          %             powcorr = nan(40,40);
          %           elseif v == 1
          %             powcorr = nan(91,91);
          %           elseif v == 12
          %             powcorr = nan(400,400);
          %           else
          %             error('!')
          %           end
          %           save(sprintf([outdir 'pp_src_powcorr_s%d_m%d_b%d_f%d_v%d.mat'],isubj,m,iblock,ifoi,v),'powcorr');
          continue
        end
        
        fprintf('Start: %d | End: %d | Dif: %d\n',data.start_of_recording,data.end_of_recording,data.end_of_recording-data.start_of_recording)
        
        %         dat = data.trial';
        %         if isempty(data.start_of_recording) && ~isempty(data.end_of_recording)
        %           if (data.end_of_recording-600*data.fsample)<1
        %             data.start_of_recording = 1;
        %           else
        %             data.start_of_recording = data.end_of_recording-600*data.fsample;
        %           end
        %         elseif ~isempty(data.start_of_recording) && isempty(data.end_of_recording)
        %           if (data.start_of_recording+600*data.fsample)>size(data.trial,2)
        %             data.end_of_recording = size(data.trial,2);
        %           else
        %             data.end_of_recording = data.start_of_recording+600*data.fsample;
        %           end
        %         elseif isempty(data.start_of_recording) && isempty(data.end_of_recording)
        %           data.start_of_recording = 5000;
        %           data.end_of_recording = 235000;
        %         end
        
        
      end
    end
  end
end

%%
for isubj = SUBJLIST
  isubj
  for iblock = 1 : 2
    for m = 1 :3
      try
        load(sprintf('~/pp/proc/pp_sens_cleandat_s%d_m%d_b%d_v%d.mat',isubj,m,iblock,1))
        
        if isempty(data.start_of_recording) && isempty(data.end_of_recording)
          empty(isubj,iblock,m) = 1;
        else
          empty(isubj,iblock,m) = 0;
        end
        
        if sum(isnan(data.trial(1,1:7000)))==0
          nans(isubj,iblock,m) = 0;
        else
          nans(isubj,iblock,m) = 1;
          dif(isubj,iblock,m) = find(~isnan(data.trial(1,:)),1,'first')-data.start_of_recording;
        end
      catch me
        nans_cnt(isubj,iblock,m) =nan;
        nans_cnt(isubj,iblock,m) = nan;
        dif_cnt(isubj,iblock,m) = nan;
      end

      try
        load(sprintf('~/pp/proc/pp_cnt_sens_cleandat_s%d_m%d_b%d_v%d.mat',isubj,m,iblock,1))
        
        if isempty(data.start_of_recording) && isempty(data.end_of_recording)
          empty_cnt(isubj,iblock,m) = 1;
        else
          empty_cnt(isubj,iblock,m) = 0;
        end
        
        if sum(isnan(data.trial(1,1:7000)))==0
          nans_cnt(isubj,iblock,m) = 0;
        else
          nans_cnt(isubj,iblock,m) = 1;
          dif_cnt(isubj,iblock,m) = find(~isnan(data.trial(1,:)),1,'first')-data.start_of_recording;
        end
      catch me
        nans_cnt(isubj,iblock,m) =nan;
        nans_cnt(isubj,iblock,m) = nan;
        dif_cnt(isubj,iblock,m) = nan;
        
      end
      
      
    end
  end
end
