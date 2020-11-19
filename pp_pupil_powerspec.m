%% PLOT RESULTS
SUBJLIST  = [4 5 6 7 8 9 10 11 12 13 15 16 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34];
addpath ~/pconn/matlab

%% LOAD RESTING STATE DATA

if ~exist('~/pp/proc/pup/pp_pupil_summary.mat')
  ord = pconn_randomization;
  for isubj = SUBJLIST
    isubj
    for m = 1 :3
      im = find(ord(isubj,:)==m);
      for iblock = 1 : 2
        try
          load(sprintf('~/pp/proc/pup/pp_pupil_diameter_cleaned_s%d_m%d_b%d.mat',isubj,im,iblock))
          %         allpup_cnt(isubj,m,iblock) = mean(pupil(1:10000,4));
          tmp=abs(fft(pupil(:,4))).^2;
          pup_rest.pow(:,isubj,m,iblock)=tmp(1:1000);
          pup_rest.var(isubj,m,iblock) = std(pupil(:,4))/mean(pupil(:,4));
          tmp = tp_dfa(pupil(:,4),[1 100],1000,0.5,20);
          pup_rest.dfa (isubj,m,iblock) = tmp.exp;
          pup_rest.mean_pup(isubj,m,iblock) = mean(pupil(1000:11000,4));
          
          clear pupil
        catch me
          pup_rest.pow(1:1000,isubj,m,iblock) = nan(1000,1);
          pup_rest.var(isubj,m,iblock) = nan;
          pup_rest.dfa (isubj,m,iblock) = nan;
          pup_rest.mean_pup(isubj,m,iblock) = nan;
        end
      end
    end
  end
  
  pup_rest.pow = nanmean(pup_rest.pow(:,SUBJLIST,:,:),4);
  pup_rest.var = nanmean(pup_rest.var(SUBJLIST,:,:),3);
  pup_rest.dfa = nanmean(pup_rest.dfa(SUBJLIST,:,:),3);
  pup_rest.mean_pup = nanmean(pup_rest.mean_pup(SUBJLIST,:,:),3);
  
  
  % LOAD COUNTING DATA
  
  ord = pconn_randomization;
  for isubj = SUBJLIST
    isubj
    for m = 1 :3
      im = find(ord(isubj,:)==m);
      for iblock = 1 : 2
        try
          load(sprintf('~/pp/proc/pup/pp_pupil_diameter_cleaned_counting_s%d_m%d_b%d.mat',isubj,im,iblock))
          tmp=abs(fft(pupil(:,4))).^2;
          pup_cnt.pow(:,isubj,m,iblock)=tmp(1:1000);
          pup_cnt.var(isubj,m,iblock) = std(pupil(:,4))/mean(pupil(:,4));
          tmp = tp_dfa(pupil(:,4),[1 100],1000,0.5,20);
          pup_cnt.dfa (isubj,m,iblock) = tmp.exp;
          pup_cnt.mean_pup(isubj,m,iblock) = mean(pupil(1000:11000,4));
          
          clear pupil
        catch me
          pup_cnt.pow(1:1000,isubj,m,iblock) = nan(1000,1);
          pup_cnt.var(isubj,m,iblock) = nan;
          pup_cnt.dfa (isubj,m,iblock) = nan;
          pup_cnt.mean_pup(isubj,m,iblock) = nan;
          
        end
      end
    end
  end
  
  pup_cnt.pow = nanmean(pup_cnt.pow(:,SUBJLIST,:,:),4);
  pup_cnt.var = nanmean(pup_cnt.var(SUBJLIST,:,:),3);
  pup_cnt.dfa = nanmean(pup_cnt.dfa(SUBJLIST,:,:),3);
  pup_cnt.mean_pup = nanmean(pup_cnt.mean_pup(SUBJLIST,:,:),3);
  
  % LOAD PRESSING DATA
  
  for isubj = SUBJLIST
    isubj
    for m = 1 :3
      im = find(ord(isubj,:)==m);
      for iblock = 1 : 2
        try
          load(sprintf('~/pp/proc/pup/pp_pupil_diameter_cleaned_button_s%d_m%d_b%d.mat',isubj,im,iblock))
          tmp=abs(fft(pupil(:,4))).^2;
          pup_bttn.pow(:,isubj,m,iblock)=tmp(1:1000);
          pup_bttn.var(isubj,m,iblock) = std(pupil(:,4))/mean(pupil(:,4));
          tmp = tp_dfa(pupil(:,4),[1 100],1000,0.5,20);
          pup_bttn.dfa (isubj,m,iblock) = tmp.exp;
          pup_bttn.mean_pup(isubj,m,iblock) = mean(pupil(1000:11000,4));
          clear pupil
        catch me
          pup_bttn.pow(1:1000,isubj,m,iblock) = nan(1000,1);
          pup_bttn.var(isubj,m,iblock) = nan;
          pup_bttn.dfa (isubj,m,iblock) = nan;
          pup_bttn.mean_pup(isubj,m,iblock) = nan;
        end
      end
    end
  end
  
  pup_bttn.pow = nanmean(pup_bttn.pow(:,SUBJLIST,:,:),4);
  pup_bttn.var = nanmean(pup_bttn.var(SUBJLIST,:,:),3);
  pup_bttn.dfa = nanmean(pup_bttn.dfa(SUBJLIST,:,:),3);
  pup_bttn.mean_pup = nanmean(pup_bttn.mean_pup(SUBJLIST,:,:),3);
  
  save('~/pp/proc/pup/pp_pupil_summary.mat','pup_rest','pup_cnt','pup_bttn')
  
else
  load('~/pp/proc/pup/pp_pupil_summary.mat')
end
%% PLOT EVERYTHING

cols = [0.7 0.7 0.7; 1 0.2 0; 0 0.2 1];

figure; set(gcf,'color','w');
subplot(2,3,1); hold on

m = nanmean(pup_rest.mean_pup);
s = nanstd(pup_rest.mean_pup)./sqrt(28);

[h_atx,p_atx] = ttest(pup_rest.mean_pup(:,1),pup_rest.mean_pup(:,2));
[h_dpz,p_dpz] = ttest(pup_rest.mean_pup(:,1),pup_rest.mean_pup(:,3));
[h_drg,p_drg] = ttest(pup_rest.mean_pup(:,2),pup_rest.mean_pup(:,3));

for i = 1 : 3
  bar(i,m(i),'facecolor',cols(i,:),'edgecolor',cols(i,:));
  line([i i],[m(i)-s(i) m(i)+s(i)],'color',cols(i,:),'linewidth',1)
end

axis([0 4 6500 9000]); tp_editplots; ylabel('Pupil diameter')
title('Rest'); set(gca,'xtick',[1 2 3],'xticklabels',{'Pbo';'Atx';'Dpz'})

% COUNTING
% -------------------
subplot(2,3,2); hold on

m = nanmean(pup_cnt.mean_pup);
s = nanstd(pup_cnt.mean_pup)./sqrt(28);

[h_atx,p_atx] = ttest(pup_cnt.mean_pup(:,1),pup_cnt.mean_pup(:,2));
[h_dpz,p_dpz] = ttest(pup_cnt.mean_pup(:,1),pup_cnt.mean_pup(:,3));
[h_drg,p_drg] = ttest(pup_cnt.mean_pup(:,2),pup_cnt.mean_pup(:,3));

for i = 1 : 3
  bar(i,m(i),'facecolor',cols(i,:),'edgecolor',cols(i,:));
  line([i i],[m(i)-s(i) m(i)+s(i)],'color',cols(i,:),'linewidth',1)
end

axis([0 4 6500 9000]); tp_editplots; ylabel('Pupil diameter')
title('Task-Counting'); set(gca,'xtick',[1 2 3],'xticklabels',{'Pbo';'Atx';'Dpz'})

% PRESSING
% -------------------
subplot(2,3,3); hold on

m = nanmean(pup_bttn.mean_pup);
s = nanstd(pup_bttn.mean_pup)./sqrt(28);

[h_atx,p_atx] = ttest(pup_bttn.mean_pup(:,1),pup_bttn.mean_pup(:,2));
[h_dpz,p_dpz] = ttest(pup_bttn.mean_pup(:,1),pup_bttn.mean_pup(:,3));
[h_drg,p_drg] = ttest(pup_bttn.mean_pup(:,2),pup_bttn.mean_pup(:,3));

for i = 1 : 3
  bar(i,m(i),'facecolor',cols(i,:),'edgecolor',cols(i,:));
  line([i i],[m(i)-s(i) m(i)+s(i)],'color',cols(i,:),'linewidth',1)
end

axis([0 4 6500 9000]); tp_editplots; ylabel('Pupil diameter')
title('Task-Button'); set(gca,'xtick',[1 2 3],'xticklabels',{'Pbo';'Atx';'Dpz'})

% correlation 



