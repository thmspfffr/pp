%% pp_pupil_pharmadiff
% compute mean pupil diameter across pharma conditions
% first determine what was recorded (area or diameter) and then apply
% zscoring (across blocks and sessions) or compare raw values.

% 04/2019
% ord   = pconn_randomization;
for isubj = 4 : 34
  for m = 1 : 3
    
    % ------------------------------
    % DETERMINE TYPE OF RECORDING FIRST
    % ------------------------------
    d_evts  = dir(sprintf('/home/tpfeffer/pconn/rawdata/edf/p%d/s%d/edf_events/*rest*',isubj,m));
    
    if length(d_evts)==0
      continue
    end
    
    fn = fopen([d_evts(1).folder '/' d_evts(1).name]);
    while 1
      tmp = fgetl(fn);
      tmp = tmp(find(~isspace(tmp)));
      if ~isempty(findstr(tmp,'PUPILAREA'))
        type = 1;
        break
      elseif ~isempty(findstr(tmp,'PUPILDIAMETER'))
        type = 2;
        break
      end
      % if neither 1 nor 2, set 0
      type = 0;
    end
    fclose(fn);
    if type == 1
      fprintf('Determined type: AREA ...\n')
    else
      fprintf('Determined type: DIAMETER ...\n')
    end
    % ------------------------------
    
    for iblock = [1 2]
      clear pupil
      %       im = find(ord(isubj,:)==m);
      try
        load(sprintf('/home/arussmann/pupil_data/pconn_pup_samples_s%d_b%d_m%d_no_blinks.mat',isubj,iblock,m))
        isubj
        if type == 1
          pupil = 256.*sqrt(pupil./pi);
        end
        save(sprintf('~/pp/proc/pp_pupil_noblinks_converted_alena_s%d_b%d_m%d_no_blinks.mat',isubj,iblock,m),'pupil')
      catch me
        
      end
    end
  end
  
end

%%
addpath ~/pconn/matlab
ord = pconn_randomization;
% ~/pp/proc/pp_pupil_noblinks_converted_alena
for isubj = 4 : 34
  d = dir(sprintf('~/pp/proc/pp_pupil_noblinks_converted_alena_s%d*no_blinks.mat',isubj));
  tmp_pup = [];
  l = 0;
  for m = 1 : 3
    im = find(ord(isubj,:)==m);
    for iblock = [1 2]
      try
       
        load(sprintf('~/pp/proc/pp_pupil_noblinks_converted_alena_s%d_b%d_m%d_no_blinks.mat',isubj,iblock,im))
        ex(m,iblock)  = 1;
        len(m,iblock) = length(pupil);

        pup(isubj,m,iblock) = mean(pupil);
        tmp_pup = [tmp_pup; pupil];
        
      catch me
        pup(isubj,m,iblock) = nan;
        ex(m,iblock)  = 0;
        len(m,iblock) = 0;
      end
    end
  end
  tmp_pup = zscore(tmp_pup);
  for m = 1 : 3
    for iblock = 1 : 2
      
      
      
      if ex(m,iblock)==1
        pup_z(isubj,m,iblock) = mean(tmp_pup(1:len(m,iblock)));
        tmp_pup(1:len(m,iblock))=[];
      else
        pup_z(isubj,m,iblock) = nan;
      end
      length(tmp_pup)
    end
  end
end

pp=nanmean(pup(SUBJLIST,:,:),3);
pp_z=nanmean(pup_z(SUBJLIST,:,:),3);

%% PLOT RESULTS
col = [0.8 0.8 0.8; 1 0.2 0; 0 0.2 1]

m = nanmean(pp);
s = nanstd(pp)/sqrt(28);
mz = nanmean(pp_z);
sz = nanstd(pp_z)/sqrt(28);

figure; set(gcf,'color','w');
subplot(1,2,1); hold on
for i = 1 : 3
bar(i,m(i),'facecolor',col(i,:),'edgecolor',col(i,:)); 
line([i i],[m(i)-s(i) m(i)+s(i)],'color',col(i,:))
end
tp_editplots; axis([0 4 6000 7500]); axis square
set(gca,'xtick',[1 2 3],'xticklabels',{'Pbo';'Atx';'Dpz'});
ylabel('Mean pupil diameter'); title('Not z-scored')

subplot(1,2,2); hold on
for i = 1 : 3
bar(i,mz(i),'facecolor',col(i,:),'edgecolor',col(i,:)); 
line([i i],[mz(i)-sz(i) mz(i)+sz(i)],'color',col(i,:))
end
tp_editplots; axis([0 4 -0.5 0.5]); axis square
set(gca,'xtick',[1 2 3],'xticklabels',{'Pbo';'Atx';'Dpz'});
ylabel('Mean pupil diameter (z-scored)'); title('Z-scored')

print(gcf,'-dpdf',sprintf('~/pupil.pdf'))

%% HOW CHANGES IN PUPIL RELATE TO CHANGES IN OTHER PARAMETER
% heart rate, behavior, fc

% 
a=100*(pp(~isnan(behav(:,2)),2)-pp(~isnan(behav(:,2)),1))./pp(~isnan(behav(:,2)),1)
% a=pp_z(~isnan(behav(:,2)),2)-pp_z(~isnan(behav(:,2)),1)

b=100*(behav1(~isnan(behav(:,2)),2)-behav1(~isnan(behav(:,2)),1))./(behav1(~isnan(behav(:,2)),1))
c=100*(hb(~isnan(behav(:,2)),2)-hb(~isnan(behav(:,2)),1))./hb(~isnan(behav(:,2)),1)

% 
% % ------------
% BEHAVIOR
% ------------
SUBJLIST1  = [4 5 6 7 8 9 10 11 12 13 15 16 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34];
para.str_behav = 'count';
behav = pconn_read_behavioral_data(SUBJLIST1,para);
behav_cnt = nanmean(behav,3)';

delta_pupil     = 100*(pp(:,2)-pp(:,1))./pp(:,1);
delta_behav_cnt = 100*(behav_cnt(:,2)-behav_cnt(:,1))./(behav_cnt(:,1));
[r,p]           = corr(delta_pupil,delta_behav_cnt);

para.str_behav = 'numb_switches';
behav = pconn_read_behavioral_data(SUBJLIST1,para);
behav_bttn = nanmean(behav,3)';

delta_pupil     = 100*(pp(~isnan(behav_bttn(:,2)),2)-pp(~isnan(behav_bttn(:,2)),1))./pp(~isnan(behav_bttn(:,2)),1);
delta_behav_bttn = 100*(behav_bttn(~isnan(behav_bttn(:,2)),2)-behav_bttn(~isnan(behav_bttn(:,2)),1))./(behav_bttn(~isnan(behav_bttn(:,2)),1));
[r,p]           = corr(delta_pupil,delta_behav_bttn);

behav = (behav_cnt+behav_bttn)./2;
behav(isnan(behav_bttn))=behav_cnt(isnan(behav_bttn));
delta_pupil1    = 100*(pp(:,2)-pp(:,1))./pp(:,1);
delta_behav     = 100*(behav(:,2)-behav(:,1))./(behav(:,1));
[r,p]           = corr(delta_pupil,delta_behav);

figure; set(gcf,'color','w');
scatter(delta_pupil,delta_behav,75,'markerfacecolor','k','markeredgecolor','w')
axis([- 42 42 -122 122])
set(gca,'XGrid','off')
lsline; tp_editplots
line([-40 40],[0 0],'linestyle',':','color','k')
line([0 0],[-120 120],'linestyle',':','color','k')
xlabel('\Delta(Pupil) [in %]');ylabel('\Delta(Switch rate) [in %]')
text(-38,100,sprintf('r = %.3f | p = %.3f',r,p))

print(gcf,'-dpdf',sprintf('~/pupil_behav.pdf',i))
%%
figure; set(gcf,'color','w');
subplot(1,2,1); 
scatter(delta_pupil1,delta_behav_cnt,75,'markerfacecolor','k','markeredgecolor','w');
axis square; axis([-40 40 -200 200])
lsline; line([0 0],[-200 200],'linestyle',':','color','k')
lsline; line([-40 40],[-0 0],'linestyle',':','color','k')
[r,p]=corr(delta_pupil1,delta_behav_cnt,'type','spearman')
xlabel('Change in switches [%]'); ylabel('Change in pupil [%]')
tp_editplots
subplot(1,2,2);
scatter(delta_pupil,delta_behav_bttn,75,'markerfacecolor','k','markeredgecolor','w')
axis square; axis([-40 40 -200 200])
lsline; line([0 0],[-200 200],'linestyle',':','color','k')
lsline; line([-40 40],[-0 0],'linestyle',':','color','k')
[r,p]=corr(delta_pupil,delta_behav_bttn,'type','spearman')
xlabel('Change in switches [%]'); ylabel('Change in pupil [%]')
tp_editplots
%%
% CHANGES IN FC
% ------
% fc = pupmod_loadpowcorr(1,1);
% fc_rest = fc(:,:,:,1,1,6);
% for i = 1 : 91
%   for j = 1 : 91
%     r(i,j) = corr(squeeze(fc_rest(i,j,:)),pp(:,1));
%   end
% end
% 
% fc_task = fc(:,:,:,1,2,6);
% for i = 1 : 91
%   for j = 1 : 91
%     r_task(i,j) = corr(squeeze(fc_task(i,j,:)),pp(:,1));
%   end
% end





behav = (behav1 + behav)./2;
