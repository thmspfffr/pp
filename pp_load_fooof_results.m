function fooof  = pp_load_fooof_results(v)

load(sprintf('~/pp/proc/pp_atlas_BNA.mat'))

% LOAD HAMBURG
% -------------
SUBJLIST = [4 5 6 7 8 9 10 11 12 13 15 16 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34];
for isubj = 1: length(SUBJLIST)
  isubj
  for iblock = 1 : 2
    try
      load(sprintf('~/pp/proc/src/pp_hh_collected_fooof_s%d_b%d_v%d.mat',SUBJLIST(isubj),iblock,v))
      load(sprintf('~/pp/proc/src/pp_hh_src_powerspectra_s%d_b%d_v%d.mat',SUBJLIST(isubj),iblock,v))
      for iff = 1 : 75
        fooof.psfit_hh(:,iff,isubj,iblock)=corr(squeeze(g(:,iff,:))',pup(~isnan(pup))');
      end
      fooof.offset_hh(:,isubj,iblock)=corr(squeeze(aper(1,:,:))',pup(~isnan(pup))');
      fooof.slope_hh(:,isubj,iblock)=corr(squeeze(aper(2,:,:))',pup(~isnan(pup))');
      for iff = 1 : 75
        fooof.psfit_df_hh(:,iff,isubj,iblock)=corr(squeeze(g(:,iff,:))',pup_df(~isnan(pup_df))');
      end
      fooof.offset_df_hh(:,isubj,iblock)=corr(squeeze(aper(1,:,:))',pup_df(~isnan(pup_df))');
      fooof.slope_df_hh(:,isubj,iblock)=corr(squeeze(aper(2,:,:))',pup_df(~isnan(pup_df))');
    catch me
      fooof.psfit_hh(:,:,isubj,iblock)=nan(246,75);
      fooof.offset_hh(:,isubj,iblock)=nan(246,1);
      fooof.slope_hh(:,isubj,iblock)=nan(246,1);
      fooof.psfit_df_hh(:,:,isubj,iblock)=nan(246,75);
      fooof.offset_df_hh(:,isubj,iblock)=nan(246,1);
      fooof.slope_df_hh(:,isubj,iblock)=nan(246,1);
      continue
    end
    % take out slope of empirical pxx
    idx = ~isnan(pxx(1,1,:));
    slopes=-aper(2,:,:).*log10(2:0.5:128)' + aper(1,:,:);
    tmp = log10(pxx(:,:,idx))-slopes;
    for iff = 1 : 253
      fooof.ps_hh_corrected(:,iff,isubj,iblock)=corr(squeeze(tmp(iff,:,:))',pup(~isnan(pup))');
      fooof.ps_hh_df_corrected(:,iff,isubj,iblock)=corr(squeeze(tmp(iff,:,:))',pup_df(~isnan(pup_df))');

    end

    % Divide pupil signal into 5% bins
    nanidx = ~isnan(pup(:)) | ~isnan(squeeze(pxx(1,1,:)));
    tmp_pup = pup(nanidx)';
    tmp_pup_dt =  pup_df(nanidx)';
    for i = 1 : 21
      thresh(i) = prctile(tmp_pup,5*(i-1));
      thresh_dt(i) = prctile(tmp_pup_dt,5*(i-1));
    end
    idx = discretize(tmp_pup,thresh);
    idx_dt = discretize(tmp_pup_dt,thresh_dt);
    tmp_pxx = pxx(:,:,nanidx);
    for i = 1 : 20
      fooof.pxx_seg_hh(:,:,i,isubj,iblock) = mean(tmp_pxx(:,:,idx==i),3);
      fooof.pxx_seg_dt_hh(:,:,i,isubj,iblock) = mean(tmp_pxx(:,:,idx_dt==i),3);
      fooof.slopes_seg_hh(:,:,i,:,isubj,iblock) = mean(aper(2,:,idx==i),3);
    end
    fooof.fxx = fxx;
    
  end
end

fooof.psfit_hh=nanmean(fooof.psfit_hh,4);
fooof.offset_hh=nanmean(fooof.offset_hh,3);
fooof.slope_hh=nanmean(fooof.slope_hh,3);
fooof.psfit_df_hh=nanmean(fooof.psfit_df_hh,4);
fooof.offset_df_hh=nanmean(fooof.offset_df_hh,3);
fooof.slope_df_hh=nanmean(fooof.slope_df_hh,3);
fooof.ps_hh_corrected = nanmean(fooof.ps_hh_corrected,4);
clear pup pup_df pxx fxx

% -------------
% LOAD GLASGOW 
% -------------
SUBJLIST=1:24; SUBJLIST([5,9]) = []; iblock = 1;

for isubj = 1: length(SUBJLIST)
  load(sprintf('~/pp/proc/src/pp_gla_collected_fooof_s%d_b%d_v%d.mat',SUBJLIST(isubj),iblock,v))
  load(sprintf('~/pp/proc/src/pp_gla_src_powerspectra_s%d_b%d_v%d.mat',SUBJLIST(isubj),iblock,v))

  isubj
  for iff = 1 : 75
    fooof.psfit_gla(:,iff,isubj,iblock)=corr(squeeze(g(:,iff,:))',pup(~isnan(pup))');
  end
  fooof.offset_gla(:,isubj,iblock)=corr(squeeze(aper(1,:,:))',pup(~isnan(pup))');
  fooof.slope_gla(:,isubj,iblock)=corr(squeeze(aper(2,:,:))',pup(~isnan(pup))');
  for iff = 1 : 75
    fooof.psfit_df_gla(:,iff,isubj,iblock)=corr(squeeze(g(:,iff,:))',pup_df(~isnan(pup_df))');
  end
  fooof.offset_df_gla(:,isubj,iblock)=corr(squeeze(aper(1,:,:))',pup_df(~isnan(pup_df))');
  fooof.slope_df_gla(:,isubj,iblock)=corr(squeeze(aper(2,:,:))',pup_df(~isnan(pup_df))');
  
  % take out slope of empirical pxx
  idx = ~isnan(pxx(1,1,:));
    slopes=-aper(2,:,:).*log10(2:0.5:128)' + aper(1,:,:);
    tmp = log10(pxx(:,:,idx))-slopes;
    for iff = 1 : 253
      fooof.ps_gla_corrected(:,iff,isubj,iblock)=corr(squeeze(tmp(iff,:,:))',pup(~isnan(pup))');
      fooof.ps_gla_df_corrected(:,iff,isubj,iblock)=corr(squeeze(tmp(iff,:,:))',pup_df(~isnan(pup_df))');

    end
    
    % Divide pupil signal into 5% bins
    nanidx = ~isnan(pup(:)) | ~isnan(squeeze(pxx(1,1,:)));
    tmp_pup = pup(nanidx)';
    tmp_pup_dt =  pup_df(nanidx)';
    for i = 1 : 21
      thresh(i) = prctile(tmp_pup,5*(i-1));
      thresh_dt(i) = prctile(tmp_pup_dt,5*(i-1));
    end
    idx = discretize(tmp_pup,thresh);
    idx_dt = discretize(tmp_pup_dt,thresh_dt);
    tmp_pxx = pxx(:,:,nanidx);
    for i = 1 : 20
      fooof.pxx_seg_gla(:,i,isubj,iblock) = mean(mean(tmp_pxx(:,:,idx==i),3),2);
      fooof.pxx_seg_dt_gla(:,i,isubj,iblock) = mean(mean(tmp_pxx(:,:,idx_dt==i),3),2);
      fooof.slopes_seg_gla(:,i,isubj,iblock) = mean(mean(aper(2,:,idx==i),3),2);
    end
end

clear pup pup_df pxx fxx
% -------------
% LOAD MUENSTER 
% -------------
SUBJLIST=1:41; SUBJLIST([4,10,12,17,19,22,27,35,38,39,40])=[]; iblock = 1;
for isubj = 1: length(SUBJLIST)
  load(sprintf('~/pp/proc/src/pp_mue_collected_fooof_s%d_b%d_v%d.mat',SUBJLIST(isubj),iblock,v))
  load(sprintf('~/pp/proc/src/pp_mue_src_powerspectra_s%d_b%d_v%d.mat',SUBJLIST(isubj),iblock,v))

  for iff = 1 : 75
    fooof.psfit_mue(:,iff,isubj,iblock)=corr(squeeze(g(:,iff,:))',pup(~isnan(pup))');
  end
  fooof.offset_mue(:,isubj,iblock)=corr(squeeze(aper(1,:,:))',pup(~isnan(pup))');
  fooof.slope_mue(:,isubj,iblock)=corr(squeeze(aper(2,:,:))',pup(~isnan(pup))');
  for iff = 1 : 75
    fooof.psfit_df_mue(:,iff,isubj,iblock)=corr(squeeze(g(:,iff,:))',pup_df(~isnan(pup_df))');
  end
  fooof.offset_df_mue(:,isubj,iblock)=corr(squeeze(aper(1,:,:))',pup_df(~isnan(pup_df))');
  fooof.slope_df_mue(:,isubj,iblock)=corr(squeeze(aper(2,:,:))',pup_df(~isnan(pup_df))');
  
  % take out slope of empirical pxx
  idx = ~isnan(pxx(1,1,:));
    slopes=-aper(2,:,:).*log10(2:0.5:128)' + aper(1,:,:);
    tmp = log10(pxx(:,:,idx))-slopes;
    for iff = 1 : 253
      fooof.ps_mue_corrected(:,iff,isubj,iblock)=corr(squeeze(tmp(iff,:,:))',pup(~isnan(pup))');
      fooof.ps_mue_df_corrected(:,iff,isubj,iblock)=corr(squeeze(tmp(iff,:,:))',pup_df(~isnan(pup_df))');

    end
    
    % Divide pupil signal into 5% bins
    nanidx = ~isnan(pup(:)) | ~isnan(squeeze(pxx(1,1,:)));
    tmp_pup = pup(nanidx)';
    tmp_pup_dt =  pup_df(nanidx)';
    for i = 1 : 21
      thresh(i) = prctile(tmp_pup,5*(i-1));
      thresh_dt(i) = prctile(tmp_pup_dt,5*(i-1));
    end
    idx = discretize(tmp_pup,thresh);
    idx_dt = discretize(tmp_pup_dt,thresh_dt);
    tmp_pxx = pxx(:,:,nanidx);
    for i = 1 : 20
      fooof.pxx_seg_mue(:,i,isubj,iblock) = mean(mean(tmp_pxx(:,:,idx==i),3),2);
      fooof.pxx_seg_dt_mue(:,i,isubj,iblock) = mean(mean(tmp_pxx(:,:,idx_dt==i),3),2);
      fooof.slopes_seg_mue(:,i,isubj,iblock) = mean(mean(aper(2,:,idx==i),3),2);
    end
end


%% LOAD POWER

% pow.ff = 2:1/(800/400):128;
% 
% SUBJLIST = 1:24; SUBJLIST([5 9])=[];
% for isubj = 1:length(SUBJLIST)
%   for iblock = 1 : 1
% 
%     load(sprintf('~/pp/proc/src/pp_gla_src_powerspectra_s%d_b%d_v%d.mat',SUBJLIST(isubj),iblock,2))
% 
%     idx = ~isnan(squeeze(pxx(1,1,:)));
%     pxx = pxx(:,:,idx);
%     [ii,i]=sort(pup(idx));
% 
%     for ivox = 1: 246
%       tmp_lo(:,ivox,isubj,iblock) = nanmean(pxx(:,ivox,i(1:floor(length(i)/3)-1)),3);
%       tmp_hi(:,ivox,isubj,iblock) = nanmean(pxx(:,ivox,i((length(i)-floor(length(i)/3)):length(i))),3);
%       tmp_me(:,ivox,isubj,iblock) = nanmean(pxx(:,ivox,i(floor(length(i)/3):(length(i)-floor(length(i)/3)-1))),3);
%     end
% 
%   end
% end
% 
% tmp_lo((pow.ff>48 & pow.ff<52) | (pow.ff>98 & pow.ff<102),:,:) = nan; 
% tmp_hi((pow.ff>48 & pow.ff<52) | (pow.ff>98 & pow.ff<102),:,:) = nan;
% tmp_me((pow.ff>48 & pow.ff<52) | (pow.ff>98 & pow.ff<102),:,:) = nan;
% 
% fooof.ps_gla(:,:,:,1)  = nanmean(tmp_lo,4);
% fooof.ps_gla(:,:,:,2)  = nanmean(tmp_me,4);
% fooof.ps_gla(:,:,:,3)  = nanmean(tmp_hi,4);
% 
% % --------------------
% % LOAD HAMBURG
% % --------------------
% 
% SUBJLIST = [4 5 6 7 8 9 10 11 12 13 15 16 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34];
% 
% for isubj =1:length(SUBJLIST)
%   isubj
%   d = dir(sprintf('~/pp/proc/src/pp_hh_src_powerspectra_s%d_b*_v%d.mat',SUBJLIST(isubj),v));
%   
%   for iblock = 1:length(d)
%     load([d(iblock).folder '/' d(iblock).name])
%     idx = ~isnan(squeeze(pxx(1,1,:)));
%     pxx = pxx(:,:,idx);
%     [ii,i]=sort(pup(idx));
% 
%     for ivox = 1: 246
%       tmp_lo(:,ivox,isubj,iblock) = nanmean(pxx(:,ivox,i(1:floor(length(i)/3)-1)),3);
%       tmp_hi(:,ivox,isubj,iblock) = nanmean(pxx(:,ivox,i((length(i)-floor(length(i)/3)):length(i))),3);
%       tmp_me(:,ivox,isubj,iblock) = nanmean(pxx(:,ivox,i(floor(length(i)/3):(length(i)-floor(length(i)/3)-1))),3);
%     end
% 
%   end
% end
% 
% tmp_lo((pow.ff>48 & pow.ff<52) | (pow.ff>98 & pow.ff<102),:,:) = nan; 
% tmp_hi((pow.ff>48 & pow.ff<52) | (pow.ff>98 & pow.ff<102),:,:) = nan;
% tmp_me((pow.ff>48 & pow.ff<52) | (pow.ff>98 & pow.ff<102),:,:) = nan;
% 
% fooof.ps_hh(:,:,:,1)  = nanmean(tmp_lo,4);
% fooof.ps_hh(:,:,:,2)  = nanmean(tmp_me,4);
% fooof.ps_hh(:,:,:,3)  = nanmean(tmp_hi,4);
% 
% % --------------------
% % LOAD MUENSTER
% % --------------------
% 
% SUBJLIST = 1 : 41; SUBJLIST([10,12,17,19,22,27,35,38,39,40])=[];
% 
% for isubj =1:length(SUBJLIST)
%   for iblock = 1 : 1
%     load(sprintf('~/pp/proc/src/pp_mue_src_powerspectra_s%d_b%d_v%d.mat',SUBJLIST(isubj),iblock,2))
% 
%     idx = ~isnan(squeeze(pxx(1,1,:)));
%     pxx = pxx(:,:,idx);
%     [ii,i]=sort(pup(idx));
% 
%     for ivox = 1: 246
% 
%       tmp_lo(:,ivox,isubj,iblock) = nanmean(pxx(:,ivox,i(1:floor(length(i)/3)-1)),3);
%       tmp_hi(:,ivox,isubj,iblock) = nanmean(pxx(:,ivox,i((length(i)-floor(length(i)/3)):length(i))),3);
%       tmp_me(:,ivox,isubj,iblock) = nanmean(pxx(:,ivox,i(floor(length(i)/3):(length(i)-floor(length(i)/3)-1))),3);
%     end
% 
%   end
% end
% 
% tmp_lo((pow.ff>48 & pow.ff<52) | (pow.ff>98 & pow.ff<102),:,:) = nan; 
% tmp_hi((pow.ff>48 & pow.ff<52) | (pow.ff>98 & pow.ff<102),:,:) = nan;
% tmp_me((pow.ff>48 & pow.ff<52) | (pow.ff>98 & pow.ff<102),:,:) = nan;
% 
% fooof.ps_mue(:,:,:,1)  = nanmean(tmp_lo,4);
% fooof.ps_mue(:,:,:,2)  = nanmean(tmp_me,4);
% fooof.ps_mue(:,:,:,3)  = nanmean(tmp_hi,4);

