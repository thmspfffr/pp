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
      fooof.pxx_seg_hh(:,:,:,isubj,iblock) = nan(size(fooof.pxx_seg_hh,1),246,14);
      fooof.pxx_seg_dt_hh(:,:,:,isubj,iblock) = nan(size(fooof.pxx_seg_hh,1),246,14);
      
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
    thresh = linspace(min(tmp_pup),max(tmp_pup),15);
    thresh_dt = linspace(min(tmp_pup_dt),max(tmp_pup_dt),15);
    idx = discretize(tmp_pup,thresh);
    idx_dt = discretize(tmp_pup_dt,thresh_dt);
    tmp_pxx = pxx(:,:,nanidx);
   
    for i = 1 : 14


      fooof.pxx_seg_hh(:,:,i,isubj,iblock) = mean(tmp_pxx(:,:,idx==i),3);
      fooof.pxx_seg_dt_hh(:,:,i,isubj,iblock) = mean(tmp_pxx(:,:,idx_dt==i),3);
      fooof.slopes_seg_hh(:,:,i,:,isubj,iblock) = mean(aper(2,:,idx==i),3);  
      fooof.slopes_seg_hh_dt(:,:,i,:,isubj) = mean(aper(2,:,idx_dt==i),3);
      fooof.offset_seg_hh(:,:,i,:,isubj,iblock) = mean(aper(1,:,idx==i),3);  
      fooof.offset_seg_hh_dt(:,:,i,:,isubj) = mean(aper(1,:,idx_dt==i),3);

    end
     if sum(sum(sum(fooof.pxx_seg_hh(:,:,:,isubj,iblock)==0)))>0
        a=1
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
% % fooof.ps_hh_corrected = nanmean(fooof.ps_hh_corrected,4);
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
    fooof.psfit_gla(:,iff,isubj)=corr(squeeze(g(:,iff,:))',pup(~isnan(pup))');
  end
  fooof.offset_gla(:,isubj)=corr(squeeze(aper(1,:,:))',pup(~isnan(pup))');
  fooof.slope_gla(:,isubj)=corr(squeeze(aper(2,:,:))',pup(~isnan(pup))');
  for iff = 1 : 75
    fooof.psfit_df_gla(:,iff,isubj)=corr(squeeze(g(:,iff,:))',pup_df(~isnan(pup_df))');
  end
  fooof.offset_df_gla(:,isubj)=corr(squeeze(aper(1,:,:))',pup_df(~isnan(pup_df))');
  fooof.slope_df_gla(:,isubj)=corr(squeeze(aper(2,:,:))',pup_df(~isnan(pup_df))');
  
%   take out slope of empirical pxx
  idx = ~isnan(pxx(1,1,:));
    slopes=-aper(2,:,:).*log10(2:0.5:128)' + aper(1,:,:);
    tmp = log10(pxx(:,:,idx))-slopes;
    for iff = 1 : 253
      fooof.ps_gla_corrected(:,iff,isubj)=corr(squeeze(tmp(iff,:,:))',pup(~isnan(pup))');
      fooof.ps_gla_df_corrected(:,iff,isubj)=corr(squeeze(tmp(iff,:,:))',pup_df(~isnan(pup_df))');

    end
    
    % Divide pupil signal into 5% bins
    nanidx = ~isnan(pup(:)) | ~isnan(squeeze(pxx(1,1,:)));
    tmp_pup = pup(nanidx)';
    tmp_pup_dt =  pup_df(nanidx)';
    thresh = linspace(min(tmp_pup),max(tmp_pup),15);
    thresh_dt = linspace(min(tmp_pup_dt),max(tmp_pup_dt),15);
    idx = discretize(tmp_pup,thresh);
    idx_dt = discretize(tmp_pup_dt,thresh_dt);
    tmp_pxx = pxx(:,:,nanidx);
   
    for i = 1 : 14


      pup_idx(i) = mean(tmp_pup(idx==i));
%       if isnan(mean(mean(mean(tmp_pxx(:,:,idx==i),3))))
%         a=1;
%       end
      fooof.pxx_seg_gla(:,:,i,isubj) = mean(tmp_pxx(:,:,idx==i),3);
      fooof.pxx_seg_dt_gla(:,:,i,isubj) = mean(tmp_pxx(:,:,idx_dt==i),3);
      fooof.slopes_seg_gla(:,:,i,:,isubj) = mean(aper(2,:,idx==i),3);
      fooof.slopes_seg_gla_dt(:,:,i,:,isubj) = mean(aper(2,:,idx_dt==i),3);
      fooof.offset_seg_gla(:,:,i,:,isubj,iblock) = mean(aper(1,:,idx==i),3);  
      fooof.offset_seg_gla_dt(:,:,i,:,isubj) = mean(aper(1,:,idx_dt==i),3);
    end
    if sum(sum(sum(fooof.pxx_seg_gla(:,:,:,isubj)==0)))>0
        a=1
      end
end

clear pup pup_df pxx fxx
% -------------
% LOAD MUENSTER 
% -------------
SUBJLIST=1:41; SUBJLIST([10,12,17,19,22,27,35,38,39,40])=[];
% 
% fooof.pxx_seg_mue = zeros(253,246,20,length(SUBJLIST),'double');
% fooof.pxx_seg_dt_mue = zeros(253,246,20,length(SUBJLIST),'double');
% 
for isubj = 1: length(SUBJLIST)
  load(sprintf('~/pp/proc/src/pp_mue_collected_fooof_s%d_b%d_v%d.mat',SUBJLIST(isubj),iblock,v))
  load(sprintf('~/pp/proc/src/pp_mue_src_powerspectra_s%d_b%d_v%d.mat',SUBJLIST(isubj),iblock,v))
% 
  for iff = 1 : 75
    fooof.psfit_mue(:,iff,isubj)=corr(squeeze(g(:,iff,:))',pup(~isnan(pup))');
  end
  fooof.offset_mue(:,isubj)=corr(squeeze(aper(1,:,:))',pup(~isnan(pup))');
  fooof.slope_mue(:,isubj)=corr(squeeze(aper(2,:,:))',pup(~isnan(pup))');
  for iff = 1 : 75
    fooof.psfit_df_mue(:,iff,isubj)=corr(squeeze(g(:,iff,:))',pup_df(~isnan(pup_df))');
  end
  fooof.offset_df_mue(:,isubj)=corr(squeeze(aper(1,:,:))',pup_df(~isnan(pup_df))');
  fooof.slope_df_mue(:,isubj)=corr(squeeze(aper(2,:,:))',pup_df(~isnan(pup_df))');
  
%   take out slope of empirical pxx
  idx = ~isnan(pxx(1,1,:));
    slopes=-aper(2,:,:).*log10(2:0.5:128)' + aper(1,:,:);
    tmp = log10(pxx(:,:,idx))-slopes;
    for iff = 1 : 253
      fooof.ps_mue_corrected(:,iff,isubj)=corr(squeeze(tmp(iff,:,:))',pup(~isnan(pup))');
      fooof.ps_mue_df_corrected(:,iff,isubj)=corr(squeeze(tmp(iff,:,:))',pup_df(~isnan(pup_df))');

    end
%     
%     % Divide pupil signal into 5% bins
    nanidx = ~isnan(pup(:)) | ~isnan(squeeze(pxx(1,1,:)));
    tmp_pup = pup(nanidx)';
    tmp_pup_dt =  pup_df(nanidx)';
    thresh = linspace(min(tmp_pup),max(tmp_pup),15);
    thresh_dt = linspace(min(tmp_pup_dt),max(tmp_pup_dt),15);
    idx = discretize(tmp_pup,thresh);
    idx_dt = discretize(tmp_pup_dt,thresh_dt);
    tmp_pxx = pxx(:,:,nanidx);
   
    for i = 1 : 14

      fooof.pxx_seg_mue(:,:,i,isubj) = mean(tmp_pxx(:,:,idx==i),3);
      fooof.pxx_seg_dt_mue(:,:,i,isubj) = mean(tmp_pxx(:,:,idx_dt==i),3);
      fooof.slopes_seg_mue(:,:,i,:,isubj) = mean(aper(2,:,idx==i),3);
      fooof.slopes_seg_mue_dt(:,:,i,:,isubj) = mean(aper(2,:,idx_dt==i),3);
      fooof.offset_seg_mue(:,:,i,:,isubj,iblock) = mean(aper(1,:,idx==i),3);  
      fooof.offset_seg_mue_dt(:,:,i,:,isubj) = mean(aper(1,:,idx_dt==i),3);
    end
end

