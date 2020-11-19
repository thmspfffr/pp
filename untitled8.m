%% pp_task_src_powcorr

% Call pupmod_all_powcorr_periphereal.m next, in order to clean
% estimated FCs from artifacts.

clear
% --------------------------------------------------------
% VERSION 12 - VOXEL LEVEL, 400 samples cortex
% --------------------------------------------------------
v               = 12;
v_postproc      = 6;
fsample         = 400;
SUBJLIST        = [4 5 6 7 8 9 10 11 12 13 15 16 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34];
allpara.filt    = 'jh_lcmv';
allpara.grid    = 'cortex_lowres';
foi_range       = 9:2:13;
para.epleng     = 60;
lpc             = 0;
timevariant     = 0;
para.wavelet    = 'bp_filt';
para.scnd_filt  = 0;
allpara.reg     = 0.05;
allpara.weigh   = 0;
allpara.tau     = nan;
% --------------------------------------------------------

if strcmp(allpara.grid,'xcoarse')
  v_grid = 2;
elseif strcmp(allpara.grid,'cortex')
  v_grid = 3;
elseif strcmp(allpara.grid,'aal')
  v_grid = 4;
elseif strcmp(allpara.grid,'medium')
  v_grid = 5;
elseif strcmp(allpara.grid,'aal_6mm')
  v_grid = 6;
elseif strcmp(allpara.grid,'aal_4mm')
  v_grid = 7;
elseif strcmp(allpara.grid,'m758_4mm')
  v_grid = 8;
elseif strcmp(allpara.grid, 'cortex_lowres')
  v_grid = 9;
elseif strcmp(allpara.grid,'genemaps')
  v_grid = 13;
elseif strcmp(allpara.grid,'genemaps_aal')
  v_grid = 14;
elseif strcmp(allpara.grid,'cortex800')
  v_grid = 16;
end

addpath /home/gnolte/meg_toolbox/toolbox/
addpath /home/gnolte/meg_toolbox/fieldtrip_utilities/
addpath /home/gnolte/meg_toolbox/toolbox_nightly/
addpath /home/gnolte/meg_toolbox/meg/
addpath /home/tpfeffer/Documents/MATLAB/fieldtrip-20130925/

outdir   = '/home/tpfeffer/pp/proc/conn/';
addpath /home/tpfeffer/pconn/matlab/

ord   = pconn_randomization;
%% LOAD DATA COMPUTE SRC TIME COURSES

for isubj = SUBJLIST
  for m = 1 : 3
    
    if ~exist(sprintf([outdir 'pp_src_powcorr_test_s%d_m%d_v%d_processing.txt'],isubj,m,v))
      system(['touch ' outdir sprintf('pp_src_powcorr_test_s%d_m%d_v%d_processing.txt',isubj,m,v)]);
    else
      continue
    end
    
    powcorr = zeros(400,400,2,length(foi_range));

    for iblock = 1:2
      
      pars = [];
      pars.sa   = sprintf('~/pconn/proc/src/pconn_sa_s%d_m%d_b%d_v%d.mat',isubj,m,iblock,v_grid);
      sa        = load(pars.sa);
      
      fprintf('Loading MEG data ...\n');
      
      try
        load(sprintf('~/pp/proc/pp_sens_cleandat_s%d_m%d_b%d_v%d.mat',isubj,m,iblock,1))
      % load(sprintf('/home/tpfeffer/pconn_cnt/proc/preproc/pconn_cnt_postpostproc_s%d_m%d_b%d_v%d.mat',isubj,m,iblock,1))
      catch me
        if ~exist(sprintf('~/pp/proc/pp_sens_cleandat_s%d_m%d_b%d_v%d.mat',isubj,m,iblock,1))
          if strcmp(allpara.grid,'cortex_lowres')
            powcorr(:,:,iblock,:) = nan(400,400,1,length(foi_range));
          end
%          	save(sprintf([outdir 'pupmod_task_src_powcorr_test_s%d_m%d_b%d_f%d_v%d.mat'],isubj,m,iblock,v),'powcorr');
          continue
        else
          error('Data corrupt?')
        end
      end
%       data.trial{1} = data.trial;
%       dat = megdata2mydata(data); 
      data.trial = data.trial(:,data.start_of_recording:data.end_of_recording-1200);
      data.trial(:,isnan(data.trial(1,:)))=[];
      dat = data.trial';
      clear data

      for ifoi = 1:length(foi_range)
            
        fprintf('Processing s%d m%d b%d f%d  ...\n', isubj,m,iblock,ifoi)

        epleng    = size(dat,1);
        epshift   = 1;
        fsample   = 400;
        segleng   = 200;
        segshift  = 200;
        f         = foi_range(ifoi);
 
        cs  = data2cs_event(dat,segleng,segshift,epleng,segleng);
        ff  = 0:fsample/segleng:fsample/2;
        idx = ff>=(f-1) & ff<= (f+1);
        
        cs = mean(real(cs(:,:,idx)),3);
        clear cs_data3 cs_data2 cs_data1
        
        para.iscs = 1;
        para.reg  = 0.05;
        filt      = pconn_beamformer(cs,sa.sa.L_coarse,para);
        
        tmp = tp_data2orthopowcorr(dat,segleng,segshift,epleng,f,fsample,filt);
        
        powcorr(:,:,iblock,ifoi) = (tmp+tmp')./2;

      end
    end
    
    save(sprintf([outdir 'pp_src_powcorr_test_s%d_m%d_v%d.mat'],isubj,m,v),'powcorr');

  end
end

error('!')

%%
v=12
clear allpowcorr
for isubj = SUBJLIST
  isubj
  for m = 1 : 3
    im = find(ord(isubj,:)==m);
%     load(sprintf([outdir 'pupmod_src_powcorr_test_s%d_m%d_v%d.mat'],isubj,im,v));
%       
%     allpowcorr(:,:,isubj,m,:,1) = single(squeeze(nanmean(powcorr,3))); clear powcorr
%     
try
    load(sprintf([outdir 'pp_task_src_powcorr_test_s%d_m%d_v%d.mat'],isubj,im,v));
      
    allpowcorr(:,:,isubj,m,:,2) = single(squeeze(nanmean(powcorr,3))); clear powcorr
catch me
  allpowcorr(:,:,isubj,m,:,2) = nan(400,400,3); clear powcorr
end
  end
end

allpowcorr = squeeze(allpowcorr(:,:,SUBJLIST,:,:,:));

%% PLOT

if ~exist('sa_meg_template','var')
  load /home/gnolte/meth/templates/mri.mat;
  load /home/gnolte/meth/templates/sa_template.mat;
  load /home/tpfeffer/pconn/proc/src/pconn_sa_s4_m1_b1_v11.mat
  grid = sa.grid_cortex_lowres;
  addpath /home/gnolte/meg_toolbox/toolbox/
  addpath /home/gnolte/meg_toolbox/fieldtrip_utilities/
  addpath /home/gnolte/meg_toolbox/toolbox_nightly/
  addpath /home/gnolte/meg_toolbox/meg/
end

close all

ifoi = 6:10; icond = 3; icont = 1;

% cmap =  hot
cmap = cbrewer('div','RdBu',128); cmap = cmap(end:-1:1,:);

par = nanmean(nanmean(nanmean(allpowcorr(:,:,:,icond,ifoi,icont),3),2),5)-nanmean(nanmean(nanmean(allpowcorr(:,:,:,1,ifoi,icont),3),2),5)%-nanmean(nanmean(p1(:,:,:,1,6),3),2);
[h] = ttest(squeeze(nanmean(nanmean(allpowcorr(:,:,:,icond,ifoi,icont),2),5)),squeeze(nanmean(nanmean(allpowcorr(:,:,:,1,ifoi,icont),2),5)),'dim',2,'alpha',0.025);%-nanmean(nanmean(p1(:,:,:,1,6),3),2)% par = par.*h;
% par(outp_atx.pval_p_atx(:,icond,ifoi)>=0.025) = 0;
par=  par.*h;
% cmap = [cmap; 0.98*ones(1,3); cmap];
para = [];
para.clim = [-max([abs(min(par)) abs(max(par))]) max([abs(min(par)) abs(max(par))])];
para.cmap = cmap;
para.grid = grid;
para.dd = 0.75;
para.fn = sprintf('~/pupmod/plots/test.png');
tp_plot_surface(par,para)

%%
for ifoi = 6:10
  
[h,~,~,s]=ttest(allpowcorr(:,:,:,3,ifoi),allpowcorr(:,:,:,1,ifoi),'dim',3)

n_pos(ifoi) = 100*sum(sum((h>0)&(s.tstat>0)))/(400*400-400);
n_neg(ifoi) = 100*sum(sum((h>0)&(s.tstat<0)))/(400*400-400);

end
%% PERMUTATION TEST

  
%   dat_res1 = single(squeeze(allpowcorr(:,:,:,[1 2],:)));
%   dat_res2 = single(squeeze(allpowcorr(:,:,:,[1 3],:)));
all_idx1 = randi(2,[size(SUBJLIST,2),1000]);

for iperm = 1 : 1000
  iperm
  % within subjects permutation test
  fprintf('Perm #%d\n',iperm);
  
  idx1 = all_idx1(:,iperm);
  idx2 = 3-idx1;
  
  %       for ifoi = 1 : 10
  %         ifoi
  for i = 1 : length(idx1)
    permdat_cnt1(:,:,i,1,:) = squeeze(allpowcorr(:,:,i,idx1(i),1:10));
    permdat_cnt1(:,:,i,2,:) = squeeze(allpowcorr(:,:,i,idx2(i),1:10));
  end
  
  [h,~,~,s]=ttest(permdat_cnt1(:,:,:,2,:),permdat_cnt1(:,:,:,1,:),'dim',3);
  
  n_pos_perm(iperm,:) = 100*sum(sum((h>0)&(s.tstat>0)))/(400*400-400);
  n_neg_perm(iperm,:) = 100*sum(sum((h>0)&(s.tstat<0)))/(400*400-400);
  
end



  












