%% pp_task_src_powcorr

clear 

% --------------------------------------------------------
% VERSION 1 - WEIGHTED AAL
% --------------------------------------------------------
% v               = 1;
% v_postproc      = 6;
% fsample         = 400;
% SUBJLIST        = [4 5 6 7 8 9 10 11 12 13 15 16 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34];
% allpara.filt    = 'jh_lcmv';
% allpara.grid    = 'aal_4mm';
% foi_range       = unique(round(2.^[1:.5:7]));
% para.segleng    = 9 ./ foi_range;
% para.bpfreq     = [foi_range-(foi_range./2)/2; foi_range+(foi_range./2)/2]';
% para.epleng     = 60;
% lpc             = 0;
% timevariant     = 0;
% para.wavelet    = 'bp_filt';
% para.scnd_filt  = 0;
% allpara.reg     = 0.05;
% allpara.weigh   = 1;
% allpara.tau     = nan;
% --------------------------------------------------------
% VERSION 12 - VOXEL LEVEL, 400 samples cortex
% --------------------------------------------------------
v               = 12;
v_postproc      = 6;
fsample         = 400;
SUBJLIST        = [4 5 6 7 8 9 10 11 12 13 15 16 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34];
allpara.filt    = 'jh_lcmv';
allpara.grid    = 'cortex_lowres';
foi_range       = unique(round(2.^[1:.5:7]));
para.segleng    = 9 ./ foi_range;
para.bpfreq     = [foi_range-(foi_range./2)/2; foi_range+(foi_range./2)/2]';
para.epleng     = 60;
lpc             = 0;
timevariant     = 0;
para.wavelet    = 'bp_filt';
para.scnd_filt  = 0;
allpara.reg     = 0.05;
allpara.weigh   = 0;
allpara.tau     = nan;
% --------------------------------------------------------
% VERSION 20 - VTPM
% --------------------------------------------------------
% v               = 20;
% fsample         = 400;
% SUBJLIST        = [4 5 6 7 8 9 10 11 12 13 15 16 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34];
% allpara.filt    = 'jh_lcmv';
% allpara.grid    = 'vtpm_4mm';
% foi_range       = unique(round(2.^[1:.5:7]));
% para.segleng    = 9 ./ foi_range;
% para.bpfreq     = [foi_range-(foi_range./2)/2; foi_range+(foi_range./2)/2]';
% para.epleng     = 60;
% lpc             = 0;
% timevariant     = 0;
% para.wavelet    = 'bp_filt';
% para.scnd_filt  = 0;
% allpara.reg     = 0.05;
% allpara.weigh   = 1;
% allpara.tau     = nan;
% --------------------------------------------------------

addpath /home/gnolte/meg_toolbox/toolbox/
addpath /home/gnolte/meg_toolbox/fieldtrip_utilities/
addpath /home/gnolte/meg_toolbox/toolbox_nightly/
addpath /home/gnolte/meg_toolbox/meg/
% addpath /home/tpfeffer/Documents/MATLAB/fieldtrip-20130925/

outdir   = '/home/tpfeffer/pp/proc/conn/';
addpath /home/tpfeffer/pconn/matlab/
run ~/Documents/MATLAB/toolboxes/NBT-NBTv0.5.3-alpha/installNBT.m
siginfo = nbt_Info;
siginfo.converted_sample_frequency = 400;

if strcmp(allpara.grid,'xcoarse')
  v_grid = 2;
elseif strcmp(allpara.grid,'aal')
  v_grid = 4;
elseif strcmp(allpara.grid,'cortex')
  v_grid = 3;
elseif strcmp(allpara.grid,'medium')
  v_grid = 5;
elseif strcmp(allpara.grid,'aal_6mm')
  v_grid = 6;
elseif strcmp(allpara.grid,'aal_4mm')
  v_grid = 7;
elseif strcmp(allpara.grid,'m758_4mm')
  v_grid = 8;
elseif strcmp(allpara.grid,'cortex_lowres')
  v_grid = 9;
elseif strcmp(allpara.grid,'vtpm_4mm')
  v_grid = 10;
elseif strcmp(allpara.grid,'genemaps')
  v_grid = 13;
elseif strcmp(allpara.grid,'genemaps_aal')
  v_grid = 14;
elseif strcmp(allpara.grid,'cortex800')
  v_grid = 16;
end


%% LOAD DATA COMPUTE SRC TIME COURSES

for isubj = SUBJLIST
  for m = 1:3
    for ifoi = 1:13
      
      if ~exist(sprintf([outdir 'pp_task_src_powcorr_s%d_m%d_f%d_v%d_processing.txt'],isubj,m,ifoi,v))
        system(['touch ' outdir sprintf('pp_task_src_powcorr_s%d_m%d_f%d_v%d_processing.txt',isubj,m,ifoi,v)]);
      else
        continue
      end
       
      fprintf('Processing s%d m%d f%d ...\n', isubj,m,ifoi)
      
      for iblock = 1:2
        
        fprintf('Loading MEG data ...\n');
        
        try
          load(sprintf('~/pp/proc/pp_cnt_sens_cleandat_s%d_m%d_b%d_v%d.mat',isubj,m,iblock,1))
        catch me
          if v==20
            powcorr = nan(46,46);
          elseif v==12
            powcorr = nan(400,400);
          elseif v==1
            powcorr = nan(91,91);
          end
          save(sprintf([outdir 'pp_task_src_powcorr_s%d_m%d_b%d_f%d_v%d.mat'],isubj,m,iblock,ifoi,v),'powcorr');
          continue
        end
          
        dat = data.trial';
        
        if isempty(data.start_of_recording) && ~isempty(data.end_of_recording)
          if (data.end_of_recording-600*data.fsample)<1
            data.start_of_recording = 1;
          else
            data.start_of_recording = data.end_of_recording-600*data.fsample;
          end
        elseif ~isempty(data.start_of_recording) && isempty(data.end_of_recording)
          if (data.start_of_recording+600*data.fsample)>size(data.trial,2)
            data.end_of_recording = size(data.trial,2);
          else
            data.end_of_recording = data.start_of_recording+600*data.fsample;
          end
        elseif isempty(data.start_of_recording) && isempty(data.end_of_recording)
          data.start_of_recording = 5000;
          data.end_of_recording = 235000; 
        else
          data.start_of_recording = 5000;
          data.end_of_recording = 235000; 
        end
        
        dat = dat(data.start_of_recording:data.end_of_recording,:);
        dat(isnan(dat(:,1)),:)=[];
        
        pars      = [];
        pars.sa   = sprintf('~/pconn_cnt/proc/src/pconn_cnt_sa_s%d_m%d_b%d_v%d.mat',isubj,m,iblock,v_grid);
        sa        = load(pars.sa);
        
        if strcmp(allpara.grid,'cortex_lowres')
          sa.sa.grid_cortex_lowres = select_chans(sa.sa.grid_cortex3000,400);
%         elseif strcmp(allpara.grid,'cortex800')
%           sa.sa.grid_cortex800 = select_chans(sa.sa.grid_cortex3000,800);
        end
        
        pars          = [];
        pars.fsample  = 400;
        
        if strcmp(para.wavelet,'bp_filt')
          pars.segleng   = round(para.segleng.*fsample); 
          pars.segshift  = round(fsample*para.segleng/2);
        else
          pars.segleng   = round(para.segleng(ifoi).*fsample);
         	pars.segshift  = round(fsample*para.segleng(ifoi)/2);
        end
        
        if ~any(size(foi_range)==1)
          pars.foi       = foi_range(ifoi,:);
        else
          pars.foi       = foi_range(ifoi);
        end
        
        pars.epleng    = size(dat,1);
        pars.epshift   = pars.epleng;
        pars.grid      = allpara.grid;
        pars.wavelet   = para.wavelet;
        pars.scnd_filt = para.scnd_filt;
        pars.filt      = allpara.filt;
        pars.tau       = allpara.tau;
        pars.reg       = allpara.reg;
        pars.bpfreq    = para.bpfreq(ifoi,:);
        pars.weigh     = allpara.weigh;
        
        % COMPUTE POWER CORRELATIONS
        if ~timevariant
          if ~lpc
            if allpara.weigh == 0
              [powcorr] = tp_powcorr_ortho(dat,pars,sa);
            else
              [powcorr] = tp_powcorr_ortho_weight(dat,pars,sa);
            end
          else
            powcorr = tp_data2lpc_jackknife(dat,pars,filt,filt);
          end
        else
          
          pars.epleng   = para.epleng*pars.fsample;
          pars.epshift  = round(pars.epleng/16);
          pars.tau      = 'all';
          
          if ~lpc
            if pars.weigh == 0
              [powcorr] = tp_powcorr_ortho(dat,pars,sa);
            else
              [powcorr] = tp_powcorr_ortho_weight(dat,pars,sa);
            end     
          else
            powcorr = tp_lpc(dat,pars,filt,filt);
          end
        end

       save(sprintf([outdir 'pp_task_src_powcorr_s%d_m%d_b%d_f%d_v%d.mat'],isubj,m,iblock,ifoi,v),'powcorr');
        
      end
    end
  end
end
  
  
error('!')


%%

    