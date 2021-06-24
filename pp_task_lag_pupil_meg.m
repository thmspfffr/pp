%% READ IN DATA AND SAVE AS MAT
% pconn_preproc_read_data

clear all
restoredefaultpath

% -------------------------------------------------------------------------
% VERSION 4 - TASK
% -------------------------------------------------------------------------
denoise   = 'no';
v         = 1;
pad       = 1200/2; % 500 ms
rest      = 0;
task      = 1;
% -------------------------------------------------------------------------

indir1  = '/home/tpfeffer/pconn/rawdata/meg/';
outdir  = '/home/tpfeffer/pconn_cnt/proc/preproc/';

addpath /home/tpfeffer/pconn_cnt/matlab/
addpath('/home/tpfeffer/Documents/MATLAB/fieldtrip-20160919/')

% addpath('/home/tpfeffer/Documents/MATLAB/fieldtrip-20130925/')

ft_defaults

%%
for m = 2
  for isubj =4
  
  indir = sprintf([indir1 'p%d/s%d/'],isubj,m);
  
%   try
    
%     if ~exist([outdir sprintf('pconn_cnt_preproc_s%d_m%d_v%d_processing.txt',isubj,m,v)]) 
% %         && ~exist([outdir sprintf('pconn_preproc_data_s%d_m%d_v%d.mat',isubj,m,v)])
%       
%       system(['touch ' outdir sprintf('pconn_cnt_preproc_s%d_m%d_v%d_processing.txt',isubj,m,v)]);
%       disp(sprintf('Processing subject s%d m%d',isubj,m));
%       
%       % no proc file, but output: create proc file
%     elseif ~exist([outdir sprintf('pconn_cnt_preproc_s%d_m%d_v%d_processing.txt',isubj,m,v)]) ...
%         && exist([outdir sprintf('pconn_cnt_preproc_data_s%d_m%d_v%d.mat',isubj,m,v)])
%       
%       system(['touch ' outdir sprintf('pconn_cnt_preproc_s%d_m%d_v%d_processing.txt',isubj,m,v)]);
%       continue
%       
%     elseif exist([outdir sprintf('pconn_cnt_preproc_s%d_m%d_v%d_processing.txt',isubj,m,v)]) ...
%         && ~exist([outdir sprintf('pconn_cnt_preproc_data_s%d_m%d_v%d.mat',isubj,m,v)])
%       
%       disp(sprintf('Processing subject s%d m%d',isubj,m));
%       
%     else
%     continue
%     end
    
    cont_dir = dir(indir);
    
    for i = 1 : length(cont_dir)
      tmp(i) = ~isempty(regexp(cont_dir(i).name,'(.ds)','once'));
    end
    
    ind = find(tmp>0);
    
    ibl(1:3) = 0;
    
    % no button presses registered for subject 10, m = 3
    if isubj == 10 && m == 3
      ind = [ind(3) ind(6)];
      cond = 3;
%     elseif isubj == 22
    end
    
    for idir = 3:8
      
      tmp_dataset      = [cont_dir(idir).name];
     
%       if isubj == 22 && m == 3
%       	cfgs.datafile    = [indir '/' tmp_dataset '/' sprintf('NK9_PConn_20141202_0%d.meg4',idir-2)];
%         cfgs.headerfile  = [indir '/' tmp_dataset '/' sprintf('NK9_PConn_20141202_0%d.res4',idir-2)];
%       else
        cfgs.datafile    = [indir '/' tmp_dataset '/' cont_dir(idir).name(1:end-2) 'meg4'];
        cfgs.headerfile  = [indir '/' tmp_dataset '/' cont_dir(idir).name(1:end-2) 'res4'];
%       end
      
      cfg = [];
      % cfg.continuous = 'yes';
      cfg.dataset    = [indir cont_dir(idir).name];
      
      if isubj == 22 && m == 3
        a=ft_read_event(cfgs.headerfile);
      else
        a=ft_read_event(cfg.dataset);
      end
      
      if ~(isubj ==10 && m==3)
        if size(cell2mat({a.value}),2)<10
            warning('Rest!');
            cond = 1;
            continue
        elseif sum(cell2mat({a.value})>50)>7
            warning('Task - button press!');
            cond = 2;
            continue
        elseif size(cell2mat({a.value}),2) > 10 && sum(cell2mat({a.value})>50)<7
            warning('Task - counting!');
            cond = 3;
  %           continue
        end
      end
      
      ibl(cond) = ibl(cond) + 1;

      ibl(cond)
%       cfg = [];
      cfg.channel   = {'EEG001';'EEG002'};
      cfg.continuous = 'yes'
      cfg.demean    = 'yes';
      cfg.precision = 'single';
      cfg.bpfreq = [1 50];
      cfg.bpfilter = 'yes';
      data          = ft_preprocessing(cfg);
      
      cfg = [];
      cfg.resamplefs = 1000;
      data = ft_resampledata(cfg,data);
      
      load ~/pp/proc/pup/pp_pupil_diameter_cleaned_counting_s4_m2_b1.mat

%       [r,lags]=xcorr(fillmissing(pupil(1:610000,2),'spline'),data.trial{1}(1,1:610000)','normalized');
      [r,lags]=xcorr(fillmissing(pupil(1:610000,3),'spline'),data.trial{1}(2,1:610000)','normalized');
%       [r,lags]=xcorr(pupil(1:610000,3),data.trial{1}(1,1:610000),'normalized');
      
      % EOG channels
%       data.trial{1}(304,:),data.trial{1}(305,:))
      
      


    end
%   catch me
%     sprintf('ERROR!');
%     save([outdir sprintf('pconn_preproc_data_s%d_m%d_v%d_err.mat',isubj,m,v)], 'me')
%   end
end
end