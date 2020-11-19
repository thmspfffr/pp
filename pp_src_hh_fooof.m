%  out = funchier_computegranger(dat,sa,para)

% Computes Granger causality as described in Dhamala et al. (2008).
% The function first computes X using fieldtrips inbuilt function.
% Then it esitimates directed interactions (granger causality) from the
% cross spectra. The measure can be computed on band-limited or broadband
% signals. The latter is the default setting, as this is used in 
% Bastos et al. (2015) Neuron and Michalareas et al. (2016) Neuron.
% -------------------------
% thms.pfffr@gmail.com, 10/2020 
% -------------------------

clear; restoredefaultpath
addpath ~/pconn/matlab

% -------------------------
% VERSION 
% -------------------------
v = 1;
% include 28 subjects, as in pfeffer et al. (2018) plos biology
SUBJLIST = [4 5 6 7 8 9 10 11 12 13 15 16 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34];
% -------------------------


outdir = '~/funchier/proc/src/';
ord    = pconn_randomization;
fs = 400;
addpath ~/Documents/MATLAB/codes/

%%
% -------------------------
for isubj = SUBJLIST
    
    % identify placebo condition (ord==1)
    im = find(ord(isubj,:)==1);
    
    for iblock = 1:2
        %
        fn = sprintf('pp_src_hh_fooof_s%d_b%d_v%d',isubj,iblock,v);
        if tp_parallel(fn,outdir,1,0)
            continue
        end
        % %
        fprintf('Processing subj%d block%d ...\n',isubj,iblock);
        
        try
            % load cleaned meg data
            load(sprintf('~/pupmod/proc/sens/pupmod_rest_sens_cleandat_s%d_m%d_b%d_v%d.mat',isubj,im,iblock,1))
            % load cleanted pupil data
%             load(sprintf('~/pp/proc/pup/pp_pupil_diameter_cleaned_s%d_m%d_b%d.mat',isubj,im,iblock))
        catch me
            continue
        end
        
        dat(:,isnan(dat(1,:)))=[];
        
        load(['/home/tpfeffer/funchier/proc/headmodel/' sprintf('funchier_headmodel_s%d_m%d_b%d_v%d.mat',isubj,im,iblock,1)]);
        
        csd_sens=data2cs_event(dat',400,400,size(dat,2),150);
        csd_sens=mean(real(csd_sens(:,:,5:end)),3);
        para          = [];
        para.reg      = 0.05;
        [filt,pow] = tp_beamformer(csd_sens,sa.L_vtpm_4mm,para);
        
        dat_src = dat'*filt;
        

        
        
        save([outdir fn '.mat'],'M')
        tp_parallel(fn,outdir,0)
    end
    
end
  
error('!')

%% LOAD AND PLOT EVERYTHING

for isubj = SUBJLIST
    isubj
%     for m = 1 : 3
%         im = find(ord(isubj,:)==m);
        for iblock = 1 : 2
            
            try
            load(sprintf('/home/tpfeffer/funchier/proc/src/funchier_computegranger_s%d_b%d_v1.mat',isubj,iblock))
            
            cgc_all(:,:,:,isubj,iblock) = (M.cgc_l + M.cgc_r)./2;
            
            catch me
                cgc_all(:,:,:,isubj,iblock) = nan(141,23,23);
                
                end
            
        end
%     end
end

cgc_all = nanmean(cgc_all(:,:,:,SUBJLIST,:),5);
     
cgc_all(:,:,:,[9 10 15 16 19]) = [];
9,10,15,16, 19
    
