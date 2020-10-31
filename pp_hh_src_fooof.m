%% pp_hh_src_fooof
% correlate bandpass filtered (or via wavelets) pupil and MEG signals

clear
restoredefaultpath

% -------------------------
% VERSION 4
% -------------------------
v = 3;
% include 28 subjects, as in pfeffer et al. (2018) plos biology
SUBJLIST = [4 5 6 7 8 9 10 11 12 13 15 16 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34];
% -------------------------

addpath ~/Documents/MATLAB/fieldtrip-20160919/
addpath ~/pconn/matlab/
load(sprintf('~/pp/proc/pp_atlas_BNA.mat'))

ft_defaults

outdir = '~/pp/proc/src/';
ord    = pconn_randomization;

freqoi=2.^(1:(1/4):7); % 2-128 Hz as per Hipp et al. (2012) Nat Neurosci

%%
% -------------------------
for isubj = SUBJLIST
    
    % identify placebo condition (ord==1)
    im = find(ord(isubj,:)==1);
    
    for iblock = 1:2
        %
        fn = sprintf('pp_hh_src_fooof_s%d_b%d_v%d',isubj,iblock,v);
        if tp_parallel(fn,outdir,1,0)
            continue
        end
        % %
        fprintf('Processing subj%d block%d ...\n',isubj,iblock);
        
        try
            % load cleaned meg data
            load(sprintf('~/pupmod/proc/sens/pupmod_rest_sens_cleandat_s%d_m%d_b%d_v%d.mat',isubj,im,iblock,1))
            % load cleanted pupil data
            load(sprintf('~/pp/proc/pup/pp_pupil_diameter_cleaned_s%d_m%d_b%d.mat',isubj,im,iblock))
        catch me
            continue
        end
        
        % bp-filter and resample pupil
        % ------
        k = 2; f_sample = 1000;
        fnq = f_sample/2;
        hil_hi = 0.005; hil_lo = 2;
        hil_Wn=[hil_hi/fnq hil_lo/fnq];
        [bhil, ahil] = butter(k, hil_Wn);
        
        pupil = filtfilt(bhil, ahil, pupil(:,4));
        pupil = resample(pupil,400,1000);
        % ------
        
        % align pupil and meg (at signal offset)
        % ------
        if size(pupil,2)>3
            pupil = pupil(end:-1:1,4);
        else
            pupil = pupil(end:-1:1);
        end
        
        %     data.trial = data.trial(:,1:data.end_of_recording);
        dat = dat(:,end:-1:1);
        
        len = min([size(pupil,1) size(dat,2)]);
        if len/400 > 600
            len = 400*600;
        end
        
        dat = dat(:,1:len);
        pupil = pupil(1:len);
        
        dat = dat(:,end:-1:1);
        pupil = pupil(end:-1:1);
        % ------
        
        % pupil shift: 930 ms from hoeks & levelt (1992)
        pup_shift = round(400*0.93);
        pupil = pupil(pup_shift:end); pupil(end+1:end+pup_shift-1)=nan;
        
        
        
        dat(:,isnan(pupil))=nan(size(dat,1),sum(isnan(pupil)));
        
        load(['/home/tpfeffer/pp/proc/src/' sprintf('pp_sa_s%d_m%d_b%d_v%d.mat',isubj,im,iblock,1)],'sa');
        
        lf = sa.L_genemaps_aal;
        
        [outp.pxx,outp.fxx]=pwelch(dat(:,~isnan(dat(1,:)))',hanning(400),0,1:1:200,400);
        
        clear csd
        
        for ifreq=1:length(freqoi)
            ifreq
            
            para          = [];
            para.freq     = freqoi(ifreq);
            para.fsample  = 400;
            para.overlap  = 0.5;
            [csd(:,:,ifreq)]=tp_compute_csd_wavelets(dat,para);
            
        end
        
        csd = nanmean(csd,3);
        
        para          = [];
        para.reg      = 0.05;
        filt = tp_beamformer(real(csd),sa.L_genemaps_aal,para);
        % --------------
        
        dat_src = dat'*filt;
        
        idx1 = isnan(dat_src(:,1));
        idx2 = isnan(pupil);
        
        dat_src(idx1 | idx2, : ) = nan;
        pupil(idx1 | idx2) = nan;
        pupil_df = diff(pupil);
        
        opt.n_win = 800; % 10s segment length, i.e., 0.1:0.1:100
        opt.n_shift = 800; % no overlap
        
        nseg=floor((size(dat_src,1)-opt.n_win)/opt.n_shift+1);
        clear pxx fxx pup pup_df
        ff = 2:0.5:40;
        
        pxx = nan(size(ff,2),size(filt,2),nseg);
        for iseg = 1 : nseg
            
            fprintf('%d / %d\n',iseg,nseg)
            seg_dat = dat_src((iseg-1)*opt.n_shift+1:(iseg-1)*opt.n_shift+opt.n_win,:);

            if any(isnan(seg_dat(:,1)))
                pup(iseg) = nan;
                pup_df(iseg)=nan;
                continue        
            end
        
            [pxx(:,:,iseg),fxx]=pwelch(seg_dat,hanning(opt.n_win),[],ff,400,'power');
            pup(iseg)  = nanmean(pupil((iseg-1)*opt.n_shift+1:(iseg-1)*opt.n_shift+opt.n_win));
            pup_df(iseg) = nanmean(pupil_df((iseg-1)*opt.n_shift+1:(iseg-1)*opt.n_shift+opt.n_win));
        end
  
        save([outdir fn '.mat'],'pxx','fxx','pup','pup_df')
        tp_parallel(fn,outdir,0)
        
        clear src_r all_nai outp
    end
end


error('!')

