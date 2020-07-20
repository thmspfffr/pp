%% pp_cnt_src_pupil_power_correlations
% correlate bandpass filtered (or via wavelets) pupil and MEG signals

clear
restoredefaultpath

SUBJLIST  = [4 5 6 7 8 9 10 11 12 13 15 16 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34];

% focus on all pupil
v = 1; % with guassian applied to pupil segments 

addpath ~/Documents/MATLAB/fieldtrip-20160919/
ft_defaults

outdir = '~/pp/proc/src/';
addpath ~/pconn/matlab/
ord = pconn_randomization;
load(sprintf('~/pp/proc/pp_atlas_BNA.mat'))
%%

for isubj = SUBJLIST
  
  im = find(ord(isubj,:)==1);
  
  for iblock = 1:2
%     
    fn = sprintf('pp_cnt_src_pupil_power_correlations_s%d_b%d_v%d',isubj,iblock,v);
    if tp_parallel(fn,outdir,1,0)
      continue
    end
    
    fprintf('Processing subj%d block%d ...\n',isubj,iblock);
    
    try
      % load meg data
      load(sprintf('~/pp/proc/pp_cnt_sens_cleandat_s%d_m%d_b%d_v%d.mat',isubj,im,iblock,1))
      % load pupil data
      load(sprintf('~/pp/proc/pup/pp_pupil_diameter_cleaned_counting_s%d_m%d_b%d.mat',isubj,im,iblock))
      
    catch me
      sprintf('%s',me.message)
      r = nan(246,25);
      save([outdir fn '.mat'],'src_r')
      continue
    end
    
    k = 2;
    f_sample = 1000;
    fnq = f_sample/2;
    hil_hi = 0.005;
    hil_lo = 2;
    hil_Wn=[hil_hi/fnq hil_lo/fnq];
    [bhil, ahil] = butter(k, hil_Wn);
    
    pupil = filtfilt(bhil, ahil, pupil(:,4));

    pupil = resample(pupil,400,1000);
    
    if size(pupil,2)>3
      pupil = pupil(end:-1:1,4);
    else
      pupil = pupil(end:-1:1);
    end
    
    data.trial = data.trial(:,1:data.end_of_recording);
    data.trial = data.trial(:,end:-1:1);
    
    len = min([size(pupil,1) size(data.trial,2)]);
    if len/400 > 600
      len = 400*600;
    end
    
    data.trial = data.trial(:,1:len);
    pupil = pupil(1:len);
    
    data.trial = data.trial(:,end:-1:1);
    pupil = pupil(end:-1:1);
    
    pup_shift = round(400*0.93); % 930s from hoeks and levelt (1992?)    
    pupil = pupil(pup_shift:end); pupil(end+1:end+pup_shift-1)=nan;
      
    data.trial(:,isnan(pupil))=nan(size(data.trial,1),sum(isnan(pupil)));
    
    % ----------------
    % THIS NEEDS TO BE CHANGED
    % ----------------
%     pupil(isnan(data.trial(1,:)))=[];
%     data.trial(:,isnan(data.trial(1,:))) = [];
    % ----------------
    % ----------------
    % ----------------
    % ----------------
    % ----------------
    
    data.avg = data.trial; data.trial = [];
    data.dimord = 'chan_time';
%     data.time = [];
    data.time =  1/data.fsample:1/data.fsample:size(data.avg,2)/data.fsample;
    
    %%% 1.3 loop over frequencies & source reconstruct
    freqoi=2.^(1:(1/4):7); % 1-128 Hz as per Hipp et al. (2012) Nat Neurosci
    srate    = 400;
    
    load(['/home/tpfeffer/pp/proc/src/' sprintf('pp_sa_s%d_m%d_b%d_v%d.mat',isubj,im,iblock,1)],'sa');

    for ifreq=1:numel(freqoi)
      ifreq
       octave = 0.5; % Frequency resolution
      % arithmetic mean
      foi_min    = 2*freqoi(ifreq)/(2^octave+1);
      foi_max    = 2*freqoi(ifreq)/(2^-octave+1);
      delta_freq = foi_max-foi_min; % 2*std in freq domain
      delta_time = 6/pi./delta_freq;
      delta_time = round(delta_time*1000)/1000;
      t_shift    = delta_time/2;
      n_win      = round(delta_time*data.fsample);
      n_shift    = round(t_shift*data.fsample);
      TAPER      = gausswin(n_win,3)'; TAPER = TAPER/sum(TAPER);
      iEXP       = exp(sqrt(-1) * ((1:n_win)-n_win/2-0.5) /data.fsample*freqoi(ifreq)*2*pi);
      KERNEL     = (TAPER.*iEXP).';
      
      nseg=floor((size(data.avg,2)-n_win)/n_shift+1);
      
      %
      for j=1:nseg
        
        dloc2=data.avg(:,(j-1)*n_shift+1:(j-1)*n_shift+n_win)';
        
        if any(isnan(dloc2(:,1)))
          warning('NaN detected')
          dataf(:,j)=nan(1,size(dloc2,2));
          pup(:,j) = nan;
          continue
        else
          dataf(:,j)=dloc2'*KERNEL;
          tmp = pupil((j-1)*n_shift+1:(j-1)*n_shift+n_win);
          pup(:,j) = mean(tmp.*gausswin(n_win,3));
        end
      end
      
      idx_valid = find(~isnan(dataf(1,:)));
      csd = (dataf(:,idx_valid)*dataf(:,idx_valid)')/size(dataf(:,idx_valid ),2);
                
      para      = [];
      v         = 1;
      para.iscs = 1;
      para.reg  = 0.05;
      
      % NOTE THAT THE LEADFIELD HAS A BAD NAME DUE TO A FALSE ASSIGNMENT IN
      % THE FUNCTION THAT CREATES THE LEADFIELD. THIS IS THE CORRECT
      % LEADFILED THOUGH
      % -------------------------------
      filt      = tp_beamformer(csd,sa.L_genemaps_aal,para);
      % -------------------------------
      
      % compute power
      src = abs((filt'*dataf).^2);
      
      % average across all vertices within a BNA region
      for igrid = 1 : max(BNA.tissue_5mm(:))
        src_pow = mean(src(BNA.tissue_5mm == igrid,:));
        src_r(igrid,ifreq) = corr(pup(idx_valid)',src_pow(idx_valid)');
      end

       clear src pup
      
    end

    save([outdir fn '.mat'],'src_r')
    tp_parallel(fn,outdir,0)
    
  end
end

error('!')

%%
v = 1;
clear all_corr
i = 0
for isubj = SUBJLIST
  for iblock = 1 : 2
     clear src_r
     try
    load(sprintf([outdir 'pp_src_pupil_power_correlations_s%d_b%d_v%d.mat'],isubj,iblock,v));
    isubj
    all_corr(:,:,isubj,iblock) = src_r;
    if ~exist('src_r','var')
      src_r = r;
      save(sprintf([outdir 'pp_src_pupil_power_correlations_s%d_b%d_v%d.mat'],isubj,iblock,v),'src_r');
    end
     catch me
       all_corr(:,:,isubj,iblock) = nan(246,25);
     continue
     end
%     all_corr(:,:,isubj,iblock) = r;

  end
end

all_corr= nanmean(all_corr(:,:,SUBJLIST,:),4);
%
%% STATISTICS (simple t-test)
addpath ~/Documents/MATLAB/cbrewer/cbrewer/
cmap = cbrewer('div','RdBu',128); cmap = cmap(end:-1:1,:)

[h,p]=ttest(zeros(size(all_corr)),all_corr,'dim',3)

figure; set(gcf,'color','w')
h=subplot(2,2,[1 3]); hold on
imagesc(nanmean(all_corr,3),[-0.05 0.05])
set(gca,'xtick',[1 5 9 13 17 21 25],'xticklabel',num2cell(freqoi([1 5 9 13 17 21 25])),'ydir','reverse')
xlabel('Frequency [Hz]'); ylabel('BNA Region')
colormap(cmap); colorbar; axis([1 25 1 246]); tp_editplots

h=subplot(2,2,[2 4]); hold on
imagesc(nanmean(all_corr,3).*(p<0.0157),[-0.05 0.05])
set(gca,'xtick',[1 5 9 13 17 21 25],'xticklabel',num2cell(freqoi([1 5 9 13 17 21 25])),'ydir','reverse')
xlabel('Frequency [Hz]'); ylabel('BNA Region')
colormap(cmap); colorbar; axis([1 25 1 246]); tp_editplots

print(gcf,'-dpdf',sprintf('~/pp/plots/pp_src_pupil_power_correlations_v%d.pdf',v))
%% PROJECT TO SURFACE

par = squeeze(nanmean(all_corr(:,freqoi>8&freqoi<18,:),2));
par = nanmean(par,2);
close all

cmap = cbrewer('div','RdBu',256); cmap = cmap(end:-1:1,:)


para = [];
para.clim = [-0.02 0.02];
para.cmap = cmap;
para.grid = BNA.centroids;
para.dd = 0.5;
para.fn = sprintf('~/pp/plots/pp_surface_v%d.png',v);
tp_plot_surface(par,para)






%%

for ifoi = 1 : 25
  
m_f(ifoi) = mean(squeeze(nanmean(all_corr(:,ifoi,:),1)));
s_f(ifoi) = nanstd(squeeze(nanmean(all_corr(:,ifoi,:),1)))/sqrt(size(all_corr,3));
[h(ifoi),p(ifoi)] = ttest(squeeze(nanmean(all_corr(:,ifoi,:),1)));

end



figure; set (gcf,'color','w')
subplot(2,2,1); hold on
shadedErrorBar([],m_f,[s_f],[]);
% plot(f,'linewidth',3);
line([1 25],[0 0],'linestyle',':','color','k')
set(gca,'xtick',[1 5 9 13 17 21 25],'xticklabel',num2cell(freqoi([1 5 9 13 17 21 25])))
xlabel('Frequency [Hz]'); ylabel('Correlation coeff.')
tp_editplots; title('Average over sensors')
axis([0 25 -0.04 0.04])
% subplot(1,2,2); hold on
% plot(t,'linewidth',3);
% line([1 25],[0 0],'linestyle',':','color','k')
% set(gca,'xtick',[1 5 9 13 17 21 25],'xticklabel',num2cell(freqoi([1 5 9 13 17 21 25])))
% xlabel('Frequency [Hz]'); ylabel('T-Value')
% tp_editplots

print(gcf,'-dpdf',sprintf('~/pp/plots/pp_sens_pupil_power_correlations_lineplot_v%d.pdf',v))

