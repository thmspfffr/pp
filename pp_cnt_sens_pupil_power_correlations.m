%% pp_sens_pupil_power_correlations
% correlate bandpass filtered (or via wavelets) pupil and MEG signals

clear
SUBJLIST  = [4 5 6 7 8 9 10 11 12 13 15 16 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34];
v = 1;
addpath ~/Documents/MATLAB/fieldtrip-20160919/
ft_defaults

outdir = '~/pp/proc/sens/';
addpath ~/pconn/matlab/
ord = pconn_randomization;

    freqoi=2.^(1:(1/4):7); % 1-128 Hz as per Hipp et al. (2012) Nat Neurosci
%%
for isubj = SUBJLIST
  
  im = find(ord(isubj,:)==1);
  
  for iblock = 1:2
% %     
    fn = sprintf('pp_cnt_sens_pupil_power_correlations_s%d_b%d_v%d',isubj,iblock,v);
    if tp_parallel(fn,outdir,1,0)
      continue
    end
    
    fprintf('Processing subj%d block%d ...\n',isubj,iblock);
    
    try
      % load meg data
      load(sprintf('~/pp/proc/pp_cnt_sens_cleandat_s%d_m%d_b%d_v%d.mat',isubj,im,iblock,v))
      % load pupil during counting
      load(sprintf('~/pp/proc/pup/pp_pupil_diameter_cleaned_counting_s%d_m%d_b%d.mat',isubj,im,iblock))
      
    catch me
      r = nan(274,25);
      save([outdir fn '.mat'],'r')
      continue
    end
%     
    k = 2;
    f_sample = 1000;
    fnq = f_sample/2;
    hil_hi = 0.005;
    hil_lo = 2;
    hil_Wn=[hil_hi/fnq hil_lo/fnq];
    [bhil, ahil] = butter(k, hil_Wn);
    
    pupil = filtfilt(bhil, ahil, pupil(:,4));
    
    pupil = resample(pupil,400,1000);
    
    data.trial = data.trial(:,1:data.end_of_recording);
    data.trial = data.trial(:,end:-1:1);
    pupil = pupil(end:-1:1);
    
    len = min([size(pupil,1) size(data.trial,2)]);
    if len/400 > 600
      len = 400*600;
    end
    
    data.trial = data.trial(:,1:len);
    pupil = pupil(1:len);
    
    data.trial = data.trial(:,end:-1:1);
    pupil = pupil(end:-1:1);
    
    % replace this later
    pupil(isnan(data.trial(1,:)))=[];
    data.trial(:,isnan(data.trial(1,:))) = [];

    data.avg = data.trial; data.trial = [];
    data.dimord = 'chan_time';
    data.time = [];
    data.time =  1/data.fsample:1/data.fsample:size(data.avg,2)/data.fsample;
    
    %%% 1.3 loop over frequencies & source reconstruct
    srate    = 400;
    for ifreq=1:numel(freqoi)
      
      cfg=[];
      cfg.method  ='wavelet';
      cfg.output  = 'fourier';
      cfg.channel = {'MEG'};
      cfg.foi     = freqoi(ifreq);
      cfg.width   = 5.83; % again, as per Hipp et al. (2012) Nat Neurosci
      tempSD      = 1./(2*pi*(cfg.foi./cfg.width)); % temporal SD in sec
      tempWd      = round(3*tempSD*srate)/srate; % set to 3 std dev, comp with 1000 Hz sampl rate
      cfg.toi     = tempWd.*(1:floor(data.time(end)./tempWd));
 
      cfg.pad='nextpow2';
      tf=ft_freqanalysis(cfg,data);
      
      % Compute mean pupil diameter for time segments
      clear pup pow
      for itoi = 1 : size(cfg.toi,2)-1
        t = cfg.toi(itoi);
        idx = find(round(data.time.*10000)==round(t*10000));
        pup(itoi) = mean(pupil(idx-(tempWd*srate-1)/2:idx+(tempWd*srate-1)/2));
      end
      % ---------------
      
      pow = abs(squeeze(tf.fourierspctrm(1,:,1,:)).^2);
      nanidx = isnan(pow(1,:));
      nanidx(end) = [];
      pup(nanidx) = [];
      pow(:,nanidx) = [];
      pow(:,end)=[];
      
      [tmp_r ,p]=corr(pup',pow');
      
      if isubj >= 32
        load(sprintf('~/pconn/proc/preproc/pconn_chanidx_s%d_m%d.mat',31,3))
      elseif isubj < 4
        load(sprintf('~/pconn/proc/preproc/pconn_chanidx_s%d_m%d.mat',4,1))
      elseif isubj == 17
        load(sprintf('~/pconn/proc/preproc/pconn_chanidx_s%d_m%d.mat',16,1))
      else
        load(sprintf('~/pconn/proc/preproc/pconn_chanidx_s%d_m%d.mat',isubj,im))
      end

      %
      r(:,ifreq) = pconn_sens_interp274(idx,tmp_r);
      %
      %
    end
    
%     plot(r); drawnow; hold on
    
    save([outdir fn '.mat'],'r')
    tp_parallel(fn,outdir,0)
    
  end
end

error('!')
%%
clear all_corr_cnt
clear all_corr
%
v = 1
for isubj = SUBJLIST
  for iblock = 1 : 2

    load(sprintf([outdir 'pp_cnt_sens_pupil_power_correlations_s%d_b%d_v%d.mat'],isubj,iblock,v));

    all_corr_cnt(:,:,isubj,iblock) = r(:,1:25); clear r
        load(sprintf([outdir 'pp_sens_pupil_power_correlations_s%d_b%d_v%d.mat'],isubj,iblock,v));

all_corr(:,:,isubj,iblock) = r(:,1:25); 
  end
end

all_corr= nanmean(all_corr(:,:,SUBJLIST,:),4);
all_corr_cnt= nanmean(all_corr_cnt(:,:,SUBJLIST,:),4);



%%
figure; set(gcf,'color','w');
freqoi=2.^(1:(1/4):7); % 1-128 Hz as per Hipp et al. (2012) Nat Neurosci
subplot(1,3,1);
imagesc(nanmean(all_corr,3),[-0.05 0.05]);
xlabel('Frequency [Hz]'); ylabel('Sensor #')
set(gca,'xtick',[1 5 9 13 17 21 25],'xticklabel',num2cell(freqoi([1 5 9 13 17 21 25])))
title('Fixation')

subplot(1,3,2);
imagesc(nanmean(all_corr_cnt,3),[-0.05 0.05]); 
xlabel('Frequency [Hz]'); ylabel('Sensor #')
set(gca,'xtick',[1 5 9 13 17 21 25],'xticklabel',num2cell(freqoi([1 5 9 13 17 21 25])))
title('Counting')

[h,p,~,s] = ttest(all_corr,all_corr_cnt,'dim',3);
h = p<fdr1(p(:),0.05);

subplot(1,3,3);
imagesc((nanmean(all_corr_cnt-all_corr,3)).*h,[-0.05 0.05]); 
xlabel('Frequency [Hz]'); ylabel('Sensor #')
set(gca,'xtick',[1 5 9 13 17 21 25],'xticklabel',num2cell(freqoi([1 5 9 13 17 21 25])))
title('Difference')

% all_corr = nanmean(all_corr(:,:,SUBJLIST,:),4);
%% 
cmap = cbrewer('div','RdBu',128); cmap = cmap(end:-1:1,:)
ifoi = 10;

load /home/tpfeffer/pconn/proc/src/pconn_sa_s26_m2_b1_v4.mat
pars            = [];
pars.markersize = 0;
pars.linewidth  = 3;
pars.cbar       = 0;
pars.scale      = [-0.05 0.05]
pars.cmap       = cmap;
pars.resolution = 600;

figure; set (gcf,'color','w')

[~,p,~,s]=ttest(zeros(274,28),squeeze(all_corr(:,ifoi,:)),'dim',2)

subplot(1,2,1);
showfield_colormap(nanmean(all_corr(:,ifoi,:),3),sa.locs_2D,pars);
subplot(1,2,2);
pars.scale      = [-5 5]
showfield_colormap(-s.tstat,sa.locs_2D,pars);

drawnow

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
