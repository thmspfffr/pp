%% pp_gla_src_pupil_power_correlations
% correlate bandpass filtered (or via wavelets) pupil and MEG signals

clear
restoredefaultpath

% -------------------------
% VERSION 4
% -------------------------
v = 3;
% include 28 subjects, as in pfeffer et al. (2018) plos biology
SUBJLIST  = 1:24;
freqoi    = 2.^(1:(1/4):7);
% -------------------------

addpath ~/Documents/MATLAB/fieldtrip-20160919/
addpath ~/pconn/matlab/
load(sprintf('~/pp/proc/pp_atlas_BNA.mat'))

ft_defaults

outdir = '~/pp/proc/src/';
ord    = pconn_randomization;

addpath /home/gnolte/meth/highlevel/
%%
% -------------------------
for isubj = 1:24
  
  clear data dat pupil pup dataf src_r
  
  for iblock = 1:1
    %
    fn = sprintf('pp_gla_src_pupil_power_correlations_s%d_b%d_v%d',isubj,iblock,v);
    if tp_parallel(fn,outdir,1,0)
      continue
    end
    %
    fprintf('Processing subj%d block%d ...\n',isubj,iblock);
    
    try
      % load pupil data
      if isubj < 10
        load(sprintf('~/pp/data_gla/fw4bt/osfstorage/data/gla01/pupil/sub0%d_gla_pp.mat',isubj));
      else
        load(sprintf('~/pp/data_gla/fw4bt/osfstorage/data/gla01/pupil/sub%d_gla_pp.mat',isubj));
      end
      
      pupil = data.trial{1}';
      f_sample = data.fsample;
      
      % load meg data
      if isubj < 10
        load(sprintf('~/pp/data_gla/fw4bt/osfstorage/data/gla01/meg/sub0%d_gla_meg.mat',isubj));
      else
        load(sprintf('~/pp/data_gla/fw4bt/osfstorage/data/gla01/meg/sub%d_gla_meg.mat',isubj));
      end
      
      [outp.pxx,outp.fxx]=pwelch(data.trial{1}',hanning(3200),0,0.1250:0.1250:200,400);
      
    catch me
      src_r = nan(246,25);
      save([outdir fn '.mat'],'src_r')
      continue
    end
%     len = min([size(data.trial{1},2) size(pupil,1)]);
%     data.trial{1} = data.trial{1}(:,1:len);
%     pupil = pupil(1:len,:);
    
    k = 2;
    fnq = f_sample/2;
    hil_hi = 0.005;
    hil_lo = 2;
    hil_Wn=[hil_hi/fnq hil_lo/fnq];
    [bhil, ahil] = butter(k, hil_Wn);
    
    pupil = filtfilt(bhil, ahil, pupil(:,4));
    
    pup_shift = round(f_sample*0.93); % 930s from hoeks and levelt (1992?)
    pupil = pupil(pup_shift:end); pupil(end+1:end+pup_shift-1)=nan;
    
%     data.trial{1}(:,isnan(pupil))=nan(size(data.trial{1},1),sum(isnan(pupil)));
    
%     tmp = data;
    data.avg = data.trial{1}'; %data.trial{1} = [];

    if isubj < 10
      load(sprintf('~/pp/data_gla/fw4bt/osfstorage/data/gla01/leadfields/sub0%d_gla_lf_BNA5mm.mat',isubj))
    else
      load(sprintf('~/pp/data_gla/fw4bt/osfstorage/data/gla01/leadfields/sub%d_gla_lf_BNA5mm.mat',isubj))
    end
   
                
    for ifreq=1:numel(freqoi)
      
      fprintf('Freq: %d\n',ifreq)
%       
      srate = 400;
      % freq analysis, pass 1, yields TF Fourier rep ##########################
      % - this is for multiplication with the ortho-normalised spatial filter
      tfcfg=[]; % this will be carried along ...
      tfcfg.method='wavelet';
      tfcfg.output='fourier';
      tfcfg.channel={'MEG'};
      tfcfg.foi=freqoi(ifreq);
      tfcfg.width=5.83; % again, as per Hipp et al. (2012) Nat Neurosci
      tempSD=1./(2*pi*(tfcfg.foi./tfcfg.width)); % temporal SD in sec
      tempWd=round(3*tempSD*srate)/srate; % set to 3 std dev, comp with 1000 Hz sampl rate
      tfcfg.toi=tempWd.*(1:floor(data.time{1}(end)./tempWd));
      % keep in mind that fieldtrip uses a proprietary setting of the gwidth
      % parameter (default = 3 cycles) internally that is independent of the
      % here requested time axis
      tfcfg.pad='nextpow2';
      tf=ft_freqanalysis(tfcfg,data);
      % reset freq to requested freq
      tf.freq=freqoi(ifreq);

%       mark time bins that coincide with relevant artifacts - outdated
%       artifPnts=sortrows([data.cfg.artfctdef.visual.artifact;
%                           data.cfg.artfctdef.blink.artifact]);

%       blinks removed via ICA, this is mostly myogenic, sq jumps, etc.
      artifPnts=data.cfg.artfctdef.visual.artifact;
      
      for iart = 1 : size(artifPnts,1)
        data.avg(artifPnts(iart,1):artifPnts(iart,2),:)=NaN;
      end
      
      tfPnts   =tf.time*srate;
% 
%       % only discard those bins that "center" on an artifact
%       %critDist =diff(tfPnts([1,2]),[],2)/2; % set critical dist to central 50%
      critDist=diff(tfPnts([1,2]),[],2); % set critical dist to central 100%
      keepBins=logical([]);
      % discard TF bins that overlap with artifacts (set zero in keepBins)
      for ibin=1:numel(tfPnts)
          keepBins(ibin)=~any(abs(tfPnts(ibin)-ceil(mean(artifPnts,2)))<critDist);   
      end

      % additionally omit edge datapoints to excl artif
      timeAx=tf.time; % original time axis
      timeKp=dsearchn(timeAx.',[3;timeAx(end)-3]); % excl ~ 1st & last 3 sec
      keepBins(1:timeKp(1)-1)=false;
      keepBins(timeKp(2)+1:end)=false;

      % compute cross-specral density by time bin, then average csd's
      csd=zeros(numel(tf.label)*[1,1]);
      % csd calc excludes artifact bins
      csdTime=tf.time(keepBins);
      csdData=tf.fourierspctrm(:,:,:,keepBins);
      for itbin=1:numel(csdTime)
          fspec=squeeze(csdData(:,:,:,itbin)).';
          for ichan=1:numel(tf.label);
              csd(:,ichan)=csd(:,ichan)+fspec(ichan)*conj(fspec);
          end
      end
      
      
      
      % TELL CHRISTIAN ABOUT THIS: 
      % csd used to be divided by numel(tf.time), but
      % it should rather be tf.time(keepBins)
      % ------------------
      csd=csd./numel(tf.time(keepBins)); % avg cross-spectral dens matrix
      csdData=[]; csdTime=[];
      % ------------------
      
     
      [KERNEL,f,opt]=tp_mkwavelet(freqoi(ifreq), .5, 400);
      KERNEL=repmat(KERNEL,1,size(data.avg,2));
        
      nseg=floor((size(data.avg,1)-opt.n_win)/opt.n_shift+1);
     
      clear dataf tp_csd
      kk = 0;
      
      for j=1:nseg
        
        dloc2=data.avg((j-1)*opt.n_shift+1:(j-1)*opt.n_shift+opt.n_win,:);
        
        if any(isnan(dloc2(:,1)))
          warning('NaN detected')
          dataf(:,j)=nan(size(dloc2,2),1);
          
          pup(:,j) = nan;
          continue
        else
          kk = kk + 1;
%           dataf(:,j) = transpose(sum(dloc2.*KERNEL));
          dataf(:,j)=dloc2'*KERNEL(:,1);
          tmp = pupil((j-1)*opt.n_shift+1:(j-1)*opt.n_shift+opt.n_win);
          pup(:,j) = mean(tmp.*gausswin(opt.n_win,3));
          if kk == 1
            tp_csd=dataf(:,j)*dataf(:,j)';
          else
            tp_csd=tp_csd+dataf(:,j)*dataf(:,j)';
          end
        end
      end
      
      tp_csd = tp_csd/kk;
          

%       [ss,ff,opt]=tp_mkwavelet(freqoi(ifreq), .5, 400);
%       
% %       fprintf('\ncomputing TFR...')
%       for isens = 1:size(data.avg,2)
%         tmp = conv(data.avg(:,isens),ss','same');
%         tfr(isens,:) = tmp(round(tfPnts(:)));
%         
%       end
%       
 
      outp.tp_sens_pow(:,ifreq) = mean(diag(abs(tp_csd)));
      outp.ck_sens_pow(:,ifreq) = mean(diag(abs(csd)));

      
%       all_csd(:,:,ifreq) = csd;
% %       idx_valid = find(~isnan(dataf(1,:))&~isnan(pup(1,:)));
%             
%       % -------------------------------
%       % beamforming
%       % --------------
%       para      = [];
%       para.iscs = 1;
%       para.reg  = 0.05;
%       [tp_filt,tp_pow]      = tp_beamformer(real(tp_csd),lf,para);
%       [ck_filt,ck_pow]      = tp_beamformer(real(csd),lf,para);
% % 
%       Lr = reshape(lf,[size(lf,1) 8799*3]); csd_noise = Lr*Lr';
%       [~,noise]      = tp_beamformer(real(csd_noise),lf,para);
%       outp.tp_src_nai(:,ifreq) = tp_pow./noise;
%       outp.ck_src_nai(:,ifreq) = ck_pow./noise;
%       
%       ck_dataf = squeeze(tf.fourierspctrm(:,:,:,1:end-1));
%       idx=intersect(find(keepBins),idx_valid);
%       % -------------------------------
%       % compute power
%       tp_src = abs(tp_filt'*dataf).^2;
%       ck_src = abs(ck_filt'*ck_dataf).^2;
% 
%       % correlate with pupil
% %       outp.tp_src_r(:,ifreq) = corr(pup(idx)',tp_src(:,idx)');
% %       outp.ck_src_r(:,ifreq) = corr(pup(idx)',ck_src(:,idx)');
%       outp.corr_meth_src(:,ifreq) = diag(corr(tp_src(:,idx)',ck_src(:,idx)'));
%       outp.
      clear src pup csd
      
    end
%     
%     ff=0.5 : 0.5 : 200;
%     for ifoi = 1:length(ff)
%       
%       ifoi      
%       segleng = 800;
%       segshift = 400;
%       epleng = size(data.trial{1},2);
%       gn_csd=data2cs_wavelet(data.avg,segleng,segshift,size(data.avg,1),ff(ifoi),400);
% %       [gn_dataf, gn_csd]=tp_data2cs_wavelet(data.trial{1}',segleng,segshift,epleng,ff(ifoi),f_sample);
%       
%       outp.gn_sens_pow(:,ifoi) = mean(diag(abs(gn_csd)));
%       
%     end
%       
    
    save([outdir fn '.mat'],'outp')
    tp_parallel(fn,outdir,0)
    
  end
end

error('!')
%%
clear all_corr all_pow
addpath /home/gnolte/meth/highlevel/
addpath ~/Documents/MATLAB/cbrewer/cbrewer/

SUBJLIST= 1 : 24;
v = 2;
clear all_corr
i = 0
for isubj = SUBJLIST
  isubj
  for iblock = 1 : 1
    clear src_r
    try
      load(sprintf([outdir 'pp_gla_src_pupil_power_correlations_s%d_b%d_v%d.mat'],isubj,iblock,v));
      isubj
      all_corr(:,:,isubj,iblock) = outp.src_r;
      all_pow(:,:,isubj,iblock) = outp.src_nai;
    catch me
      warning('!!!')
      all_corr(:,:,isubj,iblock) = nan(8799,25);
      all_pow(:,:,isubj,iblock) =  nan(8799,25);

      continue
    end
  end
end


all_corr= nanmean(all_corr(:,:,SUBJLIST,:),4);
all_pow= nanmean(all_pow(:,:,SUBJLIST,:),4);

for igrid = 1 : max(BNA.tissue_5mm(:))
  all_corr_BNA(igrid,:,:) = tanh(mean(atanh(all_corr(BNA.tissue_5mm == igrid,:,:))));
  all_pow_BNA(igrid,:,:) = mean(all_pow(BNA.tissue_5mm == igrid,:,:));
end

load /home/gnolte/meth/templates/mri.mat
cmap = cbrewer('div', 'RdBu', 256,'pchip'); cmap = cmap(end:-1:1,:);


%% PLOT SOURCE MAPS

ifoi = 10;

para.colorlimits = [-0.1 0.1]
para.colormaps{1} = jet;
para.orientation = 'axial';
para.dslice_shown = 0.5;

showmri_transp(mri,para,[BNA.grid_5mm./10 nanmean(all_corr(:,ifoi,:),3)])


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

print(gcf,'-dpdf',sprintf('~/pp/plots/pp_gla_src_pupil_power_correlations_v%d.pdf',v))

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
tp_editplots; title('Average over regions')
axis([0 25 -0.075 0.075])
% subplot(1,2,2); hold on
% plot(t,'linewidth',3);
% line([1 25],[0 0],'linestyle',':','color','k')
% set(gca,'xtick',[1 5 9 13 17 21 25],'xticklabel',num2cell(freqoi([1 5 9 13 17 21 25])))
% xlabel('Frequency [Hz]'); ylabel('T-Value')
% tp_editplots

print(gcf,'-dpdf',sprintf('~/pp/plots/pp_gla_src_pupil_power_correlations_lineplot_v%d.pdf',v))

%% 
ff=0.125 : 0.125 : 150;
p=zscore(nanmean(outp.pxx,2));

for ifoi = 1 : length(freqoi)
[~,f,opt]=tp_mkwavelet(freqoi(ifoi), .5, 400);


pow(ifoi) = mean(p(outp.fxx>f(1) & outp.fxx<f(2)),1);

end

plot(pow); hold on