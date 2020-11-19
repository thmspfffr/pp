clear

% SUBJECT10, M3, B1: no triggers

SUBJLIST  = [4 5 6 7 8 9 10 11 12 13 15 16 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34];
for isubj = SUBJLIST
  for m = 1:3
    for iblock = 2

      %       load(sprintf('/home/tpfeffer/pconn/proc/preproc/pconn_postpostproc_s%d_m%d_b%d_v1.mat',isubj,m,iblock)
      try
        load(sprintf('/home/tpfeffer/pconn/proc/preproc/pconn_trig_raw_s%d_m%d_b%d_v1.mat',isubj,m,iblock))
%       load(sprintf('/home/tpfeffer/pconn/proc/preproc/pconn_preproc_trig_proc_s%d_m%d_b%d_v1.mat',isubj,m,iblock))
      catch me
        warning(sprintf('No triggers for s%d, m%d, b%d',isubj,m,iblock))
      end
        
      load(sprintf('~/pconn/proc/preproc/pconn_preproc_artvec_s%d_m%d_b%d_v1.mat',isubj,m,iblokc));
      
      START REALIGNING MEG SIGNAL
      ----------------------------------
      art = round(art./3);
      
%       cfg.resamplefs = 400;
%       [trig] = ft_resampledata(cfg, trig);
      
      artvec = zeros(size(trig.trial{1}));
      for iart = 1 : size(art,1)
        artvec(art(iart,1):art(iart,2))=1;
      end
      bw = bwlabel(artvec);
      
      if (size(data.trial{1},2)+sum(bw>0))~=size(data.idx(:),1)
        difference = size(data.trial{1},2)+sum(bw>0)-size(data.idx(:),1);
      
        if difference > 0 && difference <= 10
          for idiff = 1 : difference
            bw(find(bw>0,1,'first'))=0;
          end
      
        elseif difference < 0
          for idiff = 1 : abs(difference)
            bw(find(bw==1,1,'last')+1)=1;
          end
      
        else
          error('wtf')
        end
      end
      
      clear idx
      dat_new = []
      data_old = data.trial{1}(1,:);
      for ibw = 1 : max(bw)
      
        if ibw == 1
          dat_len(1,:) = length(1:find(bw==ibw,1,'first')-1);
        else
          dat_len(ibw,:) = length(find(bw==ibw-1,1,'last')+1:find(bw==ibw,1,'first')-1);
        end
        art_len(ibw,:) = sum(bw==ibw);
      
      end
      
      dat_new = [];
      data_old = data.trial{1}(1,:);
      
      for i = 1 : size(dat_len,1)
        dat1 = data_old(1,1:dat_len(i));
        data_old(1:dat_len(i)) = [];
        dat_new = [dat_new dat1 nan(1,art_len(i))];
      end
      
      % FIND TRIGGERS IN DATA
      trig_idx = bwlabel(abs(trig.trial{1})>0.5);
      if max(trig_idx)==1
        if round(mean(find(trig_idx==1))) > 200000
          off = round(mean(find(trig_idx==1)));
          fprintf('S%dM%dB%d: Only offset trigger\n',isubj,m,iblock)
        else
          start = round(mean(find(trig_idx==1)));
          fprintf('S%dM%dB%d: Only start trigger\n',isubj,m,iblock)
        end
      elseif max(trig_idx)==2
        start = round(mean(find(trig_idx==1)));
        off = round(mean(find(trig_idx==2)));
        
        fprintf('S%dM%dB%d: Both triggers, diff = %d\n',isubj,m,iblock,round((off-start)./1200))
        
        if round((off-start)./1200)~=600
          start = find(trig_idx==1,1,'last')
          off = round(mean(find(trig_idx==2)));
        fprintf('S%dM%dB%d: Both triggers, corrected diff = %d\n',isubj,m,iblock,round((off-start)./1200))
        end
      else 
        warning('no triggers')
      end
      
        
      
        
      
    end
  end
end
    
    
    
    
    
    
    
    
    
    
