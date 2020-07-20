%% FIGURE OUT GENERAL DATA ORGA STUFF
% pp_organize_data_resttask

% SUBJLIST  = [4 5 6 7 8 9 10 11 12 13 15 16 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34];
addpath ~/Documents/MATLAB/fieldtrip-20160919/
ft_defaults
% for isubj = SUBJLIST
%   for m = 1:3
%     
%     length(dir(sprintf('/home/tpfeffer/pconn/rawdata/edf/p4/s1/edf_samples/*asc')))
% 
% 
%     fprintf('Processing subject %d, session%d ...\n',isubj,m)
%     
%     d = dir(sprintf('/home/tpfeffer/pconn/rawdata/edf/p%d/s%d/edf_samples/*%d_b*',isubj,m,m));
%     d_evts  = dir(sprintf('/home/tpfeffer/pconn/rawdata/edf/p%d/s%d/edf_events/*%d_b*',isubj,m,m));
%     
%     
%         
%   end
% end

%%
%     if length(d)~=length(d_evts)
%       error('!')
%     end
%     
% clear cond block
% for isubj = SUBJLIST
%   for m = 1:3
%     isubj
%    
% %     for iblock = 1 : length(d)
%       
% %       block = str2num(d(iblock).name(regexp(d(iblock).name,'-','once')+5));
%       
%       d_meg = dir(sprintf('/home/tpfeffer/pconn/rawdata/meg/p%d/s%d/*ds',isubj,m));
%       iblock = 0
%       for imeg = 1 : length(d_meg)
%         imeg
%         f = dir([d_meg(imeg).folder '/' d_meg(imeg).name '/' '*res4']);
%         a = ft_read_event([f.folder '/' f.name]);
% %         if size(cell2mat({a.value}),2)<10
% %           cond(imeg,isubj,m) = 0;
% %           block(imeg,isubj,m) = str2num(d_meg(1).name(end-3));
% %         elseif sum(cell2mat({a.value})>70)>5
%           if sum(cell2mat({a.value})>70)>5
%           cond(imeg,isubj,m) = 2;
%            iblock = iblock + 1;
%           block(imeg,isubj,m) = str2num(d_meg(1).name(end-3));
%           trig = cell2mat({a(strcmp({a.type},'UPPT002')).sample});
%           save(sprintf(['~/pconn_bttn/proc/' 'pconn_bttn_trig_samples_s%d_m%d_b%d_v%d.mat'],isubj,m,iblock,1),'trig');
% %         elseif size(cell2mat({a.value}),2) > 10 && sum(cell2mat({a.value})>70)<=5
% %           cond(imeg,isubj,m) = 1;
% %           block(imeg,isubj,m) = str2num(d_meg(1).name(end-3));
%           
%         end
%       end
%     
%   end
%   
% end

% from notes in lab book: s10 m3
% cond(:,10,3) = [0 2 1 0 2 1];
% cond(:,33,3) = [1 2 0 1 2 0];
% cond(:,30,2) = [1 2 0 1 2 0];

% save('~/pp/proc/blocks_taskrest.mat','cond')

%% OLD SCRIPT FOR EXTRACTING DATA
indir1  = '/home/tpfeffer/pconn/rawdata/meg/';
clear cond
SUBJLIST  = [4 5 6 7 8 9 10 11 12 13 15 16 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34];

for m = 1:3
for isubj = SUBJLIST

%   isubj
  indir = sprintf([indir1 'p%d/s%d/'],isubj,m);
  fprintf('\n')
%   try
    
%     if ~exist([outdir sprintf('pconn_preproc_s%d_m%d_v%d_processing.txt',isubj,m,v)]) 
% %         && ~exist([outdir sprintf('pconn_preproc_data_s%d_m%d_v%d.mat',isubj,m,v)])
%       
%       system(['touch ' outdir sprintf('pconn_preproc_s%d_m%d_v%d_processing.txt',isubj,m,v)]);
%       disp(sprintf('Processing subject s%d m%d',isubj,m));
%       
%       % no proc file, but output: create proc file
%     elseif ~exist([outdir sprintf('pconn_preproc_s%d_m%d_v%d_processing.txt',isubj,m,v)]) ...
%         && exist([outdir sprintf('pconn_preproc_data_s%d_m%d_v%d.mat',isubj,m,v)])
%       
%       system(['touch ' outdir sprintf('pconn_preproc_s%d_m%d_v%d_processing.txt',isubj,m,v)]);
%       continue
%       
%     elseif exist([outdir sprintf('pconn_preproc_s%d_m%d_v%d_processing.txt',isubj,m,v)]) ...
%         && ~exist([outdir sprintf('pconn_preproc_data_s%d_m%d_v%d.mat',isubj,m,v)])
%       
%       disp(sprintf('Processing subject s%d m%d',isubj,m));
      
%     else
%     continue
%     end
%     
    cont_dir = dir(indir);
    
    for i = 1 : length(cont_dir)
      tmp(i) = ~isempty(regexp(cont_dir(i).name,'(.ds)','once'));
    end
    
    ind = find(tmp>0);
    
    ibl(1:3) = 0;
    
    % no button presses registered for subject 10, m = 3
    if isubj == 10 && m == 3
      ind = [ind(1) ind(4)];
    end
    
    for idir = ind
%       idir
      tmp_dataset      = [cont_dir(idir).name];

        cfgs.datafile    = [indir '/' tmp_dataset '/' cont_dir(idir).name(1:end-2) 'meg4'];
        cfgs.headerfile  = [indir '/' tmp_dataset '/' cont_dir(idir).name(1:end-2) 'res4'];
      
      cfg = [];
      % cfg.continuous = 'yes';
      cfg.dataset    = [indir cont_dir(idir).name];
      
      if isubj == 22 && m == 3
        a=ft_read_event(cfgs.headerfile);
      else
        a=ft_read_event(cfg.dataset);
      end
%       
%       if sum(cell2mat({a.value})>100)<9 && sum(cell2mat({a.value})==50)==0
%           fprintf('S%d: Rest! Presses: %d\n',isubj,sum(cell2mat({a.value})>100));
%           cond(idir-2,isubj,m) = 0;
%           press(ind-2,isubj,m) = sum(cell2mat({a.value})>100);
%           
% %           continue
%       elseif sum(cell2mat({a.value})>100)>8 && sum(cell2mat({a.value})==50)>10
%           fprintf('S%d: Task - button press! Presses: %d\n',isubj,sum(cell2mat({a.value})>100));
%           cond(idir-2,isubj,m) = 2;
%           press(ind-2,isubj,m) = sum(cell2mat({a.value})>100);
% %           continue
%       elseif sum(cell2mat({a.value})>100) < 9  && sum(cell2mat({a.value})==50)>10
%           fprintf('S%d: Task - counting! Presses: %d\n',isubj,sum(cell2mat({a.value})>100));
%           cond(idir-2,isubj,m) = 1;
%           press(ind-2,isubj,m) = sum(cell2mat({a.value})>100);
% %           continue
%       else
%         warning('NOTHING FOUND')
%       end
      
       if ~(isubj ==10 && m==3)
        if size(cell2mat({a.value}),2)<10
            fprintf('S%d: Rest! Presses: %d\n',isubj,sum(cell2mat({a.value})>100));
          cond(idir-2,isubj,m) = 0;
          press(ind-2,isubj,m) = sum(cell2mat({a.value})>100);
        elseif sum(cell2mat({a.value})>50)>7
            fprintf('S%d: Task - button press! Presses: %d\n',isubj,sum(cell2mat({a.value})>100));
          cond(idir-2,isubj,m) = 2;
          press(ind-2,isubj,m) = sum(cell2mat({a.value})>100);
        elseif size(cell2mat({a.value}),2) > 10 && sum(cell2mat({a.value})>50)<7
            fprintf('S%d: Task - counting! Presses: %d\n',isubj,sum(cell2mat({a.value})>100));
          cond(idir-2,isubj,m) = 1;
          press(ind-2,isubj,m) = sum(cell2mat({a.value})>100);
  %           continue
        end
      end


    end

end
end

%  if ~(isubj ==10 && m==3)
%         if size(cell2mat({a.value}),2)<10
%             warning('Rest!');
%             cond = 1;
%             continue
%         elseif sum(cell2mat({a.value})>50)>7
%             warning('Task - button press!');
%             cond = 2;
%             continue
%         elseif size(cell2mat({a.value}),2) > 10 && sum(cell2mat({a.value})>50)<7
%             warning('Task - counting!');
%             cond = 3;
%   %           continue
%         end
%       end