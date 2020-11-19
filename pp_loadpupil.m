function allpup = pp_loadpupil(para)
% pp_loadpupil Loads mean pupil diameter within defined period of time
% allpup = pp_loadpupil(para)
% -----------------
% para.time = [t1 t2] - time interval from onset of recording
% para.subjlist - subjects to analyze (default: all)
% -----------------

clear ts_cnt ts_bttn allpup_cnt allpup_bttn
if ~isfield(para,'time')
  para.time = [3000 7000];
end

if ~isfield(para,'subjlist')
  SUBJLIST  = [4 5 6 7 8 9 10 11 12 13 15 16 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34];
else
  SUBJLIST = para.subjlist;
end

ord = pconn_randomization;

fprintf('Loading resting state ...\n');

for isubj = SUBJLIST
  for m = 1 :3
    im = find(ord(isubj,:)==m);    
    for iblock = 1 : 2
      try
        load(sprintf('~/pp/proc/pup/pp_pupil_diameter_cleaned_s%d_m%d_b%d.mat',isubj,im,iblock))
        allpup(isubj,m,iblock,1) = mean(pupil(para.time(1):para.time(2),4));
        clear pupil
      catch me
        allpup(isubj,m,iblock,1) = nan;
      end
    end
  end
end

fprintf('Loading counting ...\n');

for isubj = SUBJLIST
  isubj
  for m = 1 :3
    im = find(ord(isubj,:)==m);    
    for iblock = 1 : 2
      try
        load(sprintf('~/pp/proc/pup/pp_pupil_diameter_cleaned_counting_s%d_m%d_b%d.mat',isubj,im,iblock))
        allpup(isubj,m,iblock,2) = mean(pupil(para.time(1):para.time(2),4));
        clear pupil
      catch me
        allpup(isubj,m,iblock,2) = nan;
      end
    end
  end
end

% fprintf('Loading pressing ...\n');
% 
% 
% for isubj = SUBJLIST
%   isubj
%   for m = 1 :3
%     im = find(ord(isubj,:)==m);    
%     for iblock = 1 : 2
%       try
%         load(sprintf('~/pp/proc/pup/pp_pupil_diameter_cleaned_button_s%d_m%d_b%d.mat',isubj,im,iblock))
%         allpup(isubj,m,iblock,3) = mean(pupil(para.time(1):para.time(2),4));
% 
%         clear pupil
%       catch me
%         allpup(isubj,m,iblock,3) = nan;
%       end
%     end
%   end
% end
% 
allpup = allpup(SUBJLIST,:,:,:);
