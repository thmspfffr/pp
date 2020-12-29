function [h,p] = tp_singlethreshpermtest(par1,par2,alpha,nperm)


if ndims(par1) == 3
  
  allpar(:,:,:,1) = par1;
  allpar(:,:,:,2) = par2;
  
  perm_idx = randi(2,[size(par1,3),nperm]);

  for iperm = 1 : nperm
    
    % within subjects permutation test
    idx1 = perm_idx(:,iperm);
    idx2 = 3-idx1;
    
    for i = 1 : length(idx1)
      permdat(:,:,i,1) = allpar(:,:,i,idx1(i));
      permdat(:,:,i,2) = allpar(:,:,i,idx2(i));
    end
    
    [~,~,~,s] = ttest(permdat(:,:,:,2),permdat(:,:,:,1),'dim',3,'tail','both','alpha',alpha);
    max_tval(iperm) = max(abs(s.tstat(:)));
    
  end
  
  [~,~,~,s] = ttest(allpar(:,:,:,2),allpar(:,:,:,1),'dim',3,'tail','both','alpha',alpha);
    
  ss = sum(s.tstat >permute(repmat(max_tval',[1 size(s.tstat,1) size(s.tstat,2)]),[2 3 1]),3);
  
  p = 1-(ss./nperm);
  h = p < 0.05;

elseif ndims(par1) == 2
  
  allpar(:,:,1) = par1;
  allpar(:,:,2) = par2;
  
  perm_idx = randi(2,[size(par1,2),nperm]);

  for iperm = 1 : nperm
    
    % within subjects permutation test
    idx1 = perm_idx(:,iperm);
    idx2 = 3-idx1;
    
    for i = 1 : length(idx1)
      permdat(:,i,1) = allpar(:,i,idx1(i));
      permdat(:,i,2) = allpar(:,i,idx2(i));
    end
    
    [~,~,~,s] = ttest(permdat(:,:,2),permdat(:,:,1),'dim',2,'tail','both','alpha',alpha);
    max_tval(iperm) = max(abs(s.tstat(:)));
    
  end
  
  [~,~,~,s] = ttest(allpar(:,:,2),allpar(:,:,1),'dim',2,'tail','both','alpha',alpha);
    
  ss = sum(s.tstat>repmat(max_tval',[1 size(s.tstat,1) size(s.tstat,2)])',2);
  
  p = 1-(ss./nperm);
  h = p < 0.05;
  
  
else
  
  error('Input needs to be 2- or 3-dimensional array!')
  
end