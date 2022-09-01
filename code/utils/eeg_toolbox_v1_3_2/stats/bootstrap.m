function [p_boot,p_sig,shuffles] = bootstrap(iterations, data1, data2, shuffles, doWilcox, doSig)
%BOOTSTRAP - Perform a paired or non-paired permutation test.
% 
% Performs a permutation test comparing two sets of data, such as
% binned ERPs of recalled and not recalled words or different
% frequencies of a mean power analysis.  The data must be in the
% following format:
%
%   data1(events,bin1,bin2,bin3)
% 
% When using power, the bins are the frequencies/electrodes.  For the
% non-paired test, the number of instances/events can be different for
% data1 and data2, but the size of bin1, bin2, and bin3 must always be
% the same.
%
% To perform a paired sample test, simply take the difference of your
% pairwise samples and pass them into data1, while passing in an empty
% matrix [] to data2.
%
% The statistic tests the probability of rejecting the hypothesis that
% data1 > data2 for non-paired data or that data1 > 0 for paired data.
% Therefore, 1-p is the probability of rejecting the hypothesis that
% data2 > data1 or that data1 < 0.
%
% The proper bootstrap method is to use the same shuffle across
% channels and frequencies, so you must pass in the shuffles to be
% used for the bootstrap randomization.  You would create the shuffles
% with the randperm2 command to be the sum of the number of events in
% data1 and data2:
%
% Non-paired data:
%   shuffles = randperm2(iterations,size(data1,1)+size(data2,1));
%
% Paired data:
%   shuffles = [ones(1,fix(iterations/2)) -ones(1,fix(iterations/2))];
%   ind = randperm2(size(data1,1),iterations);
%   shuffles = shuffles(ind)';
%
% FUNCTION:
%   [p_boot,p_sig,shuffles] = bootstrap(iterations,data1,data2,shuffles,doWilcox,doSig)
%
% INPUT ARGS:
%   iterations = 1000; % the number of iterations to run
%   data1 = mpow_rec;  % data from first condition
%   data2 = mpow_not;  % data from second condition
%   shuffles = shuf;   % result of call to randperm2. (if [], computes it as above)
%   doWilcox = 1;         % pick 1 to do non-parametric test
%   instead of ttest
%   doSig = 1          % whether to return the significance of the
%   actual data as well
%
% OUTPUT ARGS:
%   p_boot(iterations,bin1,bin2,bin3); % p values for each bin for each
%                            % iteration
%   p_sig(bin1,bin2,bin3);  % actual data with same sig. test
%   shuffles- The shuffles used for the permutation
%

% 4/5/06 - PBS - Reverted to use shuffles again.
% 6/11/05 - PBS - Fixed bug in replicating shuffles for paired samples. 
% 4/1/05 - PBS - Updated to handle paired samples.



if ~exist('doWilcox','var')
  doWilcox = 1;
end
if ~exist('doSig','var')
  doSig = 1;
end

% see if we are dealing with paired data
isPaired = 0;
if isempty(data2)
  isPaired = 1;
end

if iterations > 0 & (~exist('shuffles') | size(shuffles,1)==0)
  if isPaired
    % paired data shuffles is just a random switch of sign
    shuffles = [ones(1,fix(iterations/2)) -ones(1,fix(iterations/2))];
    ind = randperm2(size(data1,1),iterations);
    shuffles = shuffles(ind)';
  else
    % not paired shuffles
    shuffles = randperm2(iterations,size(data1,1)+size(data2,1));
  end
end

% zero out the data
p_sig = single(zeros(size(data1,2),size(data1,3),size(data1,4)));
if iterations > 0
  p_boot = single(zeros(iterations,size(data1,2),size(data1,3),size(data1,4)));
else
  p_boot = [];
end

fprintf('Iterations(%d): ',iterations);


% set up the data_all
if isPaired
  % do the actual data
  if doSig
    if doWilcox == 1
      % nonparametric test
      for b1 = 1:size(data1,2)  
        for b2 = 1:size(data1,3)  
          for b3 = 1:size(data1,4)  
            % do a non-parametric
            [p_sig(b1,b2,b3)] = signrank_tail(data1(:,b1,b2,b3),[],'tail','right');
          end
        end
      end
    else
      % parametric test
      for b1 = 1:size(data1,2)  
        for b2 = 1:size(data1,3)  
          for b3 = 1:size(data1,4)  
            % do a parametric
            [h,p_sig(b1,b2,b3)] = ttest(data1(:,b1,b2,b3),0,.05,1);
          end
        end
      end      
    end
  end
  
  % do the paired data tests
  for i = 1:iterations
    fprintf('%d ',i);
    
    % randomly flip the signs of the data, use same flip for all bins
    temp_shuffle = shuffles(i,:)';
    %temp_shuffle = sign(randn(size(data1,1),1));
    temp_shuffle = temp_shuffle(:,ones(1,size(data1,2),1,1),ones(1,1,size(data1,3),1),ones(1,1,1,size(data1,4)));
    data1_rand = data1.*temp_shuffle;
    
    % see if nonparametric
    if doWilcox == 1
      % nonparametric test
      for b1 = 1:size(data1,2)  
        for b2 = 1:size(data1,3)  
          for b3 = 1:size(data1,4)  
            % do a non-parametric
            [p_boot(i,b1,b2,b3)] = signrank_tail(data1_rand(:,b1,b2,b3),[],'tail','right');
          end
        end
      end
      
    else
      % parametric test
      for b1 = 1:size(data1,2)  
        for b2 = 1:size(data1,3)  
          for b3 = 1:size(data1,4)  
            % do a parametric
            [h,p_boot(i,b1,b2,b3)] = ttest(data1_rand(:,b1,b2,b3), ...
					   0,.05,1);
	  end
        end
      end      
    end
  end
else
  % do nonpaired data tests
  
  % first do the actual data
  if doSig
    if doWilcox == 1
      for b1 = 1:size(data1,2)  
        for b2 = 1:size(data1,3)  
          for b3 = 1:size(data1,4)  
            % do a non-parametric
            [p_sig(b1,b2,b3)] = ranksum_ci(data1(:,b1,b2,b3),data2(:,b1,b2,b3),.05,1);
          end
        end
      end
    else
      % parametric test
      for b1 = 1:size(data1,2)  
        for b2 = 1:size(data1,3)      
          for b3 = 1:size(data1,4)      
            % do a ttest
            [h,p_sig(b1,b2,b3)] = ttest2(data1(:,b1,b2,b3),data2(:,b1,b2,b3),.05,1);
          end
        end
      end
    end
  end
  
  for i = 1:iterations
    fprintf('%d ',i);
    
    % randomize the data
    data_all = [data1 ; data2];
    data_all = data_all(shuffles(i,:),:,:,:);
    %data_all = data_all(randperm(size(data1,1)+size(data2,1)),:,:,:);
    data1_rand = data_all(1:size(data1,1),:,:,:);
    data2_rand = data_all(size(data1,1)+1:end,:,:,:);

    
    if doWilcox == 1
      % do a non-parametric
      % do the statistical tests
      
      for b1 = 1:size(data1,2)  
        for b2 = 1:size(data1,3)  
          for b3 = 1:size(data1,4)  
            % do a non-parametric
            [p_boot(i,b1,b2,b3),h,stats,Ws] = ranksum_ci(double(data1_rand(:,b1,b2,b3)), ...
                                                         double(data2_rand(:,b1,b2,b3)),.05,1);
          end
        end
      end
      
    else
      for b1 = 1:size(data1,2)  
        for b2 = 1:size(data1,3)      
          for b3 = 1:size(data1,4)      
            % do a ttest
            [h,p_boot(i,b1,b2,b3),ci,stats] = ttest2(double(data1_rand(:,b1,b2,b3)), ...
                                                     double(data2_rand(:,b1,b2,b3)),.05,1);
          end
        end
      end
    end
  end
end

fprintf('Done.\n');

%p_boot=squeeze(p_boot);
