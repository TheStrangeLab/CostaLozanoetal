function [b_boot,b_sig,p_boot,p_sig,shuffles,t_sig,t_boot] = regressPermute(iterations,Ydata,Xdata,shuffles)
% REGRESSPERMUTE - Do a permutation test on regression
% coefficients; this function is the regression equivalent of
% bootstrap.m.
% 
% Performs a permutation test of the regression of e.g., oscillatory
% power on various behavioral variables. The data must be in the
% following format:
%
%   Ydata(events,bin1,bin2,bin3)
%   Xdata(events,regrCoef,bin1,bin2,bin3);
% 
% When using power, the bins are the frequencies/electrodes. The
% statistic can be used to test whether any of the regression
% coefficients is statistically indistinguishable from zero, in
% which can it falls in the middle of the bootstrap distribution.
%
%
%   shuffles = [ones(1,fix(iterations/2)) -ones(1,fix(iterations/2))];
%   ind = randperm2(size(Ydata,1),iterations);
%   shuffles = shuffles(ind)';
%
% FUNCTION:
% function [b_boot,b_sig,p_sig,shuffles,t_sig,t_boot] = regressPermute(iterations,Ydata,Xdata,shuffles)
%
% INPUT ARGS:
%   iterations = 1000; % the number of iterations to run
%   Ydata = mpow_rec;  % oscillatory power for every event,
%   frequency and timebin
%   Xdata = [ones(1,size(targlure) targlure listlen];  %
%   predictors; note that for this function to work, the first
%   column should always be ones (which code for the intercept of
%   the regression)
%   shuffles = shuf;   % result of call to randperm2. (if [], computes it as above)
%
% OUTPUT ARGS: 
% b_boot(iterations,numRegress,bin1,bin2,bin3);% regression coefficients for
% each bin for each iteration
%   b_sig(numRegress,bin1,bin2,bin3); %  actual data with same
%   regression
% p_boot - p-values associated with the regression coefficients for
% the bootstrap
% p_sig - p-values associated with the actual regression coefficients
% shuffles; % The shuffles used for the permutation
% t_sig - t-statistics associated with the actual regression coefficients
% t_boot - t-statistics associated with the regression coefficients
% for the bootstrap

% zero out the data
b_sig = single(zeros(size(Xdata,2),size(Xdata,3),size(Xdata,4),size(Xdata,5)));
p_sig = single(zeros(size(Xdata,2),size(Xdata,3),size(Xdata,4),size(Xdata,5)));
t_sig = single(zeros(size(Xdata,2),size(Xdata,3),size(Xdata,4),size(Xdata,5)));
if iterations > 0
  b_boot = single(zeros(iterations,size(Xdata,2),size(Xdata,3),size(Xdata,4),size(Xdata,5)));
  p_boot = single(zeros(iterations,size(Xdata,2),size(Xdata,3),size(Xdata,4),size(Xdata,5)));
  t_boot = single(zeros(iterations,size(Xdata,2),size(Xdata,3),size(Xdata,4),size(Xdata,5)));
else
  b_boot = [];
  p_boot = [];
  t_boot = [];
end

if iterations > 0 & (~exist('shuffles') | size(shuffles,1)==0)
  shuffles = randperm2(iterations,size(Ydata,1));
end

%fprintf('Iterations(%d): ',iterations);
% do the regression for the actual data
for b1 = 1:size(Ydata,2)  
  for b2 = 1:size(Ydata,3)  
    for b3 = 1:size(Ydata,4)  
      if size(Ydata,2)>1
	[b_sig(:,b1,b2,b3),b_SE,XtXi,Rsq,Fval,Coef_stats] = mregress(squeeze(Ydata(:,b1,b2,b3)),squeeze(Xdata(:,:,b1,b2,b3)),0,0);
	p_sig(:,b1,b2,b3) = Coef_stats(:,4);
	t_sig(:,b1,b2,b3) = Coef_stats(:,3);
      else
	[b_sig,b_SE,XtXi,Rsq,Fval,Coef_stats] = mregress(Ydata,Xdata,0,0);
	p_sig = Coef_stats(:,4);	
	t_sig = Coef_stats(:,3);
      end
    end
  end
end

% do the regression bootstrapped
for i = 1:iterations
  %fprintf('%d ',i);
  for b1 = 1:size(Ydata,2)  
    for b2 = 1:size(Ydata,3)  
      for b3 = 1:size(Ydata,4) 
	thisshuffles = shuffles(i,:);
	%thisshuffles = randperm(size(Ydata,1));
	if size(Ydata,2)>1
	  [b_boot(i,:,b1,b2,b3),b_SE,XtXi,Rsq,Fval,Coef_stats]= mregress(squeeze(Ydata(thisshuffles,b1,b2,b3)),squeeze(Xdata(:,:,b1,b2,b3)),0,0);
	  p_boot(i,:,b1,b2.b3) = Coef_stats(:,4);
	  t_boot(i,:,b1,b2.b3) = Coef_stats(:,3);
	else
	  [b_boot(i,:),b_SE,XtXi,Rsq,Fval,Coef_stats]= mregress(Ydata(thisshuffles),Xdata,0,0);
	  p_boot(i,:) = Coef_stats(:,4);	
	  t_boot(i,:) = Coef_stats(:,3);
	end
      end
    end
  end   
end
%fprintf('Done.\n');
