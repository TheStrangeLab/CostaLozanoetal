function tab = F1LDF1(data,group,subject,time)
% function tab = F1LDF1(data,group,subject,time)
%
% nonparametric ANOVA for a two-way repeated measures design: two
% groups measured at two points in time. Note that group and time can
% only have 2 different IDs. Following Brunner et al (XXXX)
%
% INPUT Args:
%  data = [1 2.4 1.1 0.8...]
%  group = [1 1 1 1 1 2 2 2 2 ...] -- group ID for every datapoint
%  subject = [1 2 3 4 5 6 7 8 9...] - subject ID
%  time = [1 1 1 1 1 1 1 1 .......... 2 2 2 2 2 ] - timepoint
% 
% OUTPUT Args:
% tab = [rankGroup pvalGroupN pvalGroupT; rankTime pvalTimeN
% pvalTimeT; rankInteract pvalInteractN pvalInteractT]; where the
% postfix N or T indicates whether based on normal statistics or
% t-statistics (small-group statistic)


% 1) rank all data (all ranked together)
rdata = tiedrank(data);
% 2) sum ranks over time 1,2 and group 1,2
rmean = zeros(2,2);
timeIDs = unique(time);
groupIDs = unique(group);
% find the number of participants in each group
gr1 = find(group==groupIDs(1));
Ngrp1 = length(unique(subject(gr1)));
gr2 = find(group==groupIDs(2));
Ngrp2 = length(unique(subject(gr2)));
n = [Ngrp1 Ngrp2];
for t = 1:2
  for g = 1:2
    % sum the ranks within the group
    rmean(g,t) = sum(rdata(intersect(find(group==groupIDs(g)),find(time==timeIDs(t)))));
  end
end
rmean(1,:) = rmean(1,:)./Ngrp1;
rmean(2,:) = rmean(2,:)./Ngrp2;
% 3) group effect
sigma2 = zeros(1,2);
for g  = 1:2
  t1 = rdata(intersect(find(group==groupIDs(g)),find(time==timeIDs(1))));
  t2 = rdata(intersect(find(group==groupIDs(g)),find(time==timeIDs(2))));
  sigma2(g) = (1/(n(g)-1))*sum((t1 + t2 - repmat(rmean(g,1),n(g),1) -repmat(rmean(g,2),n(g),1)).^2);
end
UnA = (rmean(1,1)+rmean(1,2)-rmean(2,1)-rmean(2,2))/sqrt(sum(sigma2./n));
% 4) time effect
tau2 = zeros(1,2);
for g  = 1:2
  g1 = rdata(intersect(find(group==groupIDs(g)),find(time==timeIDs(1))));
  g2 = rdata(intersect(find(group==groupIDs(g)),find(time==timeIDs(2))));
  tau2(g) = (1/(n(g)-1))*sum((g1 - g2 - repmat(rmean(g,1),n(g),1) + repmat(rmean(g,2),n(g),1)).^2);
end
UnT = (rmean(1,1) - rmean(1,2) + rmean(2,1) - rmean(2,2))/sqrt(sum(tau2./n));
% 5) interaction
UnAT = (rmean(1,1)-rmean(1,2) - rmean(2,1) + rmean(2,2))/sqrt(sum(tau2./n));
% 6) compute p-value, either using the normal distribution or the t
% distribution (give both)
pAN = 2*(1-normcdf(abs(UnA)));
nA = (sum(sigma2./n)^2)/(sum((sigma2./n).^2./(n-1))); % CHECK this !!
pAT = 2*(1-tcdf(abs(UnA),nA));
nT = (sum(tau2./n)^2)/(sum((tau2./n).^2./(n-1)));
pTN = 2*(1-normcdf(abs(UnT)));
pTT = 2*(1-tcdf(abs(UnT),nT));
pATN = 2*(1-normcdf(abs(UnAT)));
pATT = 2*(1-tcdf(abs(UnAT),nT));

tab = [UnA pAN pAT; UnT pTN pTT; UnAT pATN pATT];
