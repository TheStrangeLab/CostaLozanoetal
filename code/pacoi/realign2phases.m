function [angbin1a, angbin2a2, flag, angbin2a]= realign2phases(angbin1,angbin2,x_dbin)
% angbin1/2 are trials x phase bins. The realignment is perform using the
% phase bin pointed by the vector strength of anlge1. angbin2 moves same
% amount of bins. In a second step, the phase bins of angbin2 are realigned
% towards the positive half of the x-axis
% 
% xx = linspace(-pi,pi,61);
% x_dbin = mean([xx(1:end-1);xx(2:end)],1);% the size(x_dbin)/2 must be even
% phi = linspace(-pi,pi,8);
% angbin1 = cos(x_dbin);
% 
% for s=1:8
%   angbin2 = cos(x_dbin+phi(s));
%   figure(1);  subplot(4,2,s);plot(x_dbin,[angbin1' angbin2']);title('original')
%   [ angbin1a, angbin2a2, flag(s)] = realign2phases(angbin1,angbin2,x_dbin);
%   figure(2);subplot(4,2,s);plot(x_dbin,[angbin1a' angbin2a2']);title('after realignment')
% end
%
% 13-Nov-2019 00:28:24 Diego Lozano-Soldevilla Diego Lozano-Soldevilla
%

nl = size(x_dbin,2);
zerobin = (nl/2+1);
ang_p1o = angle(mean(angbin1.*exp(sqrt(-1)*x_dbin),2));
ang_p2o = angle(mean(angbin2.*exp(sqrt(-1)*x_dbin),2));

[val,pos]=min(abs(ang_p1o-x_dbin));
nshift1 = zerobin-pos;
angbin1a = circshift(angbin1,nshift1);
angbin2a = circshift(angbin2,nshift1);

% now let's mirror right-shift angle 2 relative to y-axis (in case the
% angle of vector-strength points between [-pi 0])
ang_p2c = angle(mean(angbin2a.*exp(sqrt(-1)*x_dbin)));% current phase angle after realignment

% mirror-flip relative to y-axis
if (ang_p2o <= 0  & ang_p2o >= -pi) & (ang_p2c <= 0 & ang_p2c >= -pi) |...
   (ang_p2o <= pi & ang_p2o >= 0)   & (ang_p2c <= 0 & ang_p2c >= -pi)
  % find the closest position to the
  [val,pos]=min(abs(ang_p2c-x_dbin));
  sym_x = [-nl/2:-1 1:nl/2];
  pos_sym = find(sym_x==-sym_x(pos));
  nshift2 = pos_sym-pos;
  angbin2a2 = circshift(angbin2a,nshift2);
  flag = -1;
else
  angbin2a2 = angbin2a;
  flag = 0;
end

% clf;
% subplot(121);plot(x_dbin,angbin1, x_dbin,angbin2)
% subplot(122);plot(x_dbin,angbin1a,x_dbin,angbin2a2)
