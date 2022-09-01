function phase = ph(X);
%PH calculates the phase of a complex series without phase-fixing.
%
%  PHASE = PH(X) calculates the phase of X on a [-pi,pi] range.

%  Author:  Bijan Pesaran 04/25/98

phase=atan2(imag(X),real(X));

