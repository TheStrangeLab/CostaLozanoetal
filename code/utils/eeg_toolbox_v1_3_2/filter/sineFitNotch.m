function filtered=sineFitNotch(data,sampleRate,notchFreq)

error('JJ: this code has been acting strangely....');

if ~isvector(data)
  error('this function only takes in vectors');
end
data=data(:); %make it a column vector;

t=((1:length(data))./sampleRate)';
s=sin(notchFreq*t*2*pi);
c=cos(notchFreq*t*2*pi);


beta=[s c]\data;

fitWave=beta(1).*s+beta(2).*c;
filtered=data-fitWave;

