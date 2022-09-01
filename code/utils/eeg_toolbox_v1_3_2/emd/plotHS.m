function ifh=plotHS(d,sampleRate)
debug=true;

m=jemd(d);
[f,amp]=calcIF(m,sampleRate);

f(f<0)=nan; %don't plot negative freqs
t=(1:size(d,2))./sampleRate;
subplot(3,1,1);
plot(t,d);

subplot(3,1,2);
%f2=f;f2(amp<.5e-4)=nan;
plot(t,f);
set(gca,'yscale','log','ytick',2.^[1:8],'ylim',[1 256]);
subplot(3,1,3);
ifh=doTFPlot(f,amp,sampleRate);

%colorbar south
%keyboard

function ifh=doTFPlot(freqs,amp,sr);
times=(1:size(freqs,2))./sr;

numYBins=160;
%numYBins=20;
freqBins=logspace(log10(1),log10(256),numYBins);

ifh=zeros(numYBins,length(times));
for t=1:length(times)
  [n,freqNumbers]=histc(freqs(:,t),freqBins);
  for modeNum=1:length(freqNumbers)
    f=freqNumbers(modeNum);
    if f~=0%&amp(f,t)>100
%      ifh(f,t)=ifh(f,t)+1;
      ifh(f,t)=ifh(f,t)+amp(modeNum,t);
    end
  end
end

pcolor(times,freqBins,ifh);
set(gca,'yscale','log','ytick',2.^[1:8]);
shading flat
  
  


