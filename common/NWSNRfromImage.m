function [SNR,sigmask,noisemask] = NWSNRfromImage(img)
% assumes that signal mean is ALWAYS higher than noise mean


h = figure; imagesc(img)

% h1 = msgbox('Choose SIGNAL region','ROIs');
hS = imellipse(gca);
wait(hS);
sigmask = createMask(hS);
setColor(hS,'r')
% h2 = msgbox('Choose NOISE region','ROIs');
hN = imellipse(gca);
wait(hN);
noisemask = createMask(hN);
setColor(hN,'r')

tempsig = sigmask;
tempnoise = noisemask;

if mean2(img(tempsig)) < mean2(img(tempnoise)) % switch signal and noise
    sigmask = tempnoise;
    noisemask = tempsig;
end  

SNR = mean2(img(sigmask))/std2(img(noisemask));

close(h);

