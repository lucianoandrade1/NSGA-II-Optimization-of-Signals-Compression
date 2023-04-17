%This souce code was implemented to help the implementation of NSGA-II 
%Power Quality disturbance signals compression 
%Author: Luciano Andrade

clc;
clear all;
close all;

nroSignals = 100;
nroWavelets = 37;

tabWaveletsComp = WaveletsCompTableCreation();

pDist = LoadPatterns(nroSignals);

for k=1:nroWavelets

    wavelet = tabWaveletsComp(k).WaveletComp;

    allSignals=[];
    for i=1:nroSignals

        [C, L] = wavedec(pDist(i).DistCurve,3,wavelet);

        allSignals = [allSignals C];

    end
    
    maxQt(k) = max(allSignals);
    
    minQt(k) = min(allSignals);

end

[qTMax, Imax] = max(maxQt);
[qTMin, Imin] = min(minQt);

maxWavelet = tabWaveletsComp(Imax).WaveletComp;
minWavelet = tabWaveletsComp(Imin).WaveletComp;

disp(['Max signals value: ' num2str(qTMax) ' to ' maxWavelet]);
disp(['Min signals value: ' num2str(qTMin) ' to ' minWavelet]);

