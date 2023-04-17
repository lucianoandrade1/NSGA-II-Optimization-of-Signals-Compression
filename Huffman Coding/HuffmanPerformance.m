%This Matlab source was implemented to evaluate the performance of NSGA-II
%regarding the optimization of all parameters (wavelets, thresholds,
%scalingFactors and shiftConstants) and generation of the dictionaries used
%in the compression of power quality disturbances signals
%Author: Luciano Andrade

clc;
clear all;
close all;

disp('Loading dictionaries ...');

dicts = load('HuffSignalDictList.mat');
pop = load('HuffSignalPopulation.mat');

[a dictListLenght] = size(dicts.dictList);
[popLenght b] = size(pop.chromosome);

validDicts = cell(1,popLenght);
k=1;
for i=1:2:dictListLenght
    
    for j=1:popLenght

        if (dicts.dictList{1,i}(1)==pop.chromosome(j,1)) & (dicts.dictList{1,i}(2)==pop.chromosome(j,2))
      
            validDicts(k) = dicts.dictList(i);
            validDicts(k+1) = dicts.dictList(i+1);
           
            k = k + 2;
            
        end

    end

end

nroSignals = 100;

pDist = LoadPatterns(nroSignals);

tabWaveletsComp = WaveletsCompTableCreation();

j=1;
for k=1:2:popLenght*2

disp(['compression pop nro: ' num2str(j)]);

%wavelet
wavelet = tabWaveletsComp(abs(round(validDicts{k}(1)))).WaveletComp;

%threshold
threshold = validDicts{k}(2);

scalingFactor = validDicts{k}(3);

shiftConstant = validDicts{k}(4);

%dictionary
savedDict = validDicts{k+1};

for i=1:nroSignals
   
    [a b] = size(pDist(i).DistCurve);
    
    sSig(i)= b;
    
end


qT=[];
for i=1:nroSignals

    [C, L{i}] = wavedec(single(pDist(i).DistCurve),3,wavelet);
    
    [a b] = size(C);
    
    qTSize(i) = b;
    
    quantBool{i} = (abs(C)>threshold);
    
    quantC{i} = C.*quantBool{i};
    
    qT = [qT quantC{i}];
    
end

    qTInt = round(qT*scalingFactor);

    qTIntMax = round(3.1783 * scalingFactor) + 1;
    qTIntMin = round(-3.1705 * scalingFactor) - 1;

    dictArray = [qTIntMin:qTIntMax];

    [a b] = size(dictArray);

    qTIntArray = [dictArray qTInt];

    qTIntArrayL = length(qTIntArray);
    
    code = huffmanenco(qTIntArray,savedDict);
    
    dqTIntArray = huffmandeco(code,savedDict)
    
    dqTInt = dqTIntArray(b+1:end);
    
    dqT = [];
    
    dqT = (dqTInt ./ scalingFactor) + sign(dqTInt).*shiftConstant;
    
    nroOfBitsSegSignal = 64*sum(sSig);

    [x y] = size(code);
    
    nroOfBitsComp = 1*y;
    
    compressionRatio = nroOfBitsSegSignal/nroOfBitsComp;
   
for i=1:nroSignals
   
    if i==1 
        
        iniS=1;
        endS=qTSize(i);
        
        signalDecoded{i} = dqT(iniS:endS); 
        
    else
       
        iniS = endS + 1;
        
        endS = endS + qTSize(i);
        
        signalDecoded{i} = dqT(iniS:endS); 
        
    end
    
end

for i=1:nroSignals

    signalRec = waverec(signalDecoded{i},L{i},wavelet);
    
    mseTemp(i) = mse(pDist(i).DistCurve-signalRec);
    
end

    decodedPop(j,1) = validDicts{k}(1);%wavelet
    decodedPop(j,2) = validDicts{k}(2);%threshold

    decodedPop(j,3) = validDicts{k}(3);%scalingFactor
    decodedPop(j,4) = validDicts{k}(4);%shiftConstant

    decodedPop(j,5) = compressionRatio;
    decodedPop(j,6) = mean(mseTemp);
   
    j = j + 1;
   
end

save('HuffSignalDecoded.mat','decodedPop');

figure;
plot(decodedPop(:,5),decodedPop(:,6),'o');


