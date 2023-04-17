function [f, nroDict, dictList] = iniVariables(nroDict, dictList, pDist,tabWaveletsComp, nroWavelets,...
                          nroThresholds, nroScalingFactors, nroShiftConst, ...
                          nroTotalPop, nroPop, nroVar, nroFObj, ...
                          maxQThreshold, minQThreshold, ...
                          maxScalingFactor, minScalingFactor, ...
                          maxShiftConstant, minShiftConstant, nroSignals)

thresholds = randperm(nroThresholds)/nroThresholds;
scalingFactors = randperm(nroScalingFactors)/nroScalingFactors;
shiftConstant = randperm(nroShiftConst)/nroShiftConst;
                      
disp('Inializing variables...');

ind=1;
for i=1:nroWavelets
    for j = 1:nroThresholds
        for k = 1:nroScalingFactors
            for n=1:nroShiftConst
               
                initialPop(ind).Wavelet = tabWaveletsComp(i).WaveletComp;
                initialPop(ind).Real = tabWaveletsComp(i).Real;
                initialPop(ind).Threshold = minQThreshold + (maxQThreshold - minQThreshold)*thresholds(j);
                initialPop(ind).ScalingFactor = minScalingFactor + (maxScalingFactor - minScalingFactor)*scalingFactors(k);
                initialPop(ind).ShiftConst = minShiftConstant + (maxShiftConstant - minShiftConstant)*shiftConstant(n);

                ind = ind + 1;
            
            end
        
        end
    end
end

hash = randperm(nroTotalPop);
ini = round((nroTotalPop - nroPop) * rand(1));

j=1;
for i=ini:ini+nroPop;
   ind(j) = hash(i);
   j = j + 1;
end

for i=1:nroPop
    for j=1:nroVar
        switch j
            case 1
                f(i,j) = initialPop(ind(i)).Real;
            case 2
                f(i,j) = initialPop(ind(i)).Threshold;
            case 3
                f(i,j) = initialPop(ind(i)).ScalingFactor;
            case 4
                f(i,j) = initialPop(ind(i)).ShiftConst;
        end
    end
    
    
    
[f(i,nroVar + 1: nroVar + nroFObj), firstCount] = EvaluateObjectivesArithmetic(pDist,f(i,:), tabWaveletsComp, ...
                                                   nroFObj, nroVar, nroSignals);

%[f(i,nroVar + 1: nroVar + nroFObj), firstCount] = EvaluateObjectivesRuffman(pDist,f(i,:), tabWaveletsComp, ...
%                                                   nroFObj, nroVar, nroSignals);

count = firstCount;

wavelet = f(i,1);
threshold = f(i,2);
scalingFactor = f(i,3);
shiftConstant = f(i,4); 

compRatio = f(i,nroVar + 1);
distortion = f(i,nroVar + nroFObj);

dictList{nroDict} = [wavelet threshold scalingFactor shiftConstant compRatio distortion];
nroDict = nroDict + 1;

dictList{nroDict} = count;
nroDict = nroDict + 1;

end

end

