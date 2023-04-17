function [f, count] = EvaluateObjectivesArithmetic(signalSeg, varDec, tabWaveletsComp, nroFObj, nroVar, nroSignals)

wavelet = tabWaveletsComp(round(varDec(1))).WaveletComp;
threshold = varDec(2);
scalingFactor = varDec(3);
shiftConstant = varDec(4);

disp(['Compression: ' wavelet ' ' num2str(threshold)]);

for i=1:nroSignals
   
    [a b] = size(signalSeg(i).DistCurve);
    
    sSig(i)= b;
    
end

qT=[];
for i=1:nroSignals

    [C, L{i}] = wavedec(signalSeg(i).DistCurve,3,wavelet);
    
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
    
    qTIntArray = qTIntArray + abs(qTIntMin) + 1;
        
    qTIntArrayL = length(qTIntArray);
    
    Z = unique (qTIntArray);

    count=histc(qTIntArray,Z);
    
    code = arithenco(qTIntArray,count);

    dqTIntArray = arithdeco(code,count,qTIntArrayL);

    dqTIntArray = dqTIntArray - abs(qTIntMin) - 1;
    
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
    
    mseTemp(i) = mse(signalSeg(i).DistCurve-signalRec);
    
end

wavelet = tabWaveletsComp(round(varDec(1))).WaveletComp;
threshold = varDec(2);
scalingFactor = varDec(3);
shiftConstant = varDec(4);

    f(1) = - compressionRatio;
    f(2) = mean(mseTemp);

end
 
