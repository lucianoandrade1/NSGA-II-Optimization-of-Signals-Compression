function [f, dict] = EvaluateObjectivesRuffman(signalSeg, varDec, tabWaveletsComp, nroFObj, nroVar, nroSignals)

wavelet = tabWaveletsComp(round(varDec(1))).WaveletComp;
limiar = varDec(2);
scalingFactor = varDec(3);
shiftConstant = varDec(4);

disp(['Compression: ' wavelet ' ' num2str(limiar)]);

for i=1:nroSignals
   
    [a b] = size(signalSeg(i).DistCurve);
    
    sSig(i)= b;
    
end

qT=[];
for i=1:nroSignals

    [C, L{i}] = wavedec(signalSeg(i).DistCurve,3,wavelet);
    
    [a b] = size(C);
    
    qTSize(i) = b;
    
    quantBool{i} = (abs(C)>limiar);
    
    quantC{i} = C.*quantBool{i};
    
    qT = [qT quantC{i}];
    
end

    qTInt = round(qT*scalingFactor);

    qTIntMax = round(3.1783 * scalingFactor) + 1; %max value for Sag = 2.8826 - checkBorders
    qTIntMin = round(-3.1705 * scalingFactor) - 1; %min value for Sag = -2.8826 - checkBorders

    dictArray = [qTIntMin:qTIntMax];

    [a b] = size(dictArray);

    qTIntArray = [dictArray qTInt];
    
    qTIntArrayL = length(qTIntArray);
    
    %min(qTIntArray)
    
    Z = unique (qTIntArray);

    countElY=histc(qTIntArray,Z); %# get the count of elements

    p = countElY/numel(qTIntArray); %getting the probability distribution

    [dict,avglen] = huffmandict(Z,p); % Create dictionary.

    %disp('dictnary done');
    
    comp = huffmanenco(qTIntArray,dict); % Encode the data.
    
    %disp('encoded');
    %-----
    % MANTER ESSE CODIGO
    %dqTIntArray = huffmandeco(comp,dict);
    
    %disp('decoded');
    
    dqTInt = qTIntArray(b+1:end);%dqTInt = dqTIntArray(b+1:end);
    
    dqT = [];
    
    dqT = (dqTInt ./ scalingFactor) + sign(dqTInt).*shiftConstant;
    
    nroOfBitsSegSignal = 64*sum(sSig);

    [x y] = size(comp);
    
%     disturbanceClass = 2^3;%8 classes de disturbios
%     
%     LiSize = nroSignals * 64 * 5; %L{i} size;
    
    nroOfBitsComp = 1*y;% + disturbanceClass + LiSize;
    
    compressionRatio = nroOfBitsSegSignal/nroOfBitsComp;
    
    %save('fileDeb','qT','dqT');
    
    
%MANTER ESSE CODIGO    
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
    
    %signalRec = waverec(quantC{i},L{i},wavelet);
    
    mseTemp(i) = mse(signalSeg(i).DistCurve-signalRec);
    
end

wavelet = tabWaveletsComp(round(varDec(1))).WaveletComp;
limiar = varDec(2);
scalingFactor = varDec(3);
shiftConstant = varDec(4);

    f(1) = - compressionRatio;
    f(2) = mean(mseTemp);

%     disp([wavelet ' - ' ...
%         num2str(limiar) ' - ' ...
%         num2str(scalingFactor) ' - ' ...
%         num2str(shiftConstant) ' - ' ...
%         num2str(compressionRatio) ' - ' ...
%         num2str(mean(mseTemp))]);
    
end
 


