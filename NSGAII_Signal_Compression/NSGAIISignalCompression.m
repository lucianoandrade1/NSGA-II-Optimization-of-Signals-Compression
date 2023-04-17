%Matlab 2014 source code for Power Quality Signals Disturbances Compression
%Author: Luciano Andrade
%This souce code was adapted from original NSGA-II available in Mathworks 
%website library to Power Quality signals disturbances compressio


clc;
clear all;
close all;

set(0,'RecursionLimit',5100)

nroThresholds = 50;
nroScalingFactors = 64;
nroShiftConst = 20;

nroVar = 4;
nroFObj = 2;

nroG = 25;

maxQLimiar = 0.0700;
minQLimiar = 0.0050;

maxScalingFactor = 2^9;
minScalingFactor = 2^7;

maxShiftConstant = 1/2^8;
minShiftConstant = 0.0;

nroSignals = 25;

pDisturbbances = LoadPatterns(nroSignals);

pDist = pDisturbbances;

tabWaveletsComp = WaveletsCompTableCreation();

[a nroWavelets] = size(tabWaveletsComp);

nroTotalPop = nroWavelets*nroThresholds*nroScalingFactors*nroShiftConst;

nroPop = 100;

nroDict = 1;
dictList = cell(1, nroG*nroPop);

[chromosome, nroDict, dictList] = iniVariables(nroDict, dictList, pDist,tabWaveletsComp, nroWavelets, ...
                                              nroThresholds, nroScalingFactors, nroShiftConst,...
                                              nroTotalPop, nroPop, nroVar, nroFObj, ...
                                              maxQLimiar, minQLimiar, ...
                                              maxScalingFactor, minScalingFactor, ...
                                              maxShiftConstant, minShiftConstant, nroSignals);
                      
chromosome = non_domination_sort_mod(chromosome, nroFObj, nroVar);

for i = 1 : nroG

    disp(['--Epoch: ' num2str(i)]);

    pool = round(nroPop/2);
    tour = 2;

    parent_chromosome = tournament_selection(chromosome, pool, tour);

    mu = 20;
    mum = 5;
    [offspring_chromosome, nroDict, dictList] = genetic_operator(nroDict,dictList,parent_chromosome, ...
                                            pDist, tabWaveletsComp, ...
                                            nroFObj, nroVar, mu, mum, ...
                                            [1 minQLimiar minScalingFactor minShiftConstant], ...
                                            [37 maxQLimiar maxScalingFactor maxShiftConstant], nroSignals);

    [main_pop,temp] = size(chromosome);
    [offspring_pop,temp] = size(offspring_chromosome);

    clear temp;

    intermediate_chromosome(1:main_pop,:) = chromosome;
    intermediate_chromosome(main_pop + 1 : main_pop + offspring_pop,1 : nroFObj+nroVar) = ...
        offspring_chromosome;

    intermediate_chromosome = ...
        non_domination_sort_mod(intermediate_chromosome, nroFObj, nroVar);

    chromosome = replace_chromosome(intermediate_chromosome, nroFObj, nroVar, nroPop);

disp(num2str(nroDict));

fileName = ['ArithFlickerPopulation' num2str(i) '.mat'];

save(fileName, 'chromosome');

end

figure(1);
plot(abs(chromosome(:,5)), chromosome(:,6),'*');
title('Results of the NSGA-II for signals with Flicker.');
xlabel('Compression Ratio');
ylabel('Distortion (nmse).');

pertWavelet = zeros(nroWavelets,1);
for i=1:nroWavelets
    for j=1:nroPop
        if round(chromosome(j,1))==i
            pertWavelet(i) = pertWavelet(i) + 1;
        end
    end
end

W=[];
for i=1:nroWavelets
    
    W = [W tabWaveletsComp(i).WaveletComp ': ' num2str(pertWavelet(i)) ' \n '];
    
end

save('Wavelets.mat', 'W');

save('DictList','dictList');

