function pDisturbio = LoadPatterns(nroSignals)
 
k=1;
     for i=1:nroSignals

        path = strcat(strcat('SignalsDatabase\SignalDisturbanceType\Pattern', int2str(i)),'.mat');
       
        pDisturbio(k) = load(path);    
        
        k=k+1;
        
     end
 
end