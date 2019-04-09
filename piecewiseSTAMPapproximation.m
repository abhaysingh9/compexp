function [bestsofar, timerandom, win1, win2] = piecewiseSTAMPapproximation(A,win,comp,exp,sweep,stop)
    tic;
    bestsofar = inf;
    counter = 1;
    i = comp;   
    winlast = [];
    while(i < exp)
        completion = 0;
        if (i~=1)
            cMat = A(1:i:length(A)) ;
            cMat = [cMat; winlast];
            collection = [];
            winsel = 1:1:length(A)-win +1;
            while(length(winsel) ~= 0 &&  completion < stop)
                winr = winsel(randi([1,length(winsel)],1,1));
                winsel(winsel==winr) = []; 
                collection = [collection winr];
                temp = MASS(cMat,A(winr:winr+win-1));
                padding = win;
                if((winr/i)-win<1)
                    padding = winr/i -1;
                end
                temp(((winr/i)-padding):(winr/i)+win) = inf;
                [bestsofartemp, indexBest] = min(temp); % To be changed with min of non excluded
                if(bestsofartemp < bestsofar(counter))
                    counter=counter+1;
                    bestsofar(counter)=bestsofartemp;
                    timerandom(counter) = toc;
                    win1 = A(winr:winr+win-1);
                    win2 = cMat(indexBest:indexBest+win-1);
                end
        
                if(toc - timerandom(counter) > 0.01)
                    counter=counter+1;
                    bestsofar(counter)=bestsofar(counter - 1);
                    timerandom(counter) = toc;
                end
                completion = (length(A)-length(winsel)-win +1)/(length(A)-win +1);
            end
        end
        i = i + sweep;
        winlast = win2;
    end
    
end