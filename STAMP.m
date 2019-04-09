function [bestsofar, timerandom, win1, win2] = STAMP(A,win,comp,exp,sweep,stop)
    tic;
    toc
    cMat = [];
    bestsofar = inf;
    timec = 1;
    for i = comp:sweep:exp
        if (i~=1)
            check = cMat';
            tempee = A(1:i:length(A))+1000;
            cMat = [check, tempee'];
            cMat = cMat';
        end
    end
    winsel = 1:1:(length(A)-win+1);
    % Algorithm Random
    completion = 0;
    collection = [];
    while(length(winsel) ~= 0)
        winr = winsel(randi([1,length(winsel)],1,1));
        winsel(winsel==winr) = []; 
        collection = [collection winr];
        temp = MASS(cMat,A(winr:winr+win-1));
        [bestsofartemp, indexBest] = min(minP(A,temp,winr,win,comp,sweep,exp)); % To be changed with min of non excluded
        if(bestsofartemp < bestsofar(timec))
            timec=timec+1;
            bestsofar(timec)=bestsofartemp;
            timerandom(timec) = toc;
            win1 = A(winr:winr+win-1);
            win2 = cMat(indexBest:indexBest+win-1);
        end
        
        if(toc - timerandom(timec) > 0.01)
            timec=timec+1;
            bestsofar(timec)=bestsofar(timec - 1);
            timerandom(timec) = toc;
        end
        completion = (length(A)-length(winsel)-win +1)/(length(A)-win +1);
        if(completion>stop)
           return 
        end
    end
end


function [temp] = minP(A,temp,winr,win,comp,sweep,exp)
  %{
    len = 0; counter = 1;
    for j = comp:sweep:exp
        forbidden(counter) = winr/j + len ;
        len = length(A(1:j:end));
        counter = counter + 1;
    end
    [val, minind] = min(temp);
    while(((minind) > (forbidden - win)) && ((minind )> (forbidden - win)))
       if(val == inf) 
           return;
       end
       temp(minind) = inf;
       [val,minind] = min(temp);
    end
    
%}
    lensofar = 0;
    for j = comp:sweep:exp
        if(j~=1)
            point = lensofar + ceil(winr/j);
            paddingleft = - win + 1 ;
            paddingright = win - 1;
            if( point + paddingleft < 1)
               paddingleft = -point ; 
            end
            curlen = length(A(1:j:end));
            if(point + paddingright > length(temp))
               paddingright = length(temp) - point ;
            end
            temp(point+paddingleft+1:point+paddingright-1) = inf;
            lensofar = lensofar + curlen;
        end
    end
    %}
end