function [] = anytimeArch(A,win,comp,exp,sweep)
    % brute force

    [bsfarproper, timematproper,winc,wincm] = comexpany(A,win,comp,exp,sweep);
    [bsfarrandom, timematrandom,winr,winrm] = randomany(A,win,comp,exp,sweep);
    figure
    subplot(3,1,1)
    plot(timematproper,bsfarproper);hold on; plot(timematrandom,bsfarrandom)
    subplot(3,1,2)
    plot(zeroOneNorm(winc)); hold on; plot(zeroOneNorm(wincm));
    subplot(3,1,3)
    plot(zeroOneNorm(winr)); hold on; plot(zeroOneNorm(winrm));
end

function [bestsofar, timerandom, win1, win2] = randomany(A,win,comp,exp,sweep)
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



function [bsfar,times,win1,win2] = comexpany(A,win,comp,exp,sweep)
    tic;
    toc
    bsfar = inf;                     
    times = 0;
    warning('off','all');
    counter = 1;
    timec = 1;
    index = inf;
    for i = comp:sweep:exp 
        if(i~=1)
            [MatrixProfile, MPindex] = ...
                StompABJoin(A,A(1:i:length(A)),win);
            [MPval,indexMP] = findvaluelow(MatrixProfile,MPindex, win, i);... min(MatrixProfile);
            counter = counter + 1;
            bsfartemp = min(MPval);
            if(bsfartemp<=bsfar(timec)) 
                timec = timec+1;
                bsfar(timec) = bsfartemp;
                times(timec) = toc;
                win1 = A(indexMP:indexMP+win-1);
                temp = A(1:i:length(A));
                win2 = temp(MPindex(indexMP):MPindex(indexMP)+win-1);
                index = indexMP; 
            end
            if(toc - times(timec) > 0.01)
                timec=timec+1;
                bsfar(timec)=bsfar(timec-1);
                times(timec) = toc;
                
            end
        end
    end
end

function [val, ind] = findvaluelow(MatrixProfile,MPi,win,rate)
%{
    for i = (1:1:MatrixProfile)
        if(ceil(i/rate) > MPi(i) - win +1 && ceil(i/rate)< MPi(i)+win -1 )
            MatrixProfile(i-win+1:i+win-1) = inf;  
        end
    end
%}
    [e,indx] = min(MatrixProfile);
    while((ceil(indx/rate) > MPi(indx)-win/2+1) && (ceil(indx/rate) < (MPi(indx)+win/2-1)))
        if(e == inf) 
            val = e; ind = indx;
            return;
        end
       MatrixProfile(indx) = inf;
       [e,indx] = min(MatrixProfile);
    end
    ind = indx;
    val = MatrixProfile(ind);
end


function x = zeroOneNorm(x)
x = x - min(x(~isinf(x) & ~isnan(x)));
x = x / max(x(~isinf(x) & ~isnan(x)));

end