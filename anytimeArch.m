function [] = anytime(A,win,comp,exp,sweep)
    % brute force

    [bsfarproper, timematproper,winc,wincm] = comexpany(A,win,comp,exp,sweep);
    [bsfarrandom, timematrandom,winr,winrm] = randomany(A,win,comp,exp,sweep);
    figure
    subplot(3,1,1)
    plot(timematproper,bsfarproper, '-o');hold on; plot(timematrandom,bsfarrandom)
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
    for i = 1:length(A)-win+1
        winr = randi([1,length(winsel)],1,1);
        winsel(winr) = []; 
        temp = MASS(cMat,A(winr:winr+win-1));
        [bestsofartemp, indexBest] = min(minP(A,temp,winr,win,comp,sweep,exp)); % To be changed with min of non excluded
        if(bestsofartemp < bestsofar(timec))
            timec=timec+1;
            bestsofar(timec)=bestsofartemp;
            timerandom(timec) = toc;
            win1 = A(winr:winr+win-1);
            win2 = cMat(indexBest:indexBest+win-1);
        end
        
        if(timerandom(timec)-toc > 0.01)
            timec=timec+1;
            bestsofar(timec)=bestsofartemp;
            timerandom(timec) = toc;
        end
    end
end


function [temp] = minP(A,temp,winr,win,comp,sweep,exp)
    lensofar = 0;
    for j = comp:sweep:exp
        if(j~=1)
            point = ceil(lensofar + winr/j);
            paddingleft = -win + 1;
            paddingright = win - 1;
            if(point + paddingleft < 1)
               paddingleft = -point + 1; 
            end
            curlen = length(A(1:j:end));
            if(point + paddingright > lensofar+ curlen)
               paddingright = lensofar + curlen - point + 1;
            end
            for i = point+paddingleft:1:point+paddingright
                temp(i) = inf; 
            end
            lensofar = lensofar + curlen;
        end
    end
end



function [bsfar,times,win1,win2] = comexpany(A,win,comp,exp,sweep)
    tic;
    toc
    bsfar = inf;                     
    times = 0;
    warning('off','all');
    counter = 1;
    timec = 1;
    for i = comp:sweep:exp 
        if(i~=1)
            [MatrixProfile, MPindex] = ...
                StompABJoin(A,A(1:i:length(A)),win);
            [MPval,indexMP] = findvaluelow(MatrixProfile,MPindex, win, i);... min(MatrixProfile);
            counter = counter + 1;
            bsfartemp = min(MPval);
            if(bsfartemp<bsfar(timec)) 
                timec = timec+1;
                bsfar(timec) = bsfartemp;
                times(timec) = toc;
                win1 = A(indexMP:indexMP+win-1);
                temp = A(1:i:length(A));
                win2 = temp(MPindex(indexMP):MPindex(indexMP)+win-1);
            end
            if(times(timec)-toc > 0.01)
                bsfar(timec)=bestsofartemp;
                times(timec) = toc;
                timec=timec+1;
            end
        end
    end
end

function [val, ind] = findvaluelow(MatrixProfile,MPi,win,rate)

    [e,indx] = min(MatrixProfile);
    while(((indx/rate) > MPi(indx)-win) & ((indx/rate) < (MPi(indx)+win-1)))
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