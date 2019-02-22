function [] = anytime(A,win,comp,exp,sweep)
    % brute force
    
    [bsfarproper, timematproper] = comexpany(A,win,comp,exp,sweep);
    [bsfarrandom, timematrandom] = randomany(A,win,comp,exp,sweep);
    figure
    plot(timematproper,bsfarproper, '-o');hold on; plot(timematrandom,bsfarrandom) 
end

function [bestsofar, timerandom] = randomany(A,win,comp,exp,sweep)
    tic;
    toc
    cMat = [];
    bestsofar = inf;
    timec = 0;
    for i = comp:sweep:exp
        if (i~=1)
            check = cMat';
            tempee = A(1:i:length(A));
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
        bestsofartemp = min(temp); % To be changed with min of non excluded
        if(bestsofartemp<bestsofar)
            timec=timec+1;
            bestsofar(timec)=bestsofartemp;
            timerandom(timec) = toc;
        end
        
        if(timerandom(timec)-toc > 0.01)
            timec=timec+1;
            bestsofar(timec)=bestsofartemp;
            timerandom(timec) = toc;
        end
    end
end


function [sprofile] = consolidate(tempmin)

end



function [bsfar,times] = comexpany(A,win,comp,exp,sweep)
    tic;
    toc
    bsfar = inf;                     
    times = 0;
    warning('off','all');
    if(mod(1-comp,sweep)<0.0000000000000001)
        Profilelength = length(comp:sweep:exp)-1; 
    else 
        Profilelength = length(comp:sweep:exp); 
    end
    MatrixProfileLength = Profilelength; 
    MP = zeros(MatrixProfileLength, 1);
    MPi = zeros(MatrixProfileLength, 3);
    counter = 1;
    timec = 0;
    for i = comp:sweep:exp 
        if(i~=1)
            [MatrixProfile, MPindex] = ...
                StompABJoin(A,A(1:i:length(A)),win);
            [MP,~] = findvaluelow(MatrixProfile,MPindex, win, i);... min(MatrixProfile);
            counter = counter + 1;
            bsfartemp = min(MP);
            if(bsfartemp<bsfar) 
                timec = timec+1;
                bsfar(timec) = bsfartemp;
                times(timec) = toc;
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