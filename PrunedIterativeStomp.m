function [bsfar,times,win1,win2] = PrunedIterativeStomp(A,win,comp,exp,sweep,stop)
    tic;
    warning('off','all');
    [bsfar, times, win1,win2] = STAMP(A,win,comp,exp,sweep,stop);
    bsfartemp = bsfar(end);
    counter = length(bsfar);
    for i = comp:sweep:exp 
        if(i~=1)
            [MatrixProfile, MPindex] = ...
                StompABJoin(A,A(1:i:length(A)),win);
            [MPval,indexMP] = findvaluelow(MatrixProfile,MPindex, win, i);... min(MatrixProfile);
            bsfartemp = min(MPval);
            if(bsfartemp<=bsfar(counter)) 
                counter = counter+1;
                bsfar(counter) = bsfartemp;
                times(counter) = toc;
                win1 = A(indexMP:indexMP+win-1);
                temp = A(1:i:length(A));
                win2 = temp(MPindex(indexMP):MPindex(indexMP)+win-1);
                index = indexMP; 
            end
            if(toc - times(counter) > 0.01)
                counter=counter+1;
                bsfar(counter)=bsfar(counter-1);
                times(counter) = toc;
                
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