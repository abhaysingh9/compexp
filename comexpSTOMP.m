function [MP, MPi] = comexpSTOMP(A,win,comp,exp,sweep)
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
    for i = comp:sweep:exp 
        if(i~=1)
            [MatrixProfile, MPindex] = ...
                StompABJoin(A,A(1:i:length(A)),win);
            [MP(counter),MPi(counter,1)] = findvaluelow(MatrixProfile,MPindex, win, i);... min(MatrixProfile);
            MPi(counter,2) = MPindex(MPi(counter,1));
            MPi(counter,3) = i;
            counter = counter + 1;
                i
        end
    end
    mp = MP;
    [~,minind] = min(mp);
    figure
    subplot(5,1,1)
    plot((1:1:length(A)),A);
    title('Time series being analyzed');
    subplot(5,1,2)
    pure = (comp:sweep:exp);
    pure(pure == 1) = []; 
    plot(pure,MP)
    ylim([0 inf])
    title('Matrix sampling profile');
    subplot(5,1,3)
    plot((1:1:win),zeroOneNorm(A(MPi(minind,1):MPi(minind,1)+win-1)))
    hold on;
    temp = A(1:MPi(minind,3):length(A));
    plot((1:1:win),zeroOneNorm(temp(MPi(minind,2):MPi(minind,2)+win-1)))
    string = strcat('1st Best Motif at sampling ratio: ',num2str(MPi(minind,3)),', MASS distance:', num2str(MP(minind)));
    title(string);
    mp(minind) = inf;
    [~,minind] = min(mp);
    subplot(5,1,4)
    plot((1:1:win),zeroOneNorm(A(MPi(minind,1):MPi(minind,1)+win-1)))
    hold on;
    temp = A(1:MPi(minind,3):length(A));
    plot((1:1:win),zeroOneNorm(temp(MPi(minind,2):MPi(minind,2)+win-1)))
    string = strcat('2nd Best Motif at:',num2str(MPi(minind,3)),', distance:', num2str(MP(minind)));
    title(string);
    subplot(5,1,5)
    mp(minind) = inf;
    [~,minind] = min(mp);
    plot((1:1:win),zeroOneNorm(A(MPi(minind,1):MPi(minind,1)+win-1)))
    hold on;
    temp = A(1:MPi(minind,3):length(A));
    plot((1:1:win),zeroOneNorm(temp(MPi(minind,2):MPi(minind,2)+win-1)))
    string = strcat('3rd Best Motif at:',num2str(MPi(minind,3)),', distance:', num2str(MP(minind)));
    title(string);
    %{
    subplot(6,1,6)
    mp = MP;
    [~,maxind] = max(mp);
    plot((1:1:win),zeroOneNorm(A(MPi(maxind,1):MPi(maxind,1)+win-1)))
    hold on;
    mp(maxind) = -1;
    [~,maxind] = max(mp);
    plot((1:1:win),zeroOneNorm(A(MPi(maxind,1):MPi(maxind,1)+win-1)))
    hold on;
    mp(maxind) = -1;
    [~,maxind] = max(mp);
    plot((1:1:win),zeroOneNorm(A(MPi(maxind,1):MPi(maxind,1)+win-1)))
    title('Greatest Discords'); 
    %}
end

function x = zeroOneNorm(x)
x = x - min(x(~isinf(x) & ~isnan(x)));
x = x / max(x(~isinf(x) & ~isnan(x)));
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