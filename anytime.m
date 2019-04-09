function [] = anytimeArch(A,win,comp,exp,sweep,stop)

    [bsfarproper, timematproper,winc,wincm] = IterativeStomp(A,win,comp,exp,sweep);
    [bsfarrandom, timematrandom,winr,winrm] = PrunedIterativeStomp(A,win,comp,exp,sweep,stop);
    figure
    subplot(3,1,1)
    plot(timematproper,bsfarproper);hold on; plot(timematrandom,bsfarrandom)
    subplot(3,1,2)
    plot(zeroOneNorm(winc)); hold on; plot(zeroOneNorm(wincm));
    subplot(3,1,3)
    plot(zeroOneNorm(winr)); hold on; plot(zeroOneNorm(winrm));
end

function x = zeroOneNorm(x)
x = x - min(x(~isinf(x) & ~isnan(x)));
x = x / max(x(~isinf(x) & ~isnan(x)));

end