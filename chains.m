function [sig] = chains()
    %rwa = cumsum(randn(500,1));
    %rwb = cumsum(randn(100,1));
    %rwc = cumsum(randn(500,1));
    rwa = 10*rand(500,1)-10;
    rwb = 10*rand(100,1)-10;
    rwc = 10*rand(500,1)-10;
    fs = 512;                    % Sampling frequency (samples per second)
    dt = 1/fs;                   % seconds per sample
    StopTime = 0.50;             % seconds
    t = (0:dt:StopTime-dt)';     % seconds
    F = 60;                      % Sine wave frequency (hertz)
    data1 = awgn(5*sin(2*pi*F*t),100,'measured');
    data2 = awgn(5*sin(2*pi*F*t*0.63),100,'measured');
    %sig = [data1; data2; rwa];
    sig = [rwa; rwa(500)+data1; rwa(500)+data1(StopTime/dt)+rwb; rwa(500)+data1(StopTime/dt)+rwb(100)+data2; rwa(500)+data1(StopTime/dt)+rwb(100)+data2(StopTime/dt)+rwc];
    %plot((1:1:length(sig)),sig);
    %tic
    %[c,d]=StompABJoin(sig(1:1.000001:(length(sig))),sig , length(data1));
    %toc
    %tic
    %[a,b]=PreSCRIMPABJoin(sig,sig, 250,0.25);
    %toc
    %[MP, MPi] = comexp(sig, 250, 0.60, 0.70, 0.005);
    plot(sig)
    anytime(sig, 100, 0.60,0.66, 0.03);
end