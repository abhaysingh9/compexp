function [PreMatrixProfile, PreMPindex] = PreSCRIMPABJoin(A,B , SubsequenceLength, step)
% Compute the approximate self similarity join of time series A with
% PreSCRIMP
% Author information omitted for ICDM review.
% For details of the SCRIMP++ algorithm, see:
% "SCRIMP++: Motif Discovery at Interactive Speeds", submitted to ICDM 2018.
% Usage:
% [PreMatrixProfile, PreMPindex] = PreSCRIMP(A, SubsequenceLength, step)
% Output:
%     PreMatrixProfile: running matrix profile after PreSCRIMP is completed
%     (vector)
%     PreMPindex: running matrix profile index after PreSCRIMP is completed
%     (vector)
% Input:
%     A: input time series (vector)
%     SubsequenceLength: interested subsequence length (scalar)
%     step: s/SubsequenceLength, the step size of PreSCRIMP. For all
%     experiments in the paper, this input is set to a fixed value 0.25.

%% set trivial match exclusion zone
exclusionZone = round(SubsequenceLength/4);

%% check input
...if SubsequenceLength > length(A)/2
    ...error('Error: Time series is too short relative to desired subsequence length');
....end
if SubsequenceLength < 4
    error('Error: Subsequence length must be at least 4');
end
if length(A) == size(A, 2)
   A = A'; 
end
if length(B) == size(B, 2)
   B = B'; 
end

%% initialization
MatrixProfileLength = length(A) - SubsequenceLength + 1;
MatrixProfile = zeros(MatrixProfileLength, 1);
MPindex = zeros(MatrixProfileLength, 1);
[XA, nA, sumx2A, sumxA, meanxA, sigmax2A, sigmaxA] = ...
    fastfindNNPre(A, SubsequenceLength);
[XB, nB, sumx2B, sumxB, meanxB, sigmax2B, sigmaxB] = ...
    fastfindNNPre(B, SubsequenceLength);


%% compute the matrix profile
current_step = floor(SubsequenceLength*step);
pickedIdx = 1:current_step:MatrixProfileLength;

dotproduct = zeros(MatrixProfileLength,1);
refine_distance = inf(MatrixProfileLength,1);
orig_index = 1:MatrixProfileLength;
m = SubsequenceLength;
for i = 1:1:length(pickedIdx)
    % compute the distance profile
    idx = pickedIdx(i);
    subsequenceB = B(idx:idx+SubsequenceLength-1);
    [distanceProfile, ~, ~, ~, ~] = fastfindNN(XA, subsequenceB, nA, SubsequenceLength, ...
            sumx2A, sumxA, meanxA, sigmax2A, sigmaxA);
    distanceProfile = abs(distanceProfile);
    
    % apply exclusion zone
    exclusionZoneStart = max(1, idx-exclusionZone);
    exclusionZoneEnd = min(MatrixProfileLength, idx+exclusionZone);
    distanceProfile(exclusionZoneStart:exclusionZoneEnd) = inf;
    
    % figure out and store the neareest neighbor
    if i == 1
        MatrixProfile = distanceProfile;
        MPindex(:) = idx;
        [MatrixProfile(idx), MPindex(idx)] = min(distanceProfile);
    else
        updatePos = distanceProfile < MatrixProfile;
        MPindex(updatePos) = idx;
        MatrixProfile(updatePos) = distanceProfile(updatePos);
        [MatrixProfile(idx), MPindex(idx)] = min(distanceProfile);
    end
    
    idx_nn = MPindex(idx);
    idx_diff = idx_nn-idx;
    dotproduct(idx) = (m - MatrixProfile(idx)^2/2)*sigmaxB(idx)*sigmaxA(idx_nn) + m*meanxB(idx)*meanxA(idx_nn);
    
    endidx = min([MatrixProfileLength, idx+current_step-1, MatrixProfileLength-idx_diff]);
    dotproduct(idx+1:endidx) =  dotproduct(idx)+cumsum(B(idx+m:endidx+m-1).*A(idx_nn+m:endidx+m-1+idx_diff)-B(idx:endidx-1).*A(idx_nn:endidx-1+idx_diff));
    refine_distance(idx+1:endidx) = sqrt(abs(2*(m-(dotproduct(idx+1:endidx)-m*meanxB(idx+1:endidx).*meanxA(idx_nn+1:endidx+idx_diff))./(sigmaxB(idx+1:endidx).*sigmaxA(idx_nn+1:endidx+idx_diff)))));
    
    beginidx = max([1, idx-current_step+1, 1-idx_diff]);
    dotproduct(idx-1:-1:beginidx) = dotproduct(idx)+cumsum(B(idx-1:-1:beginidx).*A(idx_nn-1:-1:beginidx+idx_diff)-B(idx-1+m:-1:beginidx+m).*A(idx_nn-1+m:-1:beginidx+idx_diff+m));
    refine_distance(beginidx:idx-1) = sqrt(abs(2*(m-(dotproduct(beginidx:idx-1)-m*meanxB(beginidx:idx-1).*meanxA(beginidx+idx_diff:idx_nn-1))./(sigmaxB(beginidx:idx-1).*sigmaxA(beginidx+idx_diff:idx_nn-1)))));
    
    update_pos1 = find(refine_distance(beginidx:endidx)<MatrixProfile(beginidx:endidx));
    MatrixProfile(update_pos1+beginidx-1) = refine_distance(update_pos1+beginidx-1);
    MPindex(update_pos1+beginidx-1) = orig_index(update_pos1+beginidx-1)+idx_diff;
    
    update_pos2 = find(refine_distance(beginidx:endidx)<MatrixProfile(beginidx+idx_diff:endidx+idx_diff));
    MatrixProfile(update_pos2+beginidx+idx_diff-1) = refine_distance(update_pos2+beginidx-1);
    MPindex(update_pos2+beginidx+idx_diff-1) = orig_index(update_pos2+beginidx+idx_diff-1)-idx_diff;
end

PreMatrixProfile = MatrixProfile;
PreMPindex = MPindex;
%{
[~,motifApos] = min(MatrixProfile);
    figure
    subplot(3,1,1)
    plot((1:1:length(A)),A)
    title('Time series data A in blue and B in red');
    hold on;
    plot((1:1:length(B)),B)
    hold on;
    plot((motifApos:1:motifApos+SubsequenceLength-1)',(A(motifApos:motifApos+SubsequenceLength-1)),'k')
    hold on;
    plot((MPindex(motifApos):1:MPindex(motifApos)+SubsequenceLength-1)',(B(MPindex(motifApos):MPindex(motifApos)+SubsequenceLength-1)),'k')
    hold on;
    subplot(3,1,2)   
    plot((1:1:length(MatrixProfile))',MatrixProfile)
    ylim([0 inf])
    title('Matrix Profile AB Join');
    [~,minind] = min(PreMatrixProfile);
    subplot(3,1,3)
    plot((1:1:SubsequenceLength),zeroOneNorm(A(minind:minind+SubsequenceLength-1)));
    hold on;
    plot((1:1:SubsequenceLength),zeroOneNorm(B(PreMPindex(minind):PreMPindex(minind)+SubsequenceLength-1)))
    title('Best Motif Discovered');
% m is winSize
%}
function [X, n, sumx2, sumx, meanx, sigmax2, sigmax] = fastfindNNPre(x, m)
n = length(x);
%x(n+1:2*n) = 0;
X = fft(x);
cum_sumx = cumsum(x);
cum_sumx2 =  cumsum(x.^2);
sumx2 = cum_sumx2(m:n)-[0;cum_sumx2(1:n-m)];
sumx = cum_sumx(m:n)-[0;cum_sumx(1:n-m)];
meanx = sumx./m;
sigmax2 = (sumx2./m)-(meanx.^2);
sigmax = sqrt(sigmax2);

% m is winSieze
function [dist lastz dropval sumy sumy2] = fastfindNN(X, y, n, m, sumx2, sumx, meanx, sigmax2, sigmax)
%x is the data, y is the query
%y = (y-mean(y))./std(y,1);                      %Normalize the query
dropval=y(1);
y = y(end:-1:1);                                %Reverse the query
y(m+1:n) = 0;

%The main trick of getting dot products in O(n log n) time
Y = fft(y);
Z = X.*Y;
z = ifft(Z);

%compute y stats -- O(n)
sumy = sum(y);
sumy2 = sum(y.^2);
meany=sumy/m;
sigmay2 = sumy2/m-meany^2;
sigmay = sqrt(sigmay2);

%computing the distances -- O(n) time
%dist = (sumx2 - 2*sumx.*meanx + m*(meanx.^2))./sigmax2 - 2*(z(m:n) - sumy.*meanx)./sigmax + sumy2;
%dist = 1-dist./(2*m);

dist = 2*(m-(z(m:n)-m*meanx*meany)./(sigmax*sigmay));
dist = sqrt(dist);
lastz=real(z(m:n));

function x = zeroOneNorm(x)
x = x - min(x(~isinf(x) & ~isnan(x)));
x = x / max(x(~isinf(x) & ~isnan(x)));
