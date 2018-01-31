function getProbConfidence2(ks, ds, ns, dim)
% Src Code: 1409180003
% getProbConfidence2([1 1 1], [.5 .5 .4], [100 100 100], 5)
% tic;getProbConfidence2([4 3 5 3 5], [.5 .5 .4 .2 .3], [100 100 100 200 150], 5);toc
% tic;getProbConfidence2([4 3 5 3 5 3 2], [.5 .5 .4 .2 .3 .3 .7], [100 100 100 200 150 200 180], 5);toc  % 9.19 sec in desktop, 53.35 sec in laptop
% version: 1.0 beta
% authro: Yung-Kyun Noh (KAIST), 2014/09/24


if nargin < 1
    classnum = 2;
    dim = 10;
    ks = ones(1,classnum);
    ds = sqrt(sum(rand(dim,classnum).^2, 1));
    ns = floor(rand(1,classnum)*10000);
end
classnum = size(ks,2);

ps = zeros(1,classnum);


for iCurClass = 1:classnum
    curOtherClasses = [1:(iCurClass - 1) (iCurClass + 1):classnum];
    curProb = 1;
    for icntOtherClassNum = 1:classnum - 1
        totOtherClassIdxes = nchoosek(curOtherClasses, icntOtherClassNum);   % for totOtherClassIdxes #2

        for iProbClassSetIdx = 1:size(totOtherClassIdxes,1)
            OtherClassIdxesForProbCalculate = totOtherClassIdxes(iProbClassSetIdx,:);
            curClassK = ks(iCurClass);
            curClassD = ds(iCurClass);
            curClassN = ns(iCurClass);
            otherClassKs = ks(OtherClassIdxesForProbCalculate);
            otherClassDs = ds(OtherClassIdxesForProbCalculate);
            otherClassNs = ns(OtherClassIdxesForProbCalculate);
            curProb = curProb + ...
                (-1)^icntOtherClassNum*getProbCurClassLessThanProbOtherClasses(curClassK, curClassD, curClassN, ...
                otherClassKs, otherClassDs, otherClassNs, dim, [], 0, 1);
        end
        
    end
    ps(iCurClass) = curProb;
end
ks
ds
ns
dim
ps


%%%%%%%%%%%%%%%%%
function [totOtherClasses, filledNum] = getOtherClassIdxes(remainedClasses, ClassNumLeftToSelect, selectedClasses, totOtherClasses, filledNum)

if ClassNumLeftToSelect == 0
    filledNum = filledNum + 1;
    totOtherClasses(filledNum,:) = selectedClasses;
    return   % fill one more line
end

for iRemainedClassIdx = 1:size(remainedClasses,2)
    inputRemainedClasses = remainedClasses((iRemainedClassIdx + 1):end);
    intpuClassNumLeftToSelect = ClassNumLeftToSelect - 1;
    inputSelectedClasses = [selectedClasses remainedClasses(iRemainedClassIdx)];
    [totOtherClasses, filledNum] = getOtherClassIdxes(inputRemainedClasses, intpuClassNumLeftToSelect, inputSelectedClasses, totOtherClasses, filledNum);
end    


%%%%%%%%%%%%%%%%%
function sumProb = getProbCurClassLessThanProbOtherClasses(curClassK, curClassD, curClassN, ...
    otherClassKs, otherClassDs, otherClassNs, dim, iJLs, sumProb, curMultinomialCoeff)

Lval = size(otherClassKs,2) + 1;
if size(iJLs, 2) == size(otherClassKs,2)
    curClassU = curClassN*curClassD^dim;
    otherClassUs = otherClassNs.*otherClassDs.^dim;
    curClassScaledU = curClassU/sum([curClassU otherClassUs],2);
    otherClassScaledUs = otherClassUs/sum([curClassU otherClassUs],2);
    sumProb = sumProb + exp(sum([ ...
        (curClassK + 1)*log(curClassScaledU) iJLs.*log(otherClassScaledUs)], 2))* ...
        curMultinomialCoeff;
    return  % return sumProb
end
for iCurJ = 0:otherClassKs(size(iJLs,2) + 1)
% for iCurJ = otherClassKs(size(iJLs,2) + 1):-1:0
    if isempty(iJLs)
        C1 = curClassK + 0 + iCurJ;   % N in [N Choose k1, k2, ...]
        C2 = iCurJ;
    else
        C1 = curClassK + sum(iJLs, 2) + iCurJ;   % N in [N Choose k1, k2, ...]
        C2 = iCurJ;
    end
    if C2 == 0
        C1 = 1;
        C2 = 1;
    end
    curMultinomialCoeff = curMultinomialCoeff*C1/C2;
    sumProb = getProbCurClassLessThanProbOtherClasses(curClassK, curClassD, curClassN, ...
        otherClassKs, otherClassDs, otherClassNs, dim, [iJLs iCurJ], sumProb, curMultinomialCoeff);
end


