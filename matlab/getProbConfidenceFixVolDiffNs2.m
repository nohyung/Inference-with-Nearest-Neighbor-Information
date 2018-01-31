function getProbConfidenceFixVolDiffNs2(ks, Ns)
% Src Code: 1409180001
% getProbConfidenceFixVolDiffNs2([2 5 1 1 2 3 2 8 10], [1000 1005 1010 1000 1000 900 1000 1000 1000])
% tic;getProbConfidenceFixVolDiffNs2([5 3 2 5 7 10], [500 1000 1005 1000 1000 800]);toc
% version: 1.0 beta
% author: Yung-Kyun Noh (KAIST), 2014/09/24

if nargin < 1
    classnum = 2;
    ks = floor(rand(1,classnum)*10);
    Ns = floor(rand(1,classnum)*1000);
end
if nargin < 2
    disp('getProbConfidenceFixVolDiffNs(ks, Ns): Ns is not specified')
    return
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
            otherClassKs = ks(OtherClassIdxesForProbCalculate);
            curClassN = Ns(iCurClass);
            otherClassNs = Ns(OtherClassIdxesForProbCalculate);
            curProb = curProb + ...
                (-1)^icntOtherClassNum*getProbCurClassLessThanProbOtherClasses(curClassK, otherClassKs, curClassN, otherClassNs, [], 0, 1);
        end

    end
    ps(iCurClass) = curProb;
end
ks
Ns
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
function sumProb = getProbCurClassLessThanProbOtherClasses(curClassK, otherClassKs, curClassN, otherClassNs, iJLs, sumProb, curMultinomialCoeff)

if size(iJLs, 2) == size(otherClassKs,2)
    Lvals = [curClassN otherClassNs]/sum([curClassN otherClassNs], 2);
    Lconst = exp([curClassK + 1  otherClassKs - iJLs].*log(Lvals));
    if curClassK + 1 + sum(otherClassKs,2) - sum(iJLs,2) < 0
        disp('*** negative val ***')
    end
    sumProb = sumProb + prod(Lconst,2)*curMultinomialCoeff;

    return  % return sumProb
end
% for iCurJ = 0:otherClassKs(size(iJLs,2) + 1)
for iCurJ = otherClassKs(size(iJLs,2) + 1):-1:0
    if isempty(iJLs)
        C1 = curClassK + sum(otherClassKs(1:(size(iJLs,2) + 1)), 2) - 0 - iCurJ;   % N in [N Choose k1, k2, ...]
        C2 = otherClassKs(1) - iCurJ;
    else
        C1 = curClassK + sum(otherClassKs(1:(size(iJLs,2) + 1)), 2) - sum(iJLs, 2) - iCurJ;   % N in [N Choose k1, k2, ...]
        C2 = otherClassKs(size(iJLs,2) + 1) - iCurJ;
    end
    if C2 == 0
        C1 = 1;
        C2 = 1;
    end
    curMultinomialCoeff = curMultinomialCoeff*C1/C2;
    sumProb = getProbCurClassLessThanProbOtherClasses(curClassK, otherClassKs, curClassN, otherClassNs, [iJLs iCurJ], sumProb, curMultinomialCoeff);
end

