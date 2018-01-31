function getProbConfidenceFixVol2(ks)
% Src Code: 1409180002
% getProbConfidenceFixVol2([2 5 1 1 2 3 2 8 10])    % 143 sec in desktop, 835 sec in labtop
% tic;getProbConfidenceFixVol2([5 3 2 5 7 10]);toc



if nargin < 1
    classnum = 2;
    ks = floor(rand(1,classnum)*10);
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
            curProb = curProb + (-1)^icntOtherClassNum*getProbCurClassLessThanProbOtherClasses(curClassK, otherClassKs, [], 0, 1);
        end
        
    end
    ps(iCurClass) = curProb;
end
ks
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
function sumProb = getProbCurClassLessThanProbOtherClasses(curClassK, otherClassKs, iJLs, sumProb, curMultinomialCoeff)

if size(iJLs, 2) == size(otherClassKs,2)
    Lval = size(otherClassKs,2) + 1;    % number of classes
    if curClassK + 1 + sum(otherClassKs,2) - sum(iJLs,2) < 0
        disp('*** negative val ***')
    end
    sumProb = sumProb + 1/Lval^(curClassK + 1 + sum(otherClassKs,2) - sum(iJLs,2))*curMultinomialCoeff;

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
    sumProb = getProbCurClassLessThanProbOtherClasses(curClassK, otherClassKs, [iJLs iCurJ], sumProb, curMultinomialCoeff);
end

