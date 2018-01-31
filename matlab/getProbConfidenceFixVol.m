function getProbConfidenceFixVol(ks)
% getProbConfidenceFixVol([2 5 1 1 2 3 2 8 10])
% version: 1.0 beta
% authro: Yung-Kyun Noh (KAIST), 2014/08/28

if nargin < 1
    classnum = 2;
    ks = floor(rand(1,classnum)*10);
%     if size(ks, 2) ~= classnum
%         disp('****Invalide number of k-values****')
%         return
%     end
end
classnum = size(ks,2);

ps = zeros(1,classnum);


ncks = NChooseKdouble(classnum - 1);   % n choose ks
for iCurClass = 1:classnum
    curOtherClasses = [1:(iCurClass - 1) (iCurClass + 1):classnum];
    curProb = 1;
    for icntOtherClassNum = 1:classnum - 1
        totOtherClassIdxes = zeros(ncks(icntOtherClassNum + 1), icntOtherClassNum);
        totOtherClassIdxes = getOtherClassIdxes(curOtherClasses, icntOtherClassNum, [], totOtherClassIdxes, 0);

        for iProbClassSetIdx = 1:size(totOtherClassIdxes,1)
            OtherClassIdxesForProbCalculate = totOtherClassIdxes(iProbClassSetIdx,:);
            curClassK = ks(iCurClass);
            otherClassKs = ks(OtherClassIdxesForProbCalculate);
            curProb = curProb + (-1)^icntOtherClassNum*getProbCurClassLessThanProbOtherClasses(curClassK, otherClassKs, [], 0);
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
function sumProb = getProbCurClassLessThanProbOtherClasses(curClassK, otherClassKs, iJLs, sumProb)

Lval = size(otherClassKs,2) + 1;
if size(iJLs, 2) == size(otherClassKs,2)
    if curClassK + 1 + sum(otherClassKs,2) - sum(iJLs,2) < 0
        disp('*** negative val ***')
    end
    sumProb = sumProb + 1/Lval^(curClassK + 1 + sum(otherClassKs,2) - sum(iJLs,2))* ...
        multinomialCoeff(curClassK + sum(otherClassKs,2) - sum(iJLs,2),[curClassK otherClassKs - iJLs]);
    return  % return sumProb
end
for iCurJ = 0:otherClassKs(size(iJLs,2) + 1)
    sumProb = getProbCurClassLessThanProbOtherClasses(curClassK, otherClassKs, [iJLs iCurJ], sumProb);
end


% n choose ks
function multinomialCoeffVal = multinomialCoeff(nval, kvals)

karrays = zeros(1,sum(kvals,2));
curLoc = 0;
for icnt = 1:size(kvals,2)
    karrays(1,curLoc + [1:kvals(icnt)]) = 1:kvals(icnt);
    curLoc = curLoc + kvals(icnt);
end
multinomialCoeffVal = prod([1:nval]./karrays);


