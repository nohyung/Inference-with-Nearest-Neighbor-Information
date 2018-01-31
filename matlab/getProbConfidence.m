function getProbConfidence(ks, ds, ns, dim)
% getProbConfidence([1 1 1], [.5 .5 .4], [100 100 100], 5)
% getProbConfidence([1 2 5], [1 5 2], [100 200 50], 5)
% version: 1.0 beta
% authro: Yung-Kyun Noh (KAIST), 2014/08/28

if nargin < 1
    classnum = 2;
    dim = 10;
    ks = ones(1,classnum);
    ds = sqrt(sum(rand(dim,classnum).^2, 1));
    ns = floor(rand(1,classnum)*10000);
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
            curClassD = ds(iCurClass);
            curClassN = ns(iCurClass);
            otherClassKs = ks(OtherClassIdxesForProbCalculate);
            otherClassDs = ds(OtherClassIdxesForProbCalculate);
            otherClassNs = ns(OtherClassIdxesForProbCalculate);
            curProb = curProb + ...
                (-1)^icntOtherClassNum*getProbCurClassLessThanProbOtherClasses(curClassK, curClassD, curClassN, ...
                otherClassKs, otherClassDs, otherClassNs, dim, [], 0);
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
    otherClassKs, otherClassDs, otherClassNs, dim, iJLs, sumProb)

Lval = size(otherClassKs,2) + 1;
if size(iJLs, 2) == size(otherClassKs,2)
    curClassU = curClassN*curClassD^dim;
    otherClassUs = otherClassNs.*otherClassDs.^dim;
    curClassScaledU = curClassU/sum([curClassU otherClassUs],2);
    otherClassScaledUs = otherClassUs/sum([curClassU otherClassUs],2);
    sumProb = sumProb + exp(sum([ ...
        (curClassK + 1)*log(curClassScaledU) iJLs.*log(otherClassScaledUs)], 2))* ...
        multinomialCoeff(curClassK + sum(iJLs,2),[curClassK iJLs]);
    return  % return sumProb
end
for iCurJ = 0:otherClassKs(size(iJLs,2) + 1)
    sumProb = getProbCurClassLessThanProbOtherClasses(curClassK, curClassD, curClassN, ...
        otherClassKs, otherClassDs, otherClassNs, dim, [iJLs iCurJ], sumProb);
end


% n choose ks
function multinomialCoeffVal = multinomialCoeff(nval, kvals)

karrays = zeros(1,sum(kvals,2));
curLoc = 0;
for icnt = 1:size(kvals,2)
    karrays(1,curLoc + [1:kvals(icnt)]) = 1:kvals(icnt);
    curLoc = curLoc + kvals(icnt);
end
% nval
% kvals
% karrays
multinomialCoeffVal = prod([1:nval]./karrays);


