%% script for testing offset to onset lag

%what do we need - we  t_j,s (the delay between learning onset and memorisation for 
%a each pair id )

%this is defined specifically for rcrit = 0.95

%so we first we need the point of memorisation  for a subject so we do

%do we loop through pC for every trial untill we find the first r = 0.95
%but we're only doing this for nonzero sup pairs

function [modelObj] = z04_Estimate_A00(G)

xMem = getXmem(G)

 %then we get all the onset terms for each question and use them to predict
 %xStartLearning for each

 %then you have th difference term for all question types and pair ids by
 %doing...

xStartLearn = getXStrtLrn(G);

 Delta = xMem - xStartLearn;
 

%then we collect up a big matrix of delta and all the binary codes of the
%predictors into a a table

table = table(Delta,Learnt,Zero,Commute,Noncom,PairId);
 %then we model delta 

modelObj = fitlme(table,'delta ~ 1 + learnt * (zero + commute + noncom) + (1|PairId) + (1|PairId)');


return 