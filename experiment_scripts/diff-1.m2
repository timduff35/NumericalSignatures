setRandomSeed 0
path = append(path,"../")
needs "main.m2"
needs "defaults.m2"

RESULTSFILE = openOut("../experiment_data/diff-1-degree-"|toString(d))

-- override defaults here

dom = domain(d,1);
Map = diffEuclideanSigMap dom;
SlicePattern = new HashTable from {{0}=>1};
H = witnessHomotopy(dom,Map,SlicePattern);
rc = 6*d*(d-1);
load "master-experiment.m2"
close RESULTSFILE
end--
restart
for i from 2 to 6 do (
    d = i;
    << " on degree " << d << endl;
    elapsedTime load "diff-1.m2"
    )
