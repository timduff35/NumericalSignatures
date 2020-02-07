setRandomSeed 0
path = append(path,"../")
needs "main.m2"
needs "defaults.m2"

RESULTSFILE = openOut("../experiment_data/joint-3pan-degree-"|toString(d))

dom = domain(d,4);
Map = jointEuclideanSigMap dom;
SlicePattern = new HashTable from for i from 0 to 3 list {i} => 1;
H = witnessHomotopy(dom,Map,SlicePattern);
rc = 8 * d^2 * (d^2-1);

load "master-experiment.m2"
close RESULTSFILE
end--

restart
for i from 2 to 6 do (
    d = i;
    << " on degree " << d << endl;
    elapsedTime load "joint-3pan.m2"
    )

