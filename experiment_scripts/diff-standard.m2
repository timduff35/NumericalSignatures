restart
setRandomSeed 0
path = append(path,"../")
needs "main.m2"
needs "defaults.m2"

RESULTSFILE = openOut("../experiment_data/diff-standard-degree-"|toString(d))

dom = domain(d,1);
Map = diffEuclideanSigMap dom;
H = witnessHomotopy(dom,Map);
rc = 12*(d^2-d);

load "master-experiment.m2";
close RESULTSFILE

end--
restart
for i from 2 to 6 do (
    d = i;
    << " on degree " << d << endl;
    elapsedTime load "diff-standard.m2"
    )
