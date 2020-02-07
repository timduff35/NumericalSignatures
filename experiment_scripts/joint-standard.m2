setRandomSeed 0
path = append(path,"../")
needs "main.m2"
needs "defaults.m2"

RESULTSFILE = openOut("../experiment_data/joint-standard-degree-"|toString(d))

dom = domain(d,4);
Map = jointEuclideanSigMap dom;
H = witnessHomotopy(dom,Map);
rc = 12 * d * (d^3-1);

load "master-experiment.m2";
close RESULTSFILE

end--

restart
for i from 2 to 6 do (
    d = i;
    << " on degree " << d << endl;
    elapsedTime load "joint-standard.m2"
    )