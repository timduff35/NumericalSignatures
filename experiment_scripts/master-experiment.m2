monodromyResult = timing runMonodromy(H, 
    NumberOfNodes=>NNODES, NumberOfEdges=>NEDGES, Verbose=>VERBOSE, 
    ImgTol=>IMGTOL, TargetSolutionCount=>rc, NumberOfRepeats=>NREPEATS, 
    SelectEdgeAndDirection=>selectBestEdgeAndDirection, Potential=>potentialE, Gammify=>GAMMIFY);
W = last monodromyResult
monodromyTime = first monodromyResult
-- top of file contains runtime information for witness collection
RESULTSFILE << round(1,monodromyTime) << "s,"
-- do we want this assert?
--assert(if (d>2) then length preimage W == rc else length image W == rc/4)

iters = OUTERITERS*INNERITERS
totWitTime = 0.0
witSucs = 0
curves = for i from 1 to OUTERITERS list (
    C0 := if SPHERE then sphere(binomial(d+2,2), 1);--else random(RR^1,RR^(binomial(d+2,2)));
    {C0} | for j from 2 to INNERITERS list randCoordChange(RR,d,C0)
    );
witnesses = for i from 0 to OUTERITERS-1 list (
    C0 := curves#i#0;
    if VERBOSE then << "collecting " << i << "th witness set on degree " << d << endl;
    witnessResult := timing witnessCollect(C0, W, ImgTol=>IMGTOL, Verbose=>VERBOSE);
    totWitTime = totWitTime + first witnessResult;
    W0 := last witnessResult;
    witSuc := if (d>2) then (length preimage W0 == rc) else (length image W0 == rc/4);
    if witSuc then witSucs = witSucs + 1;
    W0
    )
if (OUTERITERS > 0) then (
    avgWit = sub(totWitTime/OUTERITERS,RR);
    witSucRate = sub(witSucs/OUTERITERS,RR);
    RESULTSFILE << toString(round(1,avgWit)) << "s" << DELIM;
    RESULTSFILE << toString(pct(1,witSucRate)) << "%" << "\n";
    -- data frame header names
    RESULTSFILE << "noise,track,lookup,fail,FP,FN" << "\n";
    )


for NOISELEVEL in NOISELEVELS do(
    noiseyCurves := for i from 0 to OUTERITERS-1 list {curves#i#0} | for j from 1 to INNERITERS-1 list curves#i#j + sphere(binomial(d+2,2), NOISELEVEL);
    tests = for i from 0 to OUTERITERS-1 list(
    	-- test each curve against its equivalence class (#comparisons = INNERITERS)
    	pos := for j from 0 to INNERITERS-1 list equalityTest(noiseyCurves#i#j, witnesses#i, Verbose=>VERBOSE);
	-- test each curve against the other equivalence classes (#comparisons = )
	otherClasses := remove(toList(0..OUTERITERS-1),i);
	neg := for k from 0 to NUMNEGTESTS-1 list (
	    l := random(0,INNERITERS-1);
	    j := first random otherClasses;
	    equalityTest(noiseyCurves#j#l, witnesses#i, Verbose=>VERBOSE)
	    );
    	{pos,neg}
    	);
    allPos = flatten(tests/first);
    allNeg = flatten(tests/last);
    allTests = allPos | allNeg;
    sucPos = select(allPos,status);
    sucNeg = select(allNeg,status);
    sucNeg = select(allNeg,status);
    negSucRate  =sub((# select(allNeg,status))/#allNeg,RR);
    sucTests = select(allTests,status);
    sucRate = sub((# sucTests)/#allTests,RR);
    nTests = #allTests;
    avgTrack=1000*(sum(sucTests/trackTime)/nTests); -- average track time in milliseconds
    avgLookup=1000*(sum(sucTests/lookupTime)/nTests); -- average lookup time in milliseconds
    FNR = (#select(sucPos,t-> not result t))/(#sucPos);
    FPR = (#select(sucNeg,t-> result t))/(#sucNeg);
    RESULTSFILE << NOISELEVEL << DELIM;
-*
    -- reuse if we want to format tables
    RESULTSFILE << toString(round(1,avgLookup)) << "ms" << DELIM;
    RESULTSFILE << toString(pct(1,1-sucRate)) << "%" << DELIM;
    RESULTSFILE << toString(pct(1,FPR)) << "%" << DELIM;
    RESULTSFILE << toString(pct(1,FNR)) << "%" << "\n";
*-
    RESULTSFILE << toString(sub(avgTrack,RR)) << DELIM;
    RESULTSFILE << toString(sub(avgLookup,RR)) << DELIM;
    RESULTSFILE << toString(sub(1-sucRate,RR)) << DELIM;
    RESULTSFILE << toString(sub(FPR,RR)) << DELIM;
    RESULTSFILE << toString(sub(FNR,RR)) << "\n";
    

    )
