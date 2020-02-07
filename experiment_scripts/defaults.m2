-- path tracker options
--setDefault(tStepMin=>1e-8)
--setDefault(maxCorrSteps=>2)
--setDefault(CorrectorTolerance=>1e-10)
-- experiment configuration globals
VERBOSE=false
NNODES=2
NEDGES=5
MAXDEG=6
NREPEATS=10
DELIM = "," 
IMGTOL=1e-6---points in the image at distance<IMGTOL are treated as equal (important only for witness collection. tolerance for equality test is handled by PointArray tolerance)
OUTERITERS = 10 -- how many random curves to generate
INNERITERS = 20 -- how many curves in the orbit of each random curve? used to assess sensitivity / false negative rate
NUMNEGTESTS=5 -- for monitoring the false positive rate (seems low, so we don't focus on doing many of these tests)
NOISELEVELS = {0.0} | for k from -7 to -3 list 10.0^k
SPHERE = true
DBG=0
GAMMIFY = true