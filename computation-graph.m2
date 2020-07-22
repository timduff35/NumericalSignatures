needs "main.m2"
debug SLPexpressions
needsPackage "Visualize"
needsPackage "Graphs"
PolyGraph = new Type of MutableHashTable
newPolyGraph = () -> (h := new PolyGraph; 
    h#"Graph" = digraph{}; 
    h#"#consts"=h#"#vars"=h#"#gates"=h#"#lines"=0; 
    h#"V" = new MutableList from {}; h#"nV" = 0;
    h#"E" = new MutableList from {}; h#"nE" = 0;
    h)
addLine (PolyGraph, Thing) := (h,t) -> ( h#(h#"#lines") = t; h#"#lines" = h#"#lines" + 1; )    
addUnit = method()
addUnit (Gate, PolyGraph) := (g,h) -> error "not implemented"
addUnit (InputGate, PolyGraph) := (g,h) -> if h#?g then h#g else (
    if isConstant g then (
	h#g = "C"|toString h#"#consts";
    	h#"#consts" = h#"#consts" + 1;
	)
    else (
	h#g = "X"|toString h#"#vars";
    	h#"#vars" = h#"#vars" + 1;
	);
    h#"V"#(h#"nV") = h#g;
    h#"nV" = h#"nV" + 1;
    )
addUnit (SumGate, PolyGraph) := (g,h) -> if h#?g then h#g else (
    s := between(" + ", apply(g.Inputs, gg->addUnit(gg,h)));  
    h#g = "G"|toString h#"#gates";
    h#"#gates" = h#"#gates" + 1;
    h#"V"#(h#"nV") = h#g;
    h#"nV" = h#"nV" + 1;
    scan(g.Inputs, i -> (e := {h#i, h#g}; h#"E"#(h#"nE") = e; h#"nE" = h#"nE" + 1;));
    )
addUnit (ProductGate, PolyGraph) := (g,h) -> if h#?g then h#g else (
    s := between(" * ", apply(g.Inputs, gg->addUnit(gg,h)));  
    h#g = "G"|toString h#"#gates";
    h#"#gates" = h#"#gates" + 1;
    h#"V"#(h#"nV") = h#g;
    h#"nV" = h#"nV" + 1;
    scan(g.Inputs, i -> (e := {h#i, h#g}; h#"E"#(h#"nE") = e; h#"nE" = h#"nE" + 1;));
    )
addUnit (DivideGate, PolyGraph) := (g,h) -> if h#?g then h#g else (
    (x,y) := toSequence apply(g.Inputs, gg->addUnit(gg,h));  
    h#g = "G"|toString h#"#gates";
    h#"#gates" = h#"#gates" + 1;
    h#"V"#(h#"nV") = h#g;
    h#"nV" = h#"nV" + 1;
    scan(g.Inputs, i -> (e := {h#i, h#g}; h#"E"#(h#"nE") = e; h#"nE" = h#"nE" + 1;));
    )

makePolyGraph = (inputs, outputs) -> (
    h := newPolyGraph();
    scan(inputs, g -> addUnit(g,h));
    scan(outputs, g -> addUnit(g,h));
    h#"Graph" = digraph(toList h#"V", toList h#"E");
    h#"Graph"
    )

end--
restart
needs "computation-graph.m2"
declareVariable z;
sq=z^2;
f = (sq+1)^2+sq;
printAsSLP({z}, {f})
inputs = {z}
outputs = {f}
G = makePolyGraph(inputs, outputs)
openPort"8080"
visualize G


d=6
dom = domain(d, 1)
Map = diffEuclideanSigMap dom
(I,O)=(flatten entries(vars Map|(parameters Map)_{0..binomial(2+d,2)-1}), flatten entries gateMatrix Map);
printAsSLP(I,O)
elapsedTime G=makePolyGraph(I,O);
visualize G
