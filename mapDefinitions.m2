-- these should be rewritten
-- SLPs for signature maps: feed into 'setupMonodromy' func

diffEuclideanSigMap = dom -> (
    domGS := gateSystem dom;
    f := (gateMatrix domGS)_(0,0);
    mapVars := vars domGS;
    X := mapVars_(0,0);
    Y := mapVars_(0,1);
    Z := mapVars_(0,2);
    y1 := - compress((compress diff(X, f))/(compress diff(Y,f)));
    y2 := compress(compress(diff(X, y1)) + compress(compress(diff(Y, y1)) * y1));
    y3 := compress(diff(X, y2) + diff(Y, y2) * y1);
    T1 := (1 + y1^2);
    T2 := y2;
    T2sq := T2*T2;
    subExp1 := Z^2/T1^3;
    T3 := y3 * T1 - 3 * y1 *T2sq;
    subExp2 := subExp1 * T3;
    gateSystem(
	parameters domGS, 
	mapVars, 
	transpose gateMatrix{{subExp1*T2sq, subExp2^2}}--Z^4*T3^2/T1^6}} -- toggle to get SE(2)
	)
    )

jointDiffEuclideanSigMap = dom -> (
    domGS := gateSystem dom;
    f1 := (gateMatrix domGS)_(0,0);
    f2 := (gateMatrix domGS)_(2,0);
    mapVars := vars domGS;
    X1 := mapVars_(0,0);
    Y1 := mapVars_(0,1);
    Z1 := mapVars_(0,2);
    X2 := mapVars_(0,3);
    Y2 := mapVars_(0,4);
    Z2 := mapVars_(0,5);
    fx1 := compress diff(X1, f1);
    fx2 := compress diff(X2, f2);
    fy1 := compress diff(Y1, f1);
    fy2 := compress diff(Y2, f2);
    -- implicit derivatives up to sign
    y1 := fx1 / fy1;
    y2 := fx2 / fy2;
    a := (X1*Z2-X2*Z1)/(Y2*Z1-Y1*Z2);
-*
    n1 := compress (Y2/Z2-Y1/Z1)*fx1+(X1/Z1-X2/Z2)*fy1;
    d1 := compress (X1/Z1-X2/Z2)*fx1+(Y1/Z1-Y2/Z2)*fy1;
    n2 := compress (Y2/Z2-Y1/Z1)*fx2+(X1/Z1-X2/Z2)*fy2;
    d2 := compress (X1/Z1-X2/Z2)*fx2+(Y1/Z1-Y2/Z2)*fy2;
*-
    F1 := (X1/Z1-X2/Z2)^2+(Y1/Z1-Y2/Z2)^2;
    F2 := (y1+a)/(a*y1-1);--n1 / d1;
    F3 := (y2+a)/(a*y2-1); --n2 / d2;
    gateSystem(
	parameters domGS, 
	mapVars, 
	transpose gateMatrix{{F1,F2^2,F3^2}})
    )

jointEuclideanSigMap = dom -> (
    domGS := gateSystem dom;
    pts := reshape vars domGS;
    gateSystem(
    	parameters domGS,
    	vars domGS, 
    	transpose gateMatrix{
	    apply(subsets(4,2), S -> (
		    m := transpose pts^S;
		    (m_(0,0)/m_(2,0)-m_(0,1)/m_(2,1))^2+(m_(1,0)/m_(2,0)-m_(1,1)/m_(2,1))^2
		    )
	    	)
	    }
    )  
)

projectiveJointSigMap = dom -> (
    domGS := gateSystem dom;
    mapVars := vars domGS;
    pts := reshape mapVars;
    M := (i,j,k) -> pts^{i-1,j-1,k-1}_{2,0,1};
    groupMats := {
	{(1,3,5),(1,3,4),(1,2,5)},
	{(1,3,6),(1,3,4),(1,2,6)},
	{(1,3,7),(1,3,4),(1,2,7)},
	{(1,3,8),(1,3,4),(1,2,8)},
	{(1,3,9),(1,3,4),(1,2,9)},	  
	{(2,3,5),(2,3,4),(1,2,5)},
	{(2,3,6),(2,3,4),(1,2,6)},
	{(2,3,7),(2,3,4),(1,2,7)},
	{(2,3,8),(2,3,4),(1,2,8)},
	{(2,3,9),(2,3,4),(1,2,9)}
		    };
    gateSystem(
	parameters domGS,
	mapVars, 
	transpose gateMatrix{
	    apply(groupMats,
		L-> (
		    det3(M(1,2,4)*M(L#0))/
		    det3(M(L#1)*M(L#2))
		    )^2
	    	)
	    }
	)
    )    
