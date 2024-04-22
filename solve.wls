#!/usr/bin/env wolframscript
(* ::Package:: *)

(* This file aims to calculate the frequencies with minimal interaction *)

wp = 8;(*Working precision.*)
ag = 5;(*Accuracy goal.*)
pg = 5;(*Precision goal*)

(* Constrain black brane charge: *)
$Assumptions = {qt < Sqrt[2], qt > 0};

(* Equation *)
myPhiEq:=-(((f[u]-u (f')[u]) (phi')[u])/(u f[u]))+(phi[u] (\[Omega]^2-k^2 f[u]))/(4 rH^2 u f[u]^2)+(phi'')[u];

f[u_]:= 1 - (1+qt^2) u^2+qt^2 u^3
myphiEqRS = Collect[FullSimplify[myPhiEq /. {\[Omega]->rH w,k->rH kk}],{phi''[u],phi'[u],phi[u]}];
Print[myphiEqRS/.{u->t,w->omega}]
(* Define horizon expansion*)
order=14;
phiH[u_]:=(1-u)^((I w)/(2 (-2+qt^2))) Sum[cphiH[nn] (1-u)^nn,{nn,0,order}] /. {\[Alpha]->(I w)/(2 (-2+qt^2))}
(* !!! Change !!!*)
SetDirectory["/home/finn/Documents/studienprojekt/notebooks/"]; (*/home/finn/Do[k,c]ument[e,s]/studienprojekt/notebooks/;/*)
ruleH=Get["ruleHk14"];
phiHExpansion[u_,c0_]:=SetPrecision[phiH[u]/.ruleH/.{cphiH[0]->c0,cphiH[order]->0}, 25] /. {q->qt}


(*phiH[u]/.ruleH/.{cphiH[0]->c0,cphiH[order]->0}/.{q->qt}
D[phiH[u]/.ruleH/.{cphiH[0]->c0,cphiH[order]->0}/.{q->qt},u]*)

phiHExpansion[1/2,1]/.{qt->1/10,w->1/10,kk->0}
SetPrecision[D[phiH[u]/.ruleH/.{cphiH[0]->c0,cphiH[order]->0}/.{q->qt},u], 25]/.{u->1/2,c0->1,qt->1/10,w->1/10,kk->0}
Precision[phiHExpansion[1/2,1]/.{qt->1/10,w->1/10,kk->0}]


(* Solve differential equations*)
epsilon=1/1000000000;
uBnum=epsilon;(*Numerical value of the AdS boundary.*)
uHnum=1-1/10;(*Numerical value of the horizon.*)

Clear[solphi]
solphi[omega_?NumericQ,kay_,c0_]:=
Block[
{w=omega,kk=kay},
	NDSolve[
		{
			(*The three equations*)
			0==myphiEqRS,
			phi[uHnum]==phiHExpansion[uHnum,c0],
			phi'[uHnum]==(D[phiHExpansion[u,c0],u]/.{u->uHnum})
		},
		phi[u],{u,uBnum,uHnum},
		WorkingPrecision->wp,
		AccuracyGoal->ag,
		PrecisionGoal->pg,
		MaxSteps->100000
	]
]
(* Usage: phiSol(w, k) *)
phiSol[w_?NumericQ, kk_]:=Quiet[SetPrecision[solphi[w,kk,1][[1,1,2]],wp]]
(*[[1,1,2]]*)
(* Test *)
Block[{qt = 0},
	sol = phiSol[1/2, 1/2];
	sol // Print;
	sol /. u->(uBnum) // Print;
	sol /. u->(1/2) // Print;
	ReImPlot[phi[u] /.{phi[u] -> sol}, {u, uBnum, uHnum}]]



(* Get the quasinormal modes*)
(* Usage: qnmRoutine(guess_w_real, guess_w_imaginary, k, q) *)
qnmRoutine[initWR_,initWI_,kay_,kuh_]:=Block[
	{iWR=initWR,iWI=initWI,kk=kay,qt=kuh},
	
	eqQNMVBC[wR_?NumericQ,wI_?NumericQ]:=SetPrecision[Block[{w=wR+I wI},phiSol[w,kk]/.u->uBnum],wp];
	findQNMVBC[reWi_?NumericQ,imWi_?NumericQ]:=
		FindRoot[
			{Re[eqQNMVBC[wr,wi]]==0, Im[eqQNMVBC[wr,wi]]==0},
			{wr,reWi},
			{wi,imWi},
			WorkingPrecision->wp,
		AccuracyGoal->ag,
		PrecisionGoal->pg];
	qnmVBC=SetPrecision[findQNMVBC[iWR,iWI], wp];
	Print[{qnmVBC[[1,2]],qnmVBC[[2,2]]}]
]

qnmRoutine[1,1,1,0]


(* Thermodynamics *)
Tt[qt_,m_]:=1/(2 \[Pi]) (m/(1+qt^2))^(1/4) Abs[-2+qt^2]
\[Mu]t[qt_,m_]:=Sqrt[3]/2 (m/(1+qt^2))^(1/4) qt
uHt[qt_,m_]:=((1+qt^2)/m)^(1/2)

(* Translating qt to q*)
qPure[qt_,m_]:=qt (m/(1+qt^2))^(3/4)
Plot[qPure[qt,1],{qt,0,Sqrt[2]},Frame->True,FrameLabel->{"\!\(\*OverscriptBox[\(q\), \(~\)]\)","q"}]
qPureExtremal=Sqrt[2] (m/3)^(3/4)/.{m->1}//N
qPureExtremal=Rationalize[qPureExtremal, 0]


(* Replicate paper plot Example 2*)

(* Calculate values for the initial guess *)
(*qnmRoutine[3,-3,0,0] 
listOfQNM = {};

For[ss=1,ss<40,ss++,
	qtTemp=SetPrecision[Solve[ss/40 qPureExtremal==qPure[qt,1], qt][[1,1,2]],wp];
	Print["ss: ", ss];
	(*Print[qtTemp];
	Print["Tt: ", Tt[qtTemp,1]];
	Print["\[Mu]t: ", \[Mu]t[qtTemp,1]];
	Print["uHT: ", uHt[qtTemp,1]];*)

	qnmRoutine[qnmVBC[[1,2]], qnmVBC[[2,2]], 0, qtTemp];
	(*Print[" wr=", qnmVBC[[1,2]], " wi=", qnmVBC[[2,2]], ", qnmVBC", qnmVBC];*)
	listOfQNM=Append[listOfQNM,{qnmVBC[[1,2]],qnmVBC[[2,2]]}*2/(2-qtTemp^2)];
	Print[{qnmVBC[[1,2]],qnmVBC[[2,2]]}*2/(2-qtTemp^2)];
]

ListPlot[listOfQNM,Frame->True,FrameLabel->{"Re \[Omega]/(\[Pi]T)","Im \[Omega]/(\[Pi]T)"},PlotLabel->"k=0"]*)


Print[SetPrecision[Solve[9/10 qPureExtremal==qPure[qt,1], qt][[1,1,2]],wp]];
Print[SetPrecision[Solve[39/40 qPureExtremal==qPure[qt,1], qt][[1,1,2]],wp]];
Print[SetPrecision[Solve[99/100 qPureExtremal==qPure[qt,1], qt][[1,1,2]],wp]];
Print[SetPrecision[Solve[qPureExtremal==qPure[qt,1], qt][[1,1,2]],wp]];


(* Example 6 *)

(* Set closer to critical temp, if desired *)
(*SetPrecision[Solve[39/40 qPureExtremal\[Equal]qPure[qt,1], qt][[1,1,2]],wp];*)
(*qtTemp=SetPrecision[Solve[39/40 qPureExtremal==qPure[qt,1], qt][[1,1,2]],wp];
qtTemp=0

(* Calculate values for the initial guess *)
qnmRoutine[2,-2, 0, qtTemp];
listOfQNMk = {};

calc[kTemp_] := (
	Print["kTemp: ", kTemp];
	qnmRoutine[qnmVBC[[1,2]], qnmVBC[[2,2]], kTemp, qtTemp];
	Print[qnmVBC[[2,2]], ", qnmVBC", qnmVBC];
	Print[" Scaled: ", {kTemp, qnmVBC[[2,2]]/Sqrt[uHt[qtTemp,1]]}];
	Return[{kTemp, qnmVBC[[2,2]]}];
)
listOfQNMk = Table[calc[i], {i, 0, 20, 1}];
Print[listOfQNMk];
ListPlot[listOfQNMk, Frame->True, FrameLabel->{"kTemp", "Im \[Omega]/(\[Pi]T)"}, PlotLabel->"q=0"]*)


(* Do what we acutally want: Calculate w for differnt q at critical charge *)

(* Set closer to critical temp, if desired *)
(*SetPrecision[Solve[39/40 qPureExtremal\[Equal]qPure[qt,1], qt][[1,1,2]],wp];*)
qtTemp=1/10;

(* Calculate values for the initial guess *)
qnmRoutine[2,-2, 0, qtTemp];
listOfQNMk = {};

Clear[calc];
calc[kTemp_, qtt_] := (
	(*Print["kTemp: ", kTemp];*)
	qnmRoutine[qnmVBC[[1,2]], qnmVBC[[2,2]], kTemp, qtt];
	(*Print[" Scaled: ", {kTemp, qnmVBC[[2,2]]/Sqrt[uHt[qtt,1]]}];*)
	Return[{qnmVBC[[1,2]], qnmVBC[[2,2]]}];
)
For[i=1, i < 10, i++,
	Print["i=", i];
	(*qttt = SetPrecision[Solve[10i/100 qPureExtremal==qPure[qt,1], qt][[1,1,2]],wp];*)
	qttt = i/10 * Sqrt[2]
	Print["qttt=", qttt];
	(* Reset initial guess *)
	qnmRoutine[2,-2, 0, qttt];
	
	listOfQNMk = Table[calc[i, qttt], {i, 0, 20, 1}];
	Print[listOfQNMk];
	Print[ListPlot[listOfQNMk, Frame->True, FrameLabel->{"Re^", "Im \[Omega]/(\[Pi]T)"}, PlotLabel->qttt]]
]


(* Critical scaling for (w~q)^(exp.) *)
(* For this we calculate 
	w(k=0, q=Subscript[q, ex]+\[Epsilon])
	for this \[Epsilon] is in range -1% of Subscript[q, ex] to 0*)
(*listOfQNM = {};
Clear[calc];
calc[kTemp_, qtt_] := (
	(*Print["kTemp: ", kTemp];*)
	qnmRoutine[qnmVBC[[1,2]], qnmVBC[[2,2]], kTemp, qtt];
	(*Print[" Scaled: ", {kTemp, qnmVBC[[2,2]]/Sqrt[uHt[qtt,1]]}];*)
	Return[{kTemp, qnmVBC[[2,2]]}];
)
For[i=1, i < 10, i++,
	Print["i=", i];
	qttt = (i/10) Sqrt[2];
	Print["qttt=", N[qttt]];
	(* Reset initial guess *)
	qnmRoutine[1,-1, 0, qttt];
	qnmRoutine[2,-2, 0, qttt];
	qnmRoutine[3,-3, 0, qttt];
	listOfQNM=Append[listOfQNM, {qttt, qnmVBC[[2,2]]}];
]
ListPlot[listOfQNM, Frame->True, FrameLabel->{"q", "Im \[Omega]/(\[Pi]T)"}, PlotLabel->"Scaling?"]*)


Print[listOfQNM];
ListPlot[listOfQNM]
