#!/usr/bin/env wolframscript
(* ::Package:: *)

(* This file aims to calculate the frequencies with minimal interaction *)

wp = 25;(*Working precision.*)
ag = 20;(*Accuracy goal.*)
pg = 20;(*Precision goal*)

(* Constrain black brane charge: *)
$Assumptions = {qt < Sqrt[2], qt > 0};

(* Equation *)
myPhiEq:=-(((f[u]-u (f')[u]) (phi')[u])/(u f[u]))+(phi[u] (\[Omega]^2-k^2 f[u]))/(4 rH^2 u f[u]^2)+(phi'')[u];

f[u_]:= 1 - (1+qt^2) u^2+qt^2 u^3
myphiEqRS = Collect[FullSimplify[myPhiEq /. {\[Omega]->rH w,k->rH kk}],{phi''[u],phi'[u],phi[u]}];

(* Define horizon expansion*)
order=14;
phiH[u_]:=(1-u)^((I w)/(2 (-2+qt^2))) Sum[cphiH[nn] (1-u)^nn,{nn,0,order}] /. {\[Alpha]->(I w)/(2 (-2+qt^2))}
(* !!! Change !!!*)
SetDirectory["/home/finn/Documents/studienprojekt/notebooks/"]; (*/home/finn/Do[k,c]ument[e,s]/studienprojekt/notebooks/;/*)
ruleH=Get["ruleHk14"];
phiHExpansion[u_,c0_]:=SetPrecision[phiH[u]/.ruleH/.{cphiH[0]->c0,cphiH[order]->0}, 25] /. {q->qt}

(*phiHExpansion[1/2,1]/.{qt->1/10,w->1/10,kk->0}
Precision[phiHExpansion[1/2,1]/.{qt->1/10,w->1/10,kk->0}]*)


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
phiSol[w_?NumericQ, kk_]:=SetPrecision[solphi[w,kk,1][[1,1,2]],wp]

(* Test *)
Block[
	{qt = 0},
	phiSol[1, 1] // Print
]


(* Get the quasinormal modes*)
(* Usage: qnmRoutine(guess_w_real, guess_w_imaginary, k, q) *)
qnmRoutine[initWR_,initWI_,kay_,kuh_]:=Block[
	{iWR=initWR,iWI=initWI,kk=kay,qt=kuh},
	
	eqQNMVBC[wR_?NumericQ,wI_?NumericQ]:=SetPrecision[Block[{w=wR+I wI},phiSol[w,kk]/.u->uBnum],wp];
	findQNMVBC[reWi_?NumericQ,imWi_?NumericQ]:=FindRoot[{Re[eqQNMVBC[wr,wi]]==0,Im[eqQNMVBC[wr,wi]]==0},{wr,reWi},{wi,imWi},WorkingPrecision->wp,
		AccuracyGoal->ag,
		PrecisionGoal->pg];
	qnmVBC=SetPrecision[findQNMVBC[iWR,iWI], wp];
	Print[{qnmVBC[[1,2]],qnmVBC[[2,2]]}]
]
qnmRoutine[312/100,-275/100,0,0]


(* Thermodynamics *)
Tt[qt_,m_]:=1/(2 \[Pi]) (m/(1+qt^2))^(1/4) Abs[-2+qt^2]
\[Mu]t[qt_,m_]:=Sqrt[3]/2 (m/(1+qt^2))^(1/4) qt
uHt[qt_,m_]:=((1+qt^2)/m)^(1/2)

(* Translating qt to q*)
qPure[qt_,m_]:=qt (m/(1+qt^2))^(3/4)
Plot[qPure[qt,1],{qt,0,Sqrt[2]},Frame->True,FrameLabel->{"q_tilde","q"}]
qPureExtremal=Sqrt[2] (m/3)^(3/4)/.{m->1}//N
qPureExtremal=Rationalize[qPureExtremal, 0]


Print[SetPrecision[Solve[9/10 qPureExtremal==qPure[qt,1], qt][[1,1,2]],wp]];
Print[SetPrecision[Solve[39/40 qPureExtremal==qPure[qt,1], qt][[1,1,2]],wp]];
Print[SetPrecision[Solve[99/100 qPureExtremal==qPure[qt,1], qt][[1,1,2]],wp]];
Print[SetPrecision[Solve[qPureExtremal==qPure[qt,1], qt][[1,1,2]],wp]];


(* Do what we acutally want: Calculate w for differnt k with q at critical charge *)

(* Set closer to critical temp, if desired *)
kkGlobal=0;
qnmRoutine[2,-2,0,0]

tabQNMTemp={};
qtTemp=SetPrecision[Solve[0 qPureExtremal==qPure[qt,1],qt][[1,1,2]],wp]
qtTemp // Print;
 
For[ss=0,ss<201,ss++,
	kkGlobal=ss/10;
	Print["kkGlobal", kkGlobal];
	qnmRoutine[qnmVBC[[1,2]],qnmVBC[[2,2]], qtTemp];
	Print[" wr=", qnmVBC[[1,2]], " wi=", qnmVBC[[2,2]], ", qnmVBC", qnmVBC];
	tabQNMTemp=Append[tabQNMTemp,{qnmVBC[[1,2]],qnmVBC[[2,2]]}/Sqrt[uHt[qtTemp,1]]];
]
Clear[qtTemp];


tabQNMTemp
ListPlot[tabQNMTemp,Frame->True,FrameLabel->{"Re\[Omega]","Im\[Omega]"},PlotLabel->"q=0"]



