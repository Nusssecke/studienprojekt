#!/usr/bin/env wolframscript
(* ::Package:: *)

wp=25;(*Working precision.*)
ag=20;(*Accuracy goal.*)
pg=20;(*Precision goal*)
$Assumptions={qt<Sqrt[2],qt>0}; (* Constrain black brane charge: *)
myphiEqRS = ((-kk^2 (-1+u) (-1+u (-1+qt^2 u))+w^2) phi[u])/(4 (-1+u)^2 u (1+u-qt^2 u^2)^2)+((-1+u (-1+qt^2 u)) (-1+u^2 (-1+qt^2 (-1+2 u))) Derivative[1][phi][u])/((-1+u) u (1+u-qt^2 u^2)^2)+((-1+u (-1+qt^2 u))^2 (phi^\[Prime]\[Prime])[u])/(1+u-qt^2 u^2)^2;


Quit;
(* Find near horizon expansion for non-extremal charge *)
phiIndH[u_]:=(1-u)^\[Alpha]
indexRule={\[Alpha] -> I*w/(-4+2*qt^2)};

order=3; (* Order of the series until which to expand *)
phiH[u_]:=(1-u)^\[Alpha] Sum[cphiH[nn] (1-u)^nn, {nn, 0, order}]/.indexRule

coeffphiH[nnn_] := 
	Coefficient[
		Normal[
			Series[
			Simplify[(1 - u)^(-\[Alpha] /. indexRule) myphiEqRS /. {phi -> phiH}], {u, 1, order}
			]
		], -1 + u, nnn] // Simplify

ruleH = {};

For[iii = 1, iii < order, iii++,
	Print["Step: ", iii];
	Print[coeffphiH[-2 + iii]];
	(* -1.st coefficient is first non-zero *)
	ruleH = Append[ruleH, 
				Solve[0 == coeffphiH[-2 + iii] /. ruleH, cphiH[iii]][[1, 1]] // Simplify
			]
]
 
(* Save expansion *)
SetDirectory["/home/finn/Documents/studienprojekt/"];
Save["normalExpansionRule", ruleH];

phiHExpansion[u_,c0_]:=SetPrecision[phiH[u]/.ruleH/.{cphiH[0]->c0,cphiH[order]->0},wp]/.{q->qt};
phiHExpansion[1/2,1]/.{qt->1/10,w->1/10,kk->0}






(* Find near horizon expansion for extremal charge *)
phiIndH[u_]:=e^(-I * w / (6 (1 - u)))
phiIndH[u] // Print;


order=3;  (* Order of the series until which to expand *)
phiH[u_]:= phiIndH[u] * Sum[cphiH[nn] (1-u)^nn, {nn, 0, order}]
phiH[u] // Print;


coeffphiH[nnn_] := 
	Coefficient[
		Normal[
			Series[
				Simplify[phiIndH[u] * myphiEqRS /. {phi -> phiH}], {u, 1, order}]
			], -1 + u, nnn] // Simplify


ruleH = {};

For[iii = 1, iii < order, iii++,
	Print["Step: ", iii];
	(* -4.th coefficient is first non-zero *)
	If[coeffphiH[-5 + iii] == 0, Continue[]]; (* Failsafe for if coefficient is 0*)
	ruleH = Append[ruleH, 
				Solve[0 == coeffphiH[-2 + iii] /. ruleH, cphiH[iii]][[1, 1]] // Simplify
			]
 ]
 
(* Save expansion *)
SetDirectory["/home/finn/Documents/studienprojekt/"];
Save["extremalExpansionRule", ruleH];

phiHExpansion[u_,c0_]:=SetPrecision[phiH[u]/.ruleH/.{cphiH[0]->c0,cphiH[order]->0},wp]/.{q->qt}
phiHExpansion[1/2,1]/.{qt->Sqrt[2],w->1/10,kk->0}
