(* Find near horizon expansion*)

coeffphiH[nnn_] := 
 Coefficient[
   Normal[
    Series[
     Simplify[(1 - u)^(-\[Alpha] /. indexRule) myphiEqRS /. {phi -> 
         phiH}], {u, 1, order}]], -1 + u, nnn] // Simplify
   
   ruleH = {};

For[iii = 1, iii < order, iii++,
 ruleH = Append[ruleH, 
   Solve[0 == coeffphiH[-2 + iii] /. ruleH, cphiH[iii]][[1, 1]] // 
    Simplify]
 ]