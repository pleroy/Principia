(* ::Package:: *)

(* ::Section:: *)
(*Declarations*)


(* ::Input:: *)
(*<<FunctionApproximations`*)


(* ::Input:: *)
(*ClearAll[CorrectlyRound]*)


(* ::Input:: *)
(*Get[FileNameJoin[{NotebookDirectory[],"floating_point.wl"}]]*)


(* ::Input:: *)
(*SetFloatingPointFormat[binary64]*)
(*SetRoundingMode[NearestTiesToEven]*)


(* ::Subsubsection:: *)
(*Machine Evaluation with Rounding*)


(* ::Input:: *)
(*ClearAll[machineEvaluate];*)
(*SetAttributes[machineEvaluate,HoldAll];*)
(*machineEvaluate[x:(_Plus|_Times|_Power|_?NumberQ)]:=Block[{Plus,Times,me},*)
(*ClearAttributes[Plus,Flat];*)
(*ClearAttributes[Times,Flat];*)
(*SetAttributes[me,HoldAll];*)
(*me[a_*b_+c_]:=CorrectlyRound[machineEvaluate[a]machineEvaluate[b]+machineEvaluate[c]];*)
(*me[a_ +b_]:=CorrectlyRound[machineEvaluate[a]+machineEvaluate[b]];*)
(*me[a_-b_]:=CorrectlyRound[machineEvaluate[a]-machineEvaluate[b]];*)
(*me[a_*b_]:=CorrectlyRound[machineEvaluate[a]*machineEvaluate[b]];*)
(*me[a_/b_]:=CorrectlyRound[machineEvaluate[a]/machineEvaluate[b]];*)
(*me[a_^2]:=CorrectlyRound[machineEvaluate[a ]machineEvaluate[a]];*)
(*me[a_^3]:=CorrectlyRound[machineEvaluate[a^2 ]machineEvaluate[a]];*)
(*me[a_^4]:=CorrectlyRound[machineEvaluate[a^2 ]machineEvaluate[a^2]];*)
(*me[a_?NumberQ]:=CorrectlyRound[a];*)
(*me[x]]*)


(* ::Input:: *)
(*machineEvaluate[(2x+1/3)+y^2]*)


(* ::Subsubsection:: *)
(*Machine Evaluation with Intervals*)


(* ::Input:: *)
(*ClearAll[mei2];*)
(*applyOpToIntervals[op_,{exact_Interval,approx_Interval}]:=Block[{e=op[exact],a=op[approx]},{e,a*Interval[{1-2^-53,1+2^-53}]}];*)
(*applyOpToIntervals[op_,{exact1_Interval,approx1_Interval},{exact2_Interval,approx2_Interval}]:=Block[{e=op[exact1,exact2],a=op[approx1,approx2]},{e,a*Interval[{1-2^-53,1+2^-53}]}];*)
(*applyOpToIntervals[op_,{exact1_Interval,approx1_Interval},{exact2_Interval,approx2_Interval},{exact3_Interval,approx3_Interval}]:=Block[{e=op[exact1,exact2,exact3],a=op[approx1,approx2,approx3]},{e,a*Interval[{1-2^-53,1+2^-53}]}]*)


(* ::Input:: *)
(*ClearAll[machineEvaluateInterval];SetAttributes[machineEvaluateInterval,HoldAll];machineEvaluateInterval[x:(_Plus|_Times|_Power|_Interval|_?NumberQ)]:=Block[{Plus,Times,mei},*)
(*ClearAttributes[Plus,Flat];*)
(*ClearAttributes[Times,Flat];*)
(*SetAttributes[mei,HoldAll];*)
(*mei[a_*b_+c_]:=applyOpToIntervals[#1 #2+#3&,machineEvaluateInterval[a],machineEvaluateInterval[b],machineEvaluateInterval[c]];*)
(*mei[a_ +b_]:=applyOpToIntervals[Plus,machineEvaluateInterval[a],machineEvaluateInterval[b]];*)
(*mei[a_-b_]:=applyOpToIntervals[Subtract,machineEvaluateInterval[a],machineEvaluateInterval[b]];*)
(*mei[a_*b_]:=applyOpToIntervals[Times,machineEvaluateInterval[a],machineEvaluateInterval[b]];*)
(*mei[a_/b_]:=applyOpToIntervals[Divide,machineEvaluateInterval[a],machineEvaluateInterval[b]];*)
(*mei[a_^2]:=applyOpToIntervals[#^2&,machineEvaluateInterval[a ]];*)
(*mei[a_^3]:=applyOpToIntervals[Times,machineEvaluateInterval[a^2 ],machineEvaluateInterval[a]];*)
(*mei[a_^4]:=applyOpToIntervals[#^2&,machineEvaluateInterval[a^2 ]];*)
(*mei[a_^5]:=applyOpToIntervals[Times,machineEvaluateInterval[a^4 ],machineEvaluateInterval[a]];*)
(*mei[a_Interval]:={a,a};*)
(*mei[a_?NumberQ]:=Block[{ca=CorrectlyRound[a],i},i=Interval[{ca,ca}];{i,i}];*)
(*mei[x]]*)


(* ::Input:: *)
(*ClearAll[bounds];*)
(*bounds[{exact_Interval,approx_Interval}]:={Log2[Max[Abs[Min[exact]-Min[approx]],Abs[Max[exact]-Max[approx]]]],Log2[Max[Abs[Min[exact]],Abs[Min[approx]],Abs[Max[exact]],Abs[Max[approx]]]]}*)


(* ::Input:: *)
(*i=machineEvaluateInterval[Interval[{1/3,2/3}]+Interval[{1/7,2/7}]]*)


(* ::Input:: *)
(*Min[i[[2]]]<=10/21*)


(* ::Input:: *)
(*Max[i[[2]]]>=20/21*)


(* ::Section::Closed:: *)
(*Sin*)


(* ::Text:: *)
(*We want to compute minimax polynomials that are odd and where the term in t is linear so that the computation produces two terms, t plus a correction.  Therefore, we approximate the function sinFn below.  Note that it is singular near 0 and must therefore be replaced by its limit there.*)


(* ::Input:: *)
(*Series[Sin[t],{t,0,5}]*)


(* ::Input:: *)
(*sinFn[t_]:=(Sin[t]-t)/t^3*)


(* ::Input:: *)
(*Limit[sinFn[t],t->0]*)


(* ::Input:: *)
(*Series[sinFn[t],{t,0,2}]*)


(* ::Subsection:: *)
(*Minimax Polynomials Between Table Entries*)


(* ::Text:: *)
(*Let's compute degree-1 minimax polynomials for various interval radii (i.e., various table sizes).  The error function is chosen to minimize absolute error:*)


(* ::Input:: *)
(*sinGeneralApproximationResults=Table[*)
(*{2^-n,*)
(*GeneralMiniMaxApproximation[*)
(*{t^2,If[t==0,-1/6,sinFn[t]],If[t==0,1,1/t^3]},*)
(*{t,{0,2^-n},1,0},*)
(*x,WorkingPrecision->30]},*)
(*{n,7,10}]*)


(* ::Input:: *)
(*sinGeneralApproximations=Map[*)
(*{#[[1]],Function[u, Evaluate[ u^3 (#[[2,2,1]]/.x->u^2)]]}&,*)
(*sinGeneralApproximationResults]*)


(* ::Text:: *)
(*The machine evaluation introduces some numeric noise for radius 1/1024:*)


(* ::Input:: *)
(*Map[Plot[Sin[x]-x-machineEvaluate[#[[2]][x]],{x,-#[[1]],#[[1]]},PlotRange->Full,WorkingPrecision->30]&,sinGeneralApproximations]*)


(* ::Text:: *)
(*For radii other than 1/1024, the error estimate returned by the Remez algorithm is pretty good.  (Note that FindMaximum requires tight constraints because the functions are not very smooth.)*)


(* ::Input:: *)
(*sinGeneralApproximationExtrema=Map[{#[[2,2,2]],#[[2,1,2]]}&,sinGeneralApproximationResults]*)


(* ::Input:: *)
(*sinGeneralApproximationsMachineExtrema=MapThread[*)
(*FindMaximum[{Sin[x]-x-machineEvaluate[#1[[2]][x]],x>=99#2[[2]]/100&&x<=#1[[1]]},{x,#2[[2]]},WorkingPrecision->30]&,*)
(*{sinGeneralApproximations,sinGeneralApproximationExtrema}]*)


(* ::Text:: *)
(*This gives us the number of bits available for each radius.  1/128 is not appropriate because it could only use 11 zeroes after the mantissa, but all the other radii work since they can use at least 18 bits after the mantissa:*)


(* ::Input:: *)
(*MapThread[{#1[[1]],Log[2,Abs[#2[[1]]]]}&,{sinGeneralApproximationResults,sinGeneralApproximationMachineExtrema}]*)


(* ::Text:: *)
(*For radius 1/1024, it is useful to do an exhaustive search because the function is noisy.  We see that we lose about 0.5 bits compared to the extrema found by FindMaximum:*)


(* ::Input:: *)
(*N[Log[2,Max[Abs[Table[Sin[x]-x-machineEvaluate[sinGeneralApproximations[[4,2]][x]],{x,9 sinGeneralApproximations[[4,1]]/10,sinGeneralApproximations[[4,1]],sinGeneralApproximations[[4,1]]/1*^6}]]]],30]*)


(* ::Subsection:: *)
(*Minimax  Polynomials  Near  Zero*)


(* ::Text:: *)
(*Near zero the Gal and Bachelis method is not usable, so we use a plain minimax polynomial for the entire function.  We still want it to be odd and we still want to separate out the t term so we use the sinFn function again, but this time with an error function suitable for the relative error.  We can choose either degree 1 or degree 2, with different radii.*)


(* ::Subsubsection:: *)
(*Degree  1*)


(* ::Input:: *)
(*sin0GeneralApproximationResults=Map[{#,GeneralMiniMaxApproximation[*)
(*{t^2,If[t==0,-1/6,sinFn[t]],If[t==0,1,Sin[t]/t^3]},*)
(*{t,{0,#},1,0},*)
(*x,WorkingPrecision->30]}&,{1/512,1/1024}]*)


(* ::Input:: *)
(*sin0GeneralApproximations=Map[{#[[1]],Function[u, Evaluate[ u^3 (#[[2,2,1]]/.x->u^2)]]}&,sin0GeneralApproximationResults]*)


(* ::Text:: *)
(*This time the plot of interest is for the relative error, and it is noisy for radius 1/1024:*)


(* ::Input:: *)
(*Map[Plot[1-(x+machineEvaluate[#[[2]][x]])/Sin[x],{x,-#[[1]],#[[1]]},WorkingPrecision->30,PlotRange->Full]&,sin0GeneralApproximations]*)


(* ::Text:: *)
(*In this case, FindMaximum suspiciously returns an error less than the one computed by the Remez algorithm:*)


(* ::Input:: *)
(*sin0GeneralApproximationExtrema=Map[{#[[2,2,2]],#[[2,1,2]]}&,sin0GeneralApproximationResults]*)


(* ::Input:: *)
(*sin0GeneralApproximationMachineExtrema=MapThread[*)
(*FindMaximum[{1-(x+machineEvaluate[#[[2]][x]])/Sin[x],x>=99#2[[2]]/100&&x<=#1[[1]]},{x,#2[[2]]},WorkingPrecision->30]&,*)
(*{sin0GeneralApproximations,sin0GeneralApproximationExtrema}]*)


(* ::Text:: *)
(*This gives us the number of bits available for each radius.  1/512 is not appropriate as it would only provide 16 bits after the mantissa:*)


(* ::Input:: *)
(*MapThread[{#1[[1]],Log[2,Abs[#2[[1]]]]}&,{sin0GeneralApproximationResults,sin0GeneralApproximationMachineExtrema}]*)


(* ::Text:: *)
(*For radius 1/1024, it is useful to do an exhaustive search because the function is noisy.  We see that we lose about 1 bit compared to the extrema found by FindMaximum, but we are left with 21 bits after the mantissa, which is sufficient:*)


(* ::Input:: *)
(*N[Log[2,Max[Abs[Table[1-(x+machineEvaluate[sin0GeneralApproximations[[2,2]][x]])/Sin[x],{x,9 sin0GeneralApproximations[[2,1]]/10,sin0GeneralApproximations[[2,1]],sin0GeneralApproximations[[2,1]]/1*^6}]]]],30]*)


(* ::Subsubsection:: *)
(*Degree  2*)


(* ::Text:: *)
(*Note: it is not clear if a radius of 1/512 is sufficient from the perspective of error analysis.  If it isn't we might as well go with 1/1024 and degree 1.*)


(* ::Input:: *)
(*sin0GeneralApproximationResult2=Map[{#,GeneralMiniMaxApproximation[*)
(*{t^2,If[t==0,-1/6,sinFn[t]],If[t==0,1,Sin[t]/t^3]},*)
(*{t,{0,#},2,0},*)
(*x,WorkingPrecision->30]}&,{1/512}]*)


(* ::Input:: *)
(*sin0GeneralApproximation2=Map[{#[[1]],Function[u, Evaluate[ u^3 (#[[2,2,1]]/.x->u^2)]]}&,sin0GeneralApproximationResult2]*)


(* ::Text:: *)
(*The result is extremely noisy:*)


(* ::Input:: *)
(*Map[Plot[1-(x+machineEvaluate[#[[2]][x]])/Sin[x],{x,-#[[1]],#[[1]]},WorkingPrecision->30,PlotRange->Full]&,sin0GeneralApproximation2]*)


(* ::Text:: *)
(*Therefore it doesn't make sense to search for the maximum using FindMaximum.  Instead, we just do an exhaustive search.  We find that this gives us less bits than radius 1/1024 and degree 1:*)


(* ::Input:: *)
(*N[Log[2,Max[Abs[Table[1-(x+machineEvaluate[sin0GeneralApproximation2[[1,2]][x]])/Sin[x],{x,9 sin0GeneralApproximation2[[1,1]]/10,sin0GeneralApproximation2[[1,1]],sin0GeneralApproximation2[[1,1]]/1*^6}]]]],30]*)


(* ::Section::Closed:: *)
(*Cos*)


(* ::Text:: *)
(*We want to compute minimax polynomials that are even and return 1 for 0.  Therefore, we approximate the function cosFn1 or cosFn2 below.  Note that they are singular near 0 and must therefore be replaced by their limit there.*)


(* ::Input:: *)
(*Series[Cos[t],{t,0,5}]*)


(* ::Subsection:: *)
(*Minimax  Polynomials Between Table Entries*)


(* ::Text:: *)
(*Between table entries we choose an error function that minimizes the absolute error.  We try degree 1 and degree 2, with various interval radii.*)


(* ::Subsubsection:: *)
(*Degree  1*)


(* ::Input:: *)
(*cosFn1[t_]:=(Cos[t]-1)/t^2*)


(* ::Input:: *)
(*Limit[cosFn1[t],t->0]*)


(* ::Input:: *)
(*Series[cosFn1[t],{t,0,2}]*)


(* ::Input:: *)
(*cosGeneralApproximationResults1=Table[*)
(*{2^-n,*)
(*GeneralMiniMaxApproximation[*)
(*{t^2,If[t==0,-1/2,cosFn1[t]],If[t==0,1,1/t^2]},*)
(*{t,{0,2^-n},1,0},*)
(*x,WorkingPrecision->30]},*)
(*{n,9,12}]*)


(* ::Input:: *)
(*cosGeneralApproximations1=Map[*)
(*{#[[1]],Function[u, Evaluate[ u^2(#[[2,2,1]]/.x->u^2)]]}&,*)
(*cosGeneralApproximationResults1]*)


(* ::Text:: *)
(*For radii other than 1/512, the result is extremely noisy:*)


(* ::Input:: *)
(*Map[Plot[Cos[x]-1-machineEvaluate[#[[2]][x]],{x,-#[[1]],#[[1]]},PlotRange->Full,WorkingPrecision->30]&,cosGeneralApproximations1]*)


(* ::Output:: *)
(*{Graphics[Annotation[{{{{}, {}, Annotation[{Directive[Opacity[1.], RGBColor[0.368417, 0.506779, 0.709798], AbsoluteThickness[2]], Line[CompressedData["*)
(*1:eJwB4QQe+yFib1JlAgAAAE0AAAACAAAALj8W6v//X78UTP0VIXTDu2/5I4j5*)
(*+l+/+JsQtZYHw7uwszEm8/Vfv0eDXc2hfsK7MihNYubrX7/DngMRBUfCuzcR*)
(*hNrM11+/DqoZMttGwbtA4/HKma9fvyDXvsfnUr67UofNqzNfX78J91wqZ162*)
(*u3fPhG1nvl6/INfHS2G9obvU1AvqwWFdv6kDSGOr/6o7HcBvFDccXL+oShTi*)
(*IR66O7VrxAUP3Vq/tdF2zQgjwTtgjDHw2YJZv20lSV20C8M7+JJ7iL8/WL/m*)
(*cW5a3iXDO6MO3hmY4Va/FMNdjtJYwjucSjFy04lVv8WIzgJGbsA7gmxheClJ*)
(*VL9vXTIdfQG9O3sDqndy7VK/hhtPPvZruDthgM8k1qhRv+nS8tYVFrQ7lb3l*)
(*mJxqUL8ksu3CRBqwO7rfKAysIk6/VSv8p2sVqDsjEEBCVJ5LvxsAc2xtoqE7*)
(*siqIauLjSL8luFjMk9yYO9/FsSA2Nka/snKvyMNtkDvlLJUyv7ZDv/4nyWMN*)
(*AoU7EX6pNi4BQb8vUoCz3mh3Oyw27yyl8zy/S6t+ybs2aTuDRO3QuXg3v5/R*)
(*+bds4VU7FVSukFkXMr9ZHIuJoYhBO7L2xQ/IJCq/039PATvrJDul2+WN6YQe*)
(*v6aY2ps4wes6XVEctNtjBL8qhOekcRetOrkKe9XCdQM/nmq3SAL6qTo14gEe*)
(*gVgfP+hz1PyhfvU6I/C7+TrCKT/zzQcqwCMdO0mj2c1OWDI/UaN1Cb+iPztF*)
(*TRID9bU3P0MX7cvHeVY7j1/XgDC3PD8D88qszMJpO8dOHQ1QEkE/e5z+q1i4*)
(*eDvtIRV+0ppDP/5GOgDHLYQ7dXQrYY8WRj+7IW64yoqPO9fcEFJmyEg/vYWB*)
(*XXYkmDtgeTznB0xLP8S09ks6SaE7wys3isMFTj+TjRnqhuynOyYJvOikSFA/*)
(*G+wglPTFrjsdvGtFBYhRP8P+onSZnbM7APoCqXLiUj9bUxstxOe3O/dRvV7F*)
(*JVQ/HWSU4RA9vDvbNF8bJYRVP4o3mRRce8A7cFcQESLcVj+m/MCd6wbCOxmU*)
(*5FgEHVg/4rwIxzw4wzuuW6Cn83hZP49AutpbQMM7Vz1/SMi9Wj9ugitgGp7B*)
(*O7FebSI6/Fs//37NgsRGuzv4CkMDuVVdP8lTiKG1yKo7UtE7Nh2YXj+x8Yoe*)
(*7Z2Zuwo7C8G8nV4/GQtV7PONmLvBpNpLXKNePxcGmXwonJu7MHh5YZuuXj+m*)
(*woU0ayKguw4ft4wZxV4/W1sY7lRlo7vJbDLjFfJePxvOTKjIWKm7QAgpkA5M*)
(*Xz8kH9EC7Jq0u/hx+BquUV8/YJ/+0dmstbuv28elTVdfPzI8LqIacbW7Hq9m*)
(*u4xiXz+rliUJE9a2u/xVpOYKeV8/OZCttD36uLu3ox89B6ZfP6Dm9dajxr27*)
(*bg3vx6arXz+VZNluU4G9uyZ3vlJGsV8/ql3FNWbdvruUSl1ohbxfPwbhMIYf*)
(*ir+7cvGakwPTXz/AmB6bj7bAuypbah6j2F8/p8mz4Z8twbvhxDmpQt5fP70Q*)
(*8z6VPcG7UJjYvoHpXz+c3fC0U+DBuwgCqEkh718/U6alEx1Nwru/a3fUwPRf*)
(*P9wJajJkhcK7dtVGX2D6Xz+H1b5pZvDCuy4/Fur//18/FEz9FSF0w7uG1FnU*)
(**)
(*"]]}, "Charting`Private`Tag#1"]}}, {}}, <|"HighlightElements" -> <|"Label" -> {"XYLabel"}, "Ball" -> {"InterpolatedBall"}|>, "LayoutOptions" -> <|"PanelPlotLayout" -> <||>, "PlotRange" -> {{Rational[-1, 512], Rational[1, 512]}, {-8.238932713621388*^-21, 8.153285712498936*^-21}}, "Frame" -> {{False, False}, {False, False}}, "AxesOrigin" -> {0, 2.2203545531351673`*^-16}, "ImageSize" -> {360, 360/GoldenRatio}, "Axes" -> {True, True}, "LabelStyle" -> {}, "AspectRatio" -> GoldenRatio^(-1), "DefaultStyle" -> {Directive[Opacity[1.], RGBColor[0.368417, 0.506779, 0.709798], AbsoluteThickness[2]]}, "HighlightLabelingFunctions" -> <|"CoordinatesToolOptions" -> ({Identity[Part[#, 1]], Identity[Part[#, 2]]}& ), "ScalingFunctions" -> {{Identity, Identity}, {Identity, Identity}}|>, "Primitives" -> {}, "GCFlag" -> False|>, "Meta" -> <|"DefaultHighlight" -> {"Dynamic", None}, "Index" -> {}, "Function" -> Plot, "GroupHighlight" -> False|>|>, "DynamicHighlight"], AspectRatio -> GoldenRatio^(-1), Axes -> {True, True}, AxesLabel -> {None, None}, AxesOrigin -> {0, 2.2203545531351673`*^-16}, DisplayFunction -> Identity, Frame -> {{False, False}, {False, False}}, FrameLabel -> {{None, None}, {None, None}}, FrameTicks -> {{Automatic, Automatic}, {Automatic, Automatic}}, GridLines -> {None, None}, GridLinesStyle -> Directive[GrayLevel[0.5, 0.4]], ImagePadding -> All, Method -> {"DefaultBoundaryStyle" -> Automatic, "DefaultGraphicsInteraction" -> {"Version" -> 1.2, "TrackMousePosition" -> {True, False}, "Effects" -> {"Highlight" -> {"ratio" -> 2}, "HighlightPoint" -> {"ratio" -> 2}, "Droplines" -> {"freeformCursorMode" -> True, "placement" -> {"x" -> "All", "y" -> "None"}}}}, "DefaultMeshStyle" -> AbsolutePointSize[6], "ScalingFunctions" -> None, "CoordinatesToolOptions" -> {"DisplayFunction" -> ({(Identity[#]& )[Part[#, 1]], (Identity[#]& )[Part[#, 2]]}& ), "CopiedValueFunction" -> ({(Identity[#]& )[Part[#, 1]], (Identity[#]& )[Part[#, 2]]}& )}}, PlotRange -> {{Rational[-1, 512], Rational[1, 512]}, {-8.238932713621388*^-21, 8.153285712498936*^-21}}, PlotRangeClipping -> True, PlotRangePadding -> {{Scaled[0.02], Scaled[0.02]}, {Scaled[0.05], Scaled[0.05]}}, Ticks -> {Automatic, Automatic}],Graphics[Annotation[{{{{}, {}, Annotation[{Directive[Opacity[1.], RGBColor[0.368417, 0.506779, 0.709798], AbsoluteThickness[2]], Line[CompressedData["*)
(*1:eJwB4QQe+yFib1JlAgAAAE0AAAACAAAALj8W6v//T7/jsyQtg25fu2/5I4j5*)
(*+k+/QujRNJEVXruwszEm8/VPv5BgMgJ7QGW7MihNYubrT7/hneCOmBBluzcR*)
(*hNrM10+/oLdSNSxMYLtA4/HKma9Pv+1wmOAE81u7UofNqzNfT78vyM5q7NFZ*)
(*u3fPhG1nvk6/nwI2aVCCTrvU1AvqwWFNv2wKcYYWFkA7HcBvFDccTL8qmL3S*)
(*talYO7VrxAUP3Uq/K1VomeKfYTtgjDHw2YJJvwHOU2BQsGU7+JJ7iL8/SL8C*)
(*iyFKZTViO6MO3hmY4Ua/rx2IpuXvYzucSjFy04lFv2MtqsJbxWE7gmxheClJ*)
(*RL/cGA0sbgVaO3sDqndy7UK/uPERVcRUWTthgM8k1qhBv3MRPA1dllA7lb3l*)
(*mJxqQL/RhDC+LZhMO7rfKAysIj6/AuxD6M0JSzsjEEBCVJ47v9LWm2TUOkE7*)
(*siqIauLjOL8RnrcnlYY9O9/FsSA2Nja/1EM8mXxrMTvlLJUyv7Yzv8oWtMtY*)
(*FiM7EX6pNi4BMb/SOxwNdW4NOyw27yyl8yy/OdqXNfZXFDuDRO3QuXgnvzkY*)
(*b1Wa5gU7FVSukFkXIr9YznZXD2XzOrL2xQ/IJBq/Lt7x34NT1Lql2+WN6YQO*)
(*vwFIdC2NsmG6XVEctNtj9L4XjBEU/6OLurkKe9XCdfM+AWQ8R4GBibo14gEe*)
(*gVgPPw1Ukjn1ubG6I/C7+TrCGT8LhK0+ynfgOkmj2c1OWCI/HfBHyE5I9jpF*)
(*TRID9bUnPyZqPfOCfdE6j1/XgDC3LD/ZiZJxTN0QO8dOHQ1QEjE/TOPPgisu*)
(*ITvtIRV+0pozP7cicEhWOiU7dXQrYY8WNj/2ZraGISczO9fcEFJmyDg/nIGL*)
(*EBjHOztgeTznB0w7P1cvAeiOZz07wys3isMFPj/QRaygduJKOyYJvOikSEA/*)
(*RZCaL70XUjsdvGtFBYhBP2O0/RVTd1M7APoCqXLiQj+AD7w1ampUO/dRvV7F*)
(*JUQ/vuxP1ECSWzvbNF8bJYRFP1zudAMQ6F07cFcQESLcRj87eqLCHyZiOxmU*)
(*5FgEHUg/eQinU0LSZjuuW6Cn83hJP2B4Lf5O/GM7Vz1/SMi9Sj9LF2bmx1Bi*)
(*O7FebSI6/Es/X/2B2n9kYTv4CkMDuVVNPwj8JP1nuEw7UtE7Nh2YTj8LU7TN*)
(*Sds7uwo7C8G8nU4/iJn0UB3IHbvBpNpLXKNOP1B9JY/18UO7MHh5YZuuTj/o*)
(*j/+cxYQ/uw4ft4wZxU4/ucAuYDYxS7vJbDLjFfJOPybSrd4T6ka7QAgpkA5M*)
(*Tz+jeQeq8dZUu/hx+BquUU8/VucEQGxPTruv28elTVdPP1mJTc8c01S7Hq9m*)
(*u4xiTz95UaTJx0RTu/xVpOYKeU8//zHjjrs/Vru3ox89B6ZPP6luqsrqQVm7*)
(*bg3vx6arTz9A7u+Ez3pfuyZ3vlJGsU8/ndk3eWwIWLuUSl1ohbxPP+ddr1qK*)
(*yWC7cvGakwPTTz++SzDEOAVfuypbah6j2E8/wY5oaTYmY7vhxDmpQt5PP2km*)
(*3u2nhl67UJjYvoHpTz/HIQhoNR9fuwgCqEkh708/LiSX4/kqY7u/a3fUwPRP*)
(*P8PB3m7KAWa7dtVGX2D6Tz91GaX9kKhjuy4/Fur//08/47MkLYNuX7vDLjqD*)
(**)
(*"]]}, "Charting`Private`Tag#1"]}}, {}}, <|"HighlightElements" -> <|"Label" -> {"XYLabel"}, "Ball" -> {"InterpolatedBall"}|>, "LayoutOptions" -> <|"PanelPlotLayout" -> <||>, "PlotRange" -> {{Rational[-1, 1024], Rational[1, 1024]}, {-1.4563007772160165`*^-22, 1.5101885787122222`*^-22}}, "Frame" -> {{False, False}, {False, False}}, "AxesOrigin" -> {0, 2.2204444281445765`*^-16}, "ImageSize" -> {360, 360/GoldenRatio}, "Axes" -> {True, True}, "LabelStyle" -> {}, "AspectRatio" -> GoldenRatio^(-1), "DefaultStyle" -> {Directive[Opacity[1.], RGBColor[0.368417, 0.506779, 0.709798], AbsoluteThickness[2]]}, "HighlightLabelingFunctions" -> <|"CoordinatesToolOptions" -> ({Identity[Part[#, 1]], Identity[Part[#, 2]]}& ), "ScalingFunctions" -> {{Identity, Identity}, {Identity, Identity}}|>, "Primitives" -> {}, "GCFlag" -> False|>, "Meta" -> <|"DefaultHighlight" -> {"Dynamic", None}, "Index" -> {}, "Function" -> Plot, "GroupHighlight" -> False|>|>, "DynamicHighlight"], AspectRatio -> GoldenRatio^(-1), Axes -> {True, True}, AxesLabel -> {None, None}, AxesOrigin -> {0, 2.2204444281445765`*^-16}, DisplayFunction -> Identity, Frame -> {{False, False}, {False, False}}, FrameLabel -> {{None, None}, {None, None}}, FrameTicks -> {{Automatic, Automatic}, {Automatic, Automatic}}, GridLines -> {None, None}, GridLinesStyle -> Directive[GrayLevel[0.5, 0.4]], ImagePadding -> All, Method -> {"DefaultBoundaryStyle" -> Automatic, "DefaultGraphicsInteraction" -> {"Version" -> 1.2, "TrackMousePosition" -> {True, False}, "Effects" -> {"Highlight" -> {"ratio" -> 2}, "HighlightPoint" -> {"ratio" -> 2}, "Droplines" -> {"freeformCursorMode" -> True, "placement" -> {"x" -> "All", "y" -> "None"}}}}, "DefaultMeshStyle" -> AbsolutePointSize[6], "ScalingFunctions" -> None, "CoordinatesToolOptions" -> {"DisplayFunction" -> ({(Identity[#]& )[Part[#, 1]], (Identity[#]& )[Part[#, 2]]}& ), "CopiedValueFunction" -> ({(Identity[#]& )[Part[#, 1]], (Identity[#]& )[Part[#, 2]]}& )}}, PlotRange -> {{Rational[-1, 1024], Rational[1, 1024]}, {-1.4563007772160165`*^-22, 1.5101885787122222`*^-22}}, PlotRangeClipping -> True, PlotRangePadding -> {{Scaled[0.02], Scaled[0.02]}, {Scaled[0.05], Scaled[0.05]}}, Ticks -> {Automatic, Automatic}],Graphics[Annotation[{{{{}, {}, Annotation[{Directive[Opacity[1.], RGBColor[0.368417, 0.506779, 0.709798], AbsoluteThickness[2]], Line[CompressedData["*)
(*1:eJwB4QQe+yFib1JlAgAAAE0AAAACAAAALj8W6v//P7/D7Tazyb7VOm/5I4j5*)
(*+j+/S+q1NCwpAzuwszEm8/U/v8ICzNvIUCO7MihNYubrP79uN8qUOykSOzcR*)
(*hNrM1z+/Ove30qzo57pA4/HKma8/v1yoZL51aga7UofNqzNfP79UbSBfqbsP*)
(*u3fPhG1nvj6/VzjOLRGPHDvU1AvqwWE9v4+d/bRYwQM7HcBvFDccPL+j7nYH*)
(*WKQOO7VrxAUP3Tq/4e8b92seITtgjDHw2YI5vw0GsB/bJRS7+JJ7iL8/OL8z*)
(*gmCgfzkWu6MO3hmY4Ta/2XKsHjTUGDucSjFy04k1vwJq8JkkKs26gmxheClJ*)
(*NL89vohTrS0QO3sDqndy7TK/xv2uUa6CEDthgM8k1qgxv4azIcLTJeG6lb3l*)
(*mJxqML+lGuH9pmjpOrrfKAysIi6/ZqutyOhz27ojEEBCVJ4rv8BC4lZz0fG6*)
(*siqIauLjKL/i2E75ANnyut/FsSA2Nia/b+RTa6Gb8zrlLJUyv7Yjv60+acOD*)
(*UeM6EX6pNi4BIb9slywu1fvtOiw27yyl8xy/wVz9G3XbqTqDRO3QuXgXvxnp*)
(*tygJQ9O6FVSukFkXEr/Dq61iW7HGurL2xQ/IJAq/UAmNkQR1vDql2+WN6YT+*)
(*vhJ+jrdfYpE6XVEctNtj5L4DsjR5Ai1mOrkKe9XCdeM+pwn03gc/Uzo14gEe*)
(*gVj/PuenCvnISpQ6I/C7+TrCCT+6dEDgrfeQukmj2c1OWBI/Ft/7MQpnsbpF*)
(*TRID9bUXP9XuRfk/6tg6j1/XgDC3HD/5+llOUeHdOsdOHQ1QEiE/rRvgu65Z*)
(*wTrtIRV+0pojP8NZYE+Ey+O6dXQrYY8WJj8zDBZLcrWQutfcEFJmyCg/YBbZ*)
(*RHuvADtgeTznB0wrP8K5M7MGngM7wys3isMFLj/HnR7+AA2mOiYJvOikSDA/*)
(*55aFb4C85rodvGtFBYgxP+Z/Pr1UFf26APoCqXLiMj/ADNTsX7f1uvdRvV7F*)
(*JTQ/niNGMAmR+zrbNF8bJYQ1P88xiILIgPG6cFcQESLcNj9AE/Rs2PXwOhmU*)
(*5FgEHTg/6wLgBdNHHzuuW6Cn83g5P4vlgE6+OxK7Vz1/SMi9Oj9IxT5hMXEQ*)
(*O7FebSI6/Ds/IzDuI3+EIDv4CkMDuVU9P+Hno/vtUiE7UtE7Nh2YPj9Qz7Wj*)
(*xTTyugo7C8G8nT4/k9Y7B/mSCbvBpNpLXKM+P7NjAEMoqQI7MHh5YZuuPj94*)
(*Hujl1z/zug4ft4wZxT4/eCkDaW7BDLvJbDLjFfI+P04UGfDKmxS7QAgpkA5M*)
(*Pz+bZoUDGoYQO/hx+BquUT8/nHSAuUbjDTuv28elTVc/P2Pocf9nUBC7Hq9m*)
(*u4xiPz9ymFP3TYvHuvxVpOYKeT8/TKLTCGb0FDu3ox89B6Y/P6hwZ4RVvRQ7*)
(*bg3vx6arPz8d1R28Hk4ROyZ3vlJGsT8/9Dmc0rbAE7uUSl1ohbw/P8dz7GHg*)
(*Ehc7cvGakwPTPz+mC0fmTlMOuypbah6j2D8/dH0s6fFw9LrhxDmpQt4/Pycz*)
(*wPjbvxy7UJjYvoHpPz+zv3ZSvhkJuwgCqEkh7z8/0LqhhyqfCDu/a3fUwPQ/*)
(*P6P9wy9fxRm7dtVGX2D6Pz8bfJkvpdgiuy4/Fur//z8/w+02s8m+1TqjLE0T*)
(**)
(*"]]}, "Charting`Private`Tag#1"]}}, {}}, <|"HighlightElements" -> <|"Label" -> {"XYLabel"}, "Ball" -> {"InterpolatedBall"}|>, "LayoutOptions" -> <|"PanelPlotLayout" -> <||>, "PlotRange" -> {{Rational[-1, 2048], Rational[1, 2048]}, {-7.98873038582587*^-24, 7.165015235367508*^-24}}, "Frame" -> {{False, False}, {False, False}}, "AxesOrigin" -> {0, 2.2204459609442666`*^-16}, "ImageSize" -> {360, 360/GoldenRatio}, "Axes" -> {True, True}, "LabelStyle" -> {}, "AspectRatio" -> GoldenRatio^(-1), "DefaultStyle" -> {Directive[Opacity[1.], RGBColor[0.368417, 0.506779, 0.709798], AbsoluteThickness[2]]}, "HighlightLabelingFunctions" -> <|"CoordinatesToolOptions" -> ({Identity[Part[#, 1]], Identity[Part[#, 2]]}& ), "ScalingFunctions" -> {{Identity, Identity}, {Identity, Identity}}|>, "Primitives" -> {}, "GCFlag" -> False|>, "Meta" -> <|"DefaultHighlight" -> {"Dynamic", None}, "Index" -> {}, "Function" -> Plot, "GroupHighlight" -> False|>|>, "DynamicHighlight"], AspectRatio -> GoldenRatio^(-1), Axes -> {True, True}, AxesLabel -> {None, None}, AxesOrigin -> {0, 2.2204459609442666`*^-16}, DisplayFunction -> Identity, Frame -> {{False, False}, {False, False}}, FrameLabel -> {{None, None}, {None, None}}, FrameTicks -> {{Automatic, Automatic}, {Automatic, Automatic}}, GridLines -> {None, None}, GridLinesStyle -> Directive[GrayLevel[0.5, 0.4]], ImagePadding -> All, Method -> {"DefaultBoundaryStyle" -> Automatic, "DefaultGraphicsInteraction" -> {"Version" -> 1.2, "TrackMousePosition" -> {True, False}, "Effects" -> {"Highlight" -> {"ratio" -> 2}, "HighlightPoint" -> {"ratio" -> 2}, "Droplines" -> {"freeformCursorMode" -> True, "placement" -> {"x" -> "All", "y" -> "None"}}}}, "DefaultMeshStyle" -> AbsolutePointSize[6], "ScalingFunctions" -> None, "CoordinatesToolOptions" -> {"DisplayFunction" -> ({(Identity[#]& )[Part[#, 1]], (Identity[#]& )[Part[#, 2]]}& ), "CopiedValueFunction" -> ({(Identity[#]& )[Part[#, 1]], (Identity[#]& )[Part[#, 2]]}& )}}, PlotRange -> {{Rational[-1, 2048], Rational[1, 2048]}, {-7.98873038582587*^-24, 7.165015235367508*^-24}}, PlotRangeClipping -> True, PlotRangePadding -> {{Scaled[0.02], Scaled[0.02]}, {Scaled[0.05], Scaled[0.05]}}, Ticks -> {Automatic, Automatic}],Graphics[Annotation[{{{{}, {}, Annotation[{Directive[Opacity[1.], RGBColor[0.368417, 0.506779, 0.709798], AbsoluteThickness[2]], Line[CompressedData["*)
(*1:eJwB4QQe+yFib1JlAgAAAE0AAAACAAAALj8W6v//L7+uHPlqbyHxOm/5I4j5*)
(*+i+/zLYyDyUA+zqwszEm8/Uvv4Y04uzuF/w6MihNYubrL78NhKq3jU7SOjcR*)
(*hNrM1y+/S25ej99j57pA4/HKma8vv1p2eIB9zsq6UofNqzNfL78Aj7z0EZzq*)
(*unfPhG1nvi6/1XvGJWZd9rrU1AvqwWEtvwSaZh1x0OU6HcBvFDccLL+7EEOc*)
(*PNv0OrVrxAUP3Sq/+GAmsTBuwLpgjDHw2YIpvy4eIfXltvk6+JJ7iL8/KL9D*)
(*UOtgNqjuOqMO3hmY4Sa/nNGpcLrJ0zqcSjFy04klvxgZp8t+9+i6gmxheClJ*)
(*JL8yMvkZwKbZunsDqndy7SK/5+UrY1UOljphgM8k1qghvwRk/lXI+8O6lb3l*)
(*mJxqIL8H4Tm72/zturrfKAysIh6/+Cu3LNuhw7ojEEBCVJ4bv9hDYe2Fp9A6*)
(*siqIauLjGL815WXAHl/XOt/FsSA2Nha/H4nBcNAztTrlLJUyv7YTv8Yet8i1*)
(*T8M6EX6pNi4BEb/5i/17cc26Oiw27yyl8wy/YAfv6ddLrjqDRO3QuXgHvwVR*)
(*9cAToqk6FVSukFkXAr/FBx6CxoOIurL2xQ/IJPq+oLVcSclEgTql2+WN6YTu*)
(*vpeqPc/HRGS6XVEctNtj1L6siT64QKFCOrkKe9XCddM+0NTSnsNnQjo14gEe*)
(*gVjvPg7OhRNNp0G6I/C7+TrC+T7GlBbEpY2dukmj2c1OWAI/7cjxfEXpqTpF*)
(*TRID9bUHP83FYlg0xas6j1/XgDC3DD9trCO4rPazOsdOHQ1QEhE/EN/kjwCB*)
(*wzrtIRV+0poTP58WNefnqc+6dXQrYY8WFj9HxhAZhlKmOtfcEFJmyBg/GV3e*)
(*t5T3xDpgeTznB0wbP5C4ehsrNaS6wys3isMFHj9307ATwqPPuiYJvOikSCA/*)
(*a58DMDYu5LodvGtFBYghP++eXF3yF9I6APoCqXLiIj/4ikSR2vfjOvdRvV7F*)
(*JSQ/JNpLj2Ux6jrbNF8bJYQlP5eGZ2Q9EOa6cFcQESLcJj9r6M7Mbvr+OhmU*)
(*5FgEHSg/P4L+zr6k1bquW6Cn83gpP5HsxmT+G/A6Vz1/SMi9Kj9vfl0dcXnh*)
(*OrFebSI6/Cs/JGb+QML68zr4CkMDuVUtP2hCj7cGBfU6UtE7Nh2YLj/3IOiy*)
(*2Dzyugo7C8G8nS4/8EKaNXmF/TrBpNpLXKMuP0/x3dgpsJs6MHh5YZuuLj+2*)
(*6DIwz+rlOg4ft4wZxS4/LyKjHqBy7brJbDLjFfIuP8mB3QKFXu+6QAgpkA5M*)
(*Lz97ySWgzZD2uvhx+BquUS8/kjFwoSpX7rqv28elTVcvP4Dmy7pxev46Hq9m*)
(*u4xiLz9/L48M9/nsOvxVpOYKeS8/6SRJZi26+Dq3ox89B6YvP/UvqO1Az+i6*)
(*bg3vx6arLz9Mbq31kBnSOiZ3vlJGsS8/zBVQKrfA/TqUSl1ohbwvPwwBtpQ7*)
(*wdS6cvGakwPTLz+f/oxwDmP0Oipbah6j2C8/C9asbvrP3DrhxDmpQt4vP86X*)
(*qMV/TO66UJjYvoHpLz9nYXh2RfvMOggCqEkh7y8/QoHZmwjf+Dq/a3fUwPQv*)
(*P5aAtWoYycC6dtVGX2D6Lz8af+kV+oHiOi4/Fur//y8/rhz5am8h8TrCmGVO*)
(**)
(*"]]}, "Charting`Private`Tag#1"]}}, {}}, <|"HighlightElements" -> <|"Label" -> {"XYLabel"}, "Ball" -> {"InterpolatedBall"}|>, "LayoutOptions" -> <|"PanelPlotLayout" -> <||>, "PlotRange" -> {{Rational[-1, 4096], Rational[1, 4096]}, {-1.1666161206489614`*^-24, 1.6015381522847428`*^-24}}, "Frame" -> {{False, False}, {False, False}}, "AxesOrigin" -> {0, 2.2204460360462934`*^-16}, "ImageSize" -> {360, 360/GoldenRatio}, "Axes" -> {True, True}, "LabelStyle" -> {}, "AspectRatio" -> GoldenRatio^(-1), "DefaultStyle" -> {Directive[Opacity[1.], RGBColor[0.368417, 0.506779, 0.709798], AbsoluteThickness[2]]}, "HighlightLabelingFunctions" -> <|"CoordinatesToolOptions" -> ({Identity[Part[#, 1]], Identity[Part[#, 2]]}& ), "ScalingFunctions" -> {{Identity, Identity}, {Identity, Identity}}|>, "Primitives" -> {}, "GCFlag" -> False|>, "Meta" -> <|"DefaultHighlight" -> {"Dynamic", None}, "Index" -> {}, "Function" -> Plot, "GroupHighlight" -> False|>|>, "DynamicHighlight"], AspectRatio -> GoldenRatio^(-1), Axes -> {True, True}, AxesLabel -> {None, None}, AxesOrigin -> {0, 2.2204460360462934`*^-16}, DisplayFunction -> Identity, Frame -> {{False, False}, {False, False}}, FrameLabel -> {{None, None}, {None, None}}, FrameTicks -> {{Automatic, Automatic}, {Automatic, Automatic}}, GridLines -> {None, None}, GridLinesStyle -> Directive[GrayLevel[0.5, 0.4]], ImagePadding -> All, Method -> {"DefaultBoundaryStyle" -> Automatic, "DefaultGraphicsInteraction" -> {"Version" -> 1.2, "TrackMousePosition" -> {True, False}, "Effects" -> {"Highlight" -> {"ratio" -> 2}, "HighlightPoint" -> {"ratio" -> 2}, "Droplines" -> {"freeformCursorMode" -> True, "placement" -> {"x" -> "All", "y" -> "None"}}}}, "DefaultMeshStyle" -> AbsolutePointSize[6], "ScalingFunctions" -> None, "CoordinatesToolOptions" -> {"DisplayFunction" -> ({(Identity[#]& )[Part[#, 1]], (Identity[#]& )[Part[#, 2]]}& ), "CopiedValueFunction" -> ({(Identity[#]& )[Part[#, 1]], (Identity[#]& )[Part[#, 2]]}& )}}, PlotRange -> {{Rational[-1, 4096], Rational[1, 4096]}, {-1.1666161206489614`*^-24, 1.6015381522847428`*^-24}}, PlotRangeClipping -> True, PlotRangePadding -> {{Scaled[0.02], Scaled[0.02]}, {Scaled[0.05], Scaled[0.05]}}, Ticks -> {Automatic, Automatic}]}*)


(* ::Text:: *)
(*In this case, the errors returned by FindMaximum for radii other than 1/512 are not to be believed:*)


(* ::Input:: *)
(*cosGeneralApproximationExtrema1=Map[{#[[2,2,2]],#[[2,1,2]]}&,cosGeneralApproximationResults1]*)


(* ::Input:: *)
(*cosGeneralApproximationMachineExtrema1=MapThread[*)
(*FindMaximum[{Cos[x]-1-machineEvaluate[#[[2]][x]],x>=99#2[[2]]/100&&x<=#1[[1]]},{x,#2[[2]]},WorkingPrecision->30]&,*)
(*{cosGeneralApproximations1,cosGeneralApproximationExtrema1}]*)


(* ::Text:: *)
(*To obtain the number of bits we have to use exhaustive search:*)


(* ::Input:: *)
(*MapThread[{#1[[1]],N[Log[2,Max[Abs[Table[Cos[x]-1-machineEvaluate[#2[[2]][x]],{x,9 #2[[1]]/10,#2[[1]],#2[[1]]/1*^6}]]]],30]}&,{cosGeneralApproximationResults1,cosGeneralApproximations1}]*)


(* ::Text:: *)
(*A radius of 1/1024 gives us 72 bits, which covers the 18 bits we want to have in the accurate values.  Remember that for the same radius the Sin approximation has 84 bits.*)


(* ::Subsubsection:: *)
(*Degree  2*)


(* ::Input:: *)
(*cosFn2[t_]:=(Cos[t]-1+t^2/2)/t^4*)


(* ::Input:: *)
(*Limit[cosFn2[t],t->0]*)


(* ::Input:: *)
(*Series[cosFn2[t],{t,0,2}]*)


(* ::Input:: *)
(*cosGeneralApproximationResults2=Table[*)
(*{2^-n,*)
(*GeneralMiniMaxApproximation[*)
(*{t^2,If[t==0,1/24,cosFn2[t]],If[t==0,1,1/t^4]},*)
(*{t,{0,2^-n},1,0},*)
(*x,WorkingPrecision->40]},*)
(*{n,7,8}]*)


(* ::Input:: *)
(*cosGeneralApproximations2=Map[*)
(*{#[[1]],Function[u, Evaluate[u^4(HornerForm[#[[2,2,1]]]/.x->u^2)]]}&,*)
(*cosGeneralApproximationResults2]*)


(* ::Text:: *)
(*In this case the polynomials don't exhibit numeric noise:*)


(* ::Input:: *)
(*Map[Plot[Cos[x]-1+x^2/2-machineEvaluate[#[[2]][x]],{x,-#[[1]],#[[1]]},PlotRange->Full,WorkingPrecision->30]&,cosGeneralApproximations2]*)


(* ::Output:: *)
(*{Graphics[Annotation[{{{{}, {}, Annotation[{Directive[Opacity[1.], RGBColor[0.368417, 0.506779, 0.709798], AbsoluteThickness[2]], Line[CompressedData["*)
(*1:eJwB4QQe+yFib1JlAgAAAE0AAAACAAAALj8W6v//f78Mu4wEm+c/O2/5I4j5*)
(*+n+/cz+szMpGPzuwszEm8/V/vxVw8RcLqj47MihNYubrf7+nBVunpms9OzcR*)
(*hNrM13+/h2ffx9IAOztA4/HKma9/v8g3hmS+ZzY7UofNqzNff78UPfhEzO8r*)
(*O3fPhG1nvn6/wi1gIxI4vbrU1AvqwWF9v0KpWj5NAjS7HcBvFDccfL9xCLkF*)
(*aQA9u7VrxAUP3Xq/nmR5v7nWP7tgjDHw2YJ5v0gFqrr/vj67+JJ7iL8/eL9x*)
(*NrUeloE7u6MO3hmY4Xa/UKFo2o3qNrucSjFy04l1v61GknOdMTK7gmxheClJ*)
(*dL89DL64wSksu3sDqndy7XK/2iPUV/RvJLthgM8k1qhxvy3cHbgFFx27lb3l*)
(*mJxqcL/bf0kOfQcUu7rfKAysIm6/t3XffYBzCbsjEEBCVJ5rv/JH6bpCsv+6*)
(*siqIauLjaL+DCh1mZ8jxut/FsSA2Nma/i7Xc6/Gv4rrlLJUyv7Zjv5n5DhdQ*)
(*2NK6EX6pNi4BYb8EBChxD/m/uiw27yyl81y/by8Lc2nfqLqDRO3QuXhXvz+K*)
(*18zGy4y6FVSukFkXUr/2ci6GiYhourL2xQ/IJEq/2yztljQfPLql2+WN6YQ+*)
(*vzg66cVgKPK5XVEctNtjJL+kOdAgOE1aubkKe9XCdSM/DYJWgjBkU7k14gEe*)
(*gVg/P5vPzZ9mRvW5I/C7+TrCST/hbC3TXbk5ukmj2c1OWFI/8vaueVi3arpF*)
(*TRID9bVXP9M5NoIulo66j1/XgDC3XD+Nl1x5U7inusdOHQ1QEmE/s5A5L7ZX*)
(*wLrtIRV+0ppjP5HoWpIwQNK6dXQrYY8WZj+oRZwlzxniutfcEFJmyGg/rGIe*)
(*I/Fc8bpgeTznB0xrP9XxTbXtuP26wys3isMFbj+t0vqWmfAIuyYJvOikSHA/*)
(*W6/NYsAzE7sdvGtFBYhxP3gnstgjCxy7APoCqXLicj8dN41HPjgku/dRvV7F*)
(*JXQ/MpV/wStQK7vbNF8bJYR1P5Qt7kWsHTK7cFcQESLcdj+JTl6FB9g2uxmU*)
(*5FgEHXg/wl59Wh8VO7uuW6Cn83h5P4xRyidtqz67Vz1/SMi9ej/QELKEHeE/*)
(*u7FebSI6/Hs/vR2inQ+CPbv4CkMDuVV9P7nAQnBIfzS7UtE7Nh2Yfj9rFWwc*)
(*77AHuwo7C8G8nX4/ninq+f1eBLvBpNpLXKN+P62Lv6U+IwG7MHh5YZuufj+R*)
(*p7Y3rcr0ug4ft4wZxX4/CUCiKEUi2jrJbDLjFfJ+PwWM9C0gExA7QAgpkA5M*)
(*fz9bPhA1zzcoO/hx+BquUX8/58lqfEBGKTuv28elTVd/PzVEdwfrWio7Hq9m*)
(*u4xifz9QykMKlpMsO/xVpOYKeX8/4j7SWFqTMDu3ox89B6Z/PxKlRUJhWzU7*)
(*bg3vx6arfz/Us7gxyfg1OyZ3vlJGsX8/4rOtVraZNjuUSl1ohbx/P1GO9ADT*)
(*3Dc7cvGakwPTfz/Nhx8rj3c6Oypbah6j2H8/MsQwCcYeOzvhxDmpQt5/P8OS*)
(*5E3cyjs7UJjYvoHpfz+OvX/23x89OwgCqEkh738/E4OKGODOPTu/a3fUwPR/*)
(*PymC/UmGfj47dtVGX2D6fz8om0m+2i8/Oy4/Fur//38/DLuMBJvnPzsyX06a*)
(**)
(*"]]}, "Charting`Private`Tag#1"]}}, {}}, <|"HighlightElements" -> <|"Label" -> {"XYLabel"}, "Ball" -> {"InterpolatedBall"}|>, "LayoutOptions" -> <|"PanelPlotLayout" -> <||>, "PlotRange" -> {{Rational[-1, 128], Rational[1, 128]}, {-2.6369985773252234`*^-23, 2.6390956844366166`*^-23}}, "Frame" -> {{False, False}, {False, False}}, "AxesOrigin" -> {0, 2.2204457562388255`*^-16}, "ImageSize" -> {360, 360/GoldenRatio}, "Axes" -> {True, True}, "LabelStyle" -> {}, "AspectRatio" -> GoldenRatio^(-1), "DefaultStyle" -> {Directive[Opacity[1.], RGBColor[0.368417, 0.506779, 0.709798], AbsoluteThickness[2]]}, "HighlightLabelingFunctions" -> <|"CoordinatesToolOptions" -> ({Identity[Part[#, 1]], Identity[Part[#, 2]]}& ), "ScalingFunctions" -> {{Identity, Identity}, {Identity, Identity}}|>, "Primitives" -> {}, "GCFlag" -> False|>, "Meta" -> <|"DefaultHighlight" -> {"Dynamic", None}, "Index" -> {}, "Function" -> Plot, "GroupHighlight" -> False|>|>, "DynamicHighlight"], AspectRatio -> GoldenRatio^(-1), Axes -> {True, True}, AxesLabel -> {None, None}, AxesOrigin -> {0, 2.2204457562388255`*^-16}, DisplayFunction -> Identity, Frame -> {{False, False}, {False, False}}, FrameLabel -> {{None, None}, {None, None}}, FrameTicks -> {{Automatic, Automatic}, {Automatic, Automatic}}, GridLines -> {None, None}, GridLinesStyle -> Directive[GrayLevel[0.5, 0.4]], ImagePadding -> All, Method -> {"DefaultBoundaryStyle" -> Automatic, "DefaultGraphicsInteraction" -> {"Version" -> 1.2, "TrackMousePosition" -> {True, False}, "Effects" -> {"Highlight" -> {"ratio" -> 2}, "HighlightPoint" -> {"ratio" -> 2}, "Droplines" -> {"freeformCursorMode" -> True, "placement" -> {"x" -> "All", "y" -> "None"}}}}, "DefaultMeshStyle" -> AbsolutePointSize[6], "ScalingFunctions" -> None, "CoordinatesToolOptions" -> {"DisplayFunction" -> ({(Identity[#]& )[Part[#, 1]], (Identity[#]& )[Part[#, 2]]}& ), "CopiedValueFunction" -> ({(Identity[#]& )[Part[#, 1]], (Identity[#]& )[Part[#, 2]]}& )}}, PlotRange -> {{Rational[-1, 128], Rational[1, 128]}, {-2.6369985773252234`*^-23, 2.6390956844366166`*^-23}}, PlotRangeClipping -> True, PlotRangePadding -> {{Scaled[0.02], Scaled[0.02]}, {Scaled[0.05], Scaled[0.05]}}, Ticks -> {Automatic, Automatic}],Graphics[Annotation[{{{{}, {}, Annotation[{Directive[Opacity[1.], RGBColor[0.368417, 0.506779, 0.709798], AbsoluteThickness[2]], Line[CompressedData["*)
(*1:eJwB4QQe+yFib1JlAgAAAE0AAAACAAAALj8W6v//b78WrP4NfgnAOm/5I4j5*)
(*+m+/ajg/nJpFvzqwszEm8/Vvvyjj2ZPtuL46MihNYubrb7+womeSVG69OjcR*)
(*hNrM12+/IJy+C+7FujpA4/HKma9vvwZzjzBbYLY6UofNqzNfb78MTSmceSas*)
(*OnfPhG1nvm6/1imiQtr1QTrU1AvqwWFtv+dGdxh+IbS6HcBvFDccbL8S74r7*)
(*Xhe9urVrxAUP3Wq//kmPqtDMv7pgjDHw2YJpv83H44wgtb66+JJ7iL8/aL/T*)
(*r+RvLoC7uqMO3hmY4Wa/hAXtUM/0trqcSjFy04llv/GFgGJmN7K6gmxheClJ*)
(*ZL98U2lPzBysunsDqndy7WK/0J6TKB5kpLphgM8k1qhhv7w3EYugIZ26lb3l*)
(*mJxqYL+vwsKn1gSUurrfKAysIl6/cLT1wfhzibojEEBCVJ5bv/Kis0ODrH+6*)
(*siqIauLjWL/ktjcFQblxut/FsSA2Nla/zI2xPpqwYrrlLJUyv7ZTvwNAepKu*)
(*1VK6EX6pNi4BUb+lvkpDj+k/uiw27yyl80y/U735FrHMKLqDRO3QuXhHv5EB*)
(*q0+zqAy6FVSukFkXQr+x0lPfMULoubL2xQ/IJDq/1VYzpMQ8vbml2+WN6YQu*)
(*v2SMu1Jw2nO5XVEctNtjFL81+fdFSe/guLkKe9XCdRM/WBe7giwp37g14gEe*)
(*gVgvPyimkfjp9XK5I/C7+TrCOT/DIqLTAHG5uUmj2c1OWEI/WNVJlV0y6rlF*)
(*TRID9bVHP3RJk/U72Q66j1/XgDC3TD9BFagLicAnusdOHQ1QElE/T8VGf09n*)
(*QLrtIRV+0ppTPyigBNITLFK6dXQrYY8WVj/2QHmCgCNiutfcEFJmyFg/ot7V*)
(*3pNZcbpgeTznB0xbP1flmYkj1326wys3isMFXj/np+v+BQyJuiYJvOikSGA/*)
(*wkkwIx8yk7odvGtFBYhhP8BJQVKtBpy6APoCqXLiYj9yUh1LQECkuvdRvV7F*)
(*JWQ/Avor2slQq7rbNF8bJYRlPxSrKrabD7K6cFcQESLcZj9uZbVK7NS2uhmU*)
(*5FgEHWg/3WwkDD0Ju7quW6Cn83hpP4xVeF1pvr66Vz1/SMi9aj9xL56pK9G/*)
(*urFebSI6/Gs/UPSUDB9qvbr4CkMDuVVtP1C3Ed7cj7S6UtE7Nh2Ybj+jyqai*)
(*JkOIugo7C8G8nW4/sIaSmuragrrBpNpLXKNuP2QbhVVaioG6MHh5YZuubj89*)
(*7BNdlNh2ug4ft4wZxW4/wN9dJWUjZDrJbDLjFfJuP5mnzFHmy5A6QAgpkA5M*)
(*bz/dbH3ydUyoOvhx+BquUW8/4DDrxhUiqTqv28elTVdvP4jtsm8zfKo6Hq9m*)
(*u4xibz/o7dFfejesOvxVpOYKeW8/8XxRSDWfsDq3ox89B6ZvPzRqo1YtV7U6*)
(*bg3vx6arbz88M/tMIfm1OiZ3vlJGsW8/0ZB79afGtjqUSl1ohbxvP4KNh2Bg*)
(*4Lc6cvGakwPTbz/Mo618gJ+6Oipbah6j2G8/CK0HLfgPuzrhxDmpQt5vPwEZ*)
(*6GzfwLs6UJjYvoHpbz/8wOHE0+W8OggCqEkh728/zmWRTLuyvTq/a3fUwPRv*)
(*PzMFnLIbbL46dtVGX2D6bz/CbEAkbge/Oi4/Fur//28/Fqz+DX4JwDqgZGWh*)
(**)
(*"]]}, "Charting`Private`Tag#1"]}}, {}}, <|"HighlightElements" -> <|"Label" -> {"XYLabel"}, "Ball" -> {"InterpolatedBall"}|>, "LayoutOptions" -> <|"PanelPlotLayout" -> <||>, "PlotRange" -> {{Rational[-1, 256], Rational[1, 256]}, {-1.0280650589930651`*^-25, 1.0363719845969535`*^-25}}, "Frame" -> {{False, False}, {False, False}}, "AxesOrigin" -> {0, 2.220446048107562*^-16}, "ImageSize" -> {360, 360/GoldenRatio}, "Axes" -> {True, True}, "LabelStyle" -> {}, "AspectRatio" -> GoldenRatio^(-1), "DefaultStyle" -> {Directive[Opacity[1.], RGBColor[0.368417, 0.506779, 0.709798], AbsoluteThickness[2]]}, "HighlightLabelingFunctions" -> <|"CoordinatesToolOptions" -> ({Identity[Part[#, 1]], Identity[Part[#, 2]]}& ), "ScalingFunctions" -> {{Identity, Identity}, {Identity, Identity}}|>, "Primitives" -> {}, "GCFlag" -> False|>, "Meta" -> <|"DefaultHighlight" -> {"Dynamic", None}, "Index" -> {}, "Function" -> Plot, "GroupHighlight" -> False|>|>, "DynamicHighlight"], AspectRatio -> GoldenRatio^(-1), Axes -> {True, True}, AxesLabel -> {None, None}, AxesOrigin -> {0, 2.220446048107562*^-16}, DisplayFunction -> Identity, Frame -> {{False, False}, {False, False}}, FrameLabel -> {{None, None}, {None, None}}, FrameTicks -> {{Automatic, Automatic}, {Automatic, Automatic}}, GridLines -> {None, None}, GridLinesStyle -> Directive[GrayLevel[0.5, 0.4]], ImagePadding -> All, Method -> {"DefaultBoundaryStyle" -> Automatic, "DefaultGraphicsInteraction" -> {"Version" -> 1.2, "TrackMousePosition" -> {True, False}, "Effects" -> {"Highlight" -> {"ratio" -> 2}, "HighlightPoint" -> {"ratio" -> 2}, "Droplines" -> {"freeformCursorMode" -> True, "placement" -> {"x" -> "All", "y" -> "None"}}}}, "DefaultMeshStyle" -> AbsolutePointSize[6], "ScalingFunctions" -> None, "CoordinatesToolOptions" -> {"DisplayFunction" -> ({(Identity[#]& )[Part[#, 1]], (Identity[#]& )[Part[#, 2]]}& ), "CopiedValueFunction" -> ({(Identity[#]& )[Part[#, 1]], (Identity[#]& )[Part[#, 2]]}& )}}, PlotRange -> {{Rational[-1, 256], Rational[1, 256]}, {-1.0280650589930651`*^-25, 1.0363719845969535`*^-25}}, PlotRangeClipping -> True, PlotRangePadding -> {{Scaled[0.02], Scaled[0.02]}, {Scaled[0.05], Scaled[0.05]}}, Ticks -> {Automatic, Automatic}]}*)


(* ::Text:: *)
(*And we see that the error estimate from the Remez algorithm is very close to the one obtained by FindMaximum:*)


(* ::Input:: *)
(*cosGeneralApproximationExtrema2=Map[{#[[2,2,2]],#[[2,1,2]]}&,cosGeneralApproximationResults2]*)


(* ::Input:: *)
(*cosGeneralApproximationMachineExtrema2=MapThread[*)
(*FindMaximum[{Cos[x]-1+x^2/2-machineEvaluate[#[[2]][x]],x>=99#2[[2]]/100&&x<=#1[[1]]},{x,#2[[2]]},WorkingPrecision->30]&,*)
(*{cosGeneralApproximations2,cosGeneralApproximationExtrema2}]*)


(* ::Text:: *)
(*This gives us the number of bits available for each radius.  1/128 is not appropriate because it could only use 11 zeroes after the mantissa, but all the other radii work since they can use at least 18 bits after the mantissa:*)


(* ::Input:: *)
(*MapThread[{#1[[1]],Log[2,Abs[#2[[1]]]]}&,{cosGeneralApproximationResults2,cosGeneralApproximationMachineExtrema2}]*)


(* ::Text:: *)
(*Both radii give us enough bits, but we have to remember that 1/128 does not work for the Sin function.  So if we use degree 2 for Cos, we have to use a radius of 1/256.*)


(* ::Section:: *)
(*Error Analysis*)


(* ::Subsection:: *)
(*Stehle\:0301-Zimmermann Polynomials*)


(* ::Input:: *)
(*szHMax=CorrectlyRound[2^-10.977`30,RoundingMode->Toward0]*)


(* ::Input:: *)
(*szHMin=-szHMax*)


(* ::Input:: *)
(*szHInterval=Interval[{szHMin,szHMax}]*)


(* ::Input:: *)
(*machineEvaluateInterval[szHMax]//bounds//N*)


(* ::Input:: *)
(*machineEvaluateInterval[szHMax^2]//bounds//N*)


(* ::Input:: *)
(*machineEvaluateInterval[szHMax^3]//bounds//N*)


(* ::Input:: *)
(*a5=1/5!;*)


(* ::Input:: *)
(*machineEvaluateInterval[a5 szHMax^2]//bounds//N*)


(* ::Input:: *)
(*a3=1/3!;*)


(* ::Input:: *)
(*machineEvaluateInterval[a5 szHMax^2-a3]//bounds//N*)


(* ::Input:: *)
(*machineEvaluateInterval[szHMax^3(a5 szHMax^2-a3)]//bounds//N*)


(* ::Input:: *)
(*a4=1/4!;*)


(* ::Input:: *)
(*machineEvaluateInterval[a4 szHMax^2]//bounds//N*)


(* ::Input:: *)
(*machineEvaluateInterval[a4 szHMax^2-1/2]//bounds//N*)


(* ::Input:: *)
(*machineEvaluateInterval[szHMax^2(a4 szHMax^2-1/2)]//bounds//N*)


(* ::Section::Closed:: *)
(*Cutoff Points*)


(* ::Text::RGBColor[1, 0, 0]:: *)
(*Work in progress!*)


(* ::Input:: *)
(*SetRoundingMode[TowardPositiveInfinity]*)


(* ::Input:: *)
(*multiTableCutoffs=Table[CorrectlyRound[ArcSin[2^-n]],{n,1,10}]*)


(* ::Input:: *)
(*multiTableCutoffs//N*)


(* ::Input:: *)
(*N[Sin[multiTableCutoffs],30]*)


(* ::Input:: *)
(*HexLiteral[multiTableCutoffs]*)


(* ::Input:: *)
(*cutoffProbabilities=(Drop[Prepend[multiTableCutoffs,\[Pi]/4],-1]-multiTableCutoffs)/(\[Pi]/4)*)


(* ::Input:: *)
(*cutoffExpectation=Total[MapIndexed[#1*#2[[1]]&,cutoffProbabilities]]*)


(* ::Input:: *)
(*N[{cutoffExpectation,cutoffProbabilities}]*)


(* ::Text:: *)
(*Comparisons with the nested loop are a bit more expensive:*)


(* ::Input:: *)
(*nlComparisons=Total[MapIndexed[#1*#2[[1]]&,cutoffProbabilities[[3;;]]]]*)


(* ::Input:: *)
(*nlProbability=Total[cutoffProbabilities[[3;;]]]*)


(* ::Input:: *)
(*nlCutoffExpectation=nlComparisons/nlProbability*)


(* ::Input:: *)
(*N[{nlCutoffExpectation,nlComparisons,nlProbability}]*)


(* ::Input:: *)
(*nlCost=1+nlComparisons+(1-nlProbability)//N*)


(* ::Input:: *)
(*Table[(ArcSin[2^-(n-1+1)]-ArcSin[2^-(n+1)])*2^(9+n),{n,1,8}]//N*)
