(* ::Package:: *)

(* ::Section:: *)
(*Declarations*)


(* ::Input:: *)
(*<<FunctionApproximations`*)


(* ::Input:: *)
(*Get[FileNameJoin[{NotebookDirectory[],"floating_point.wl"}]]*)


(* ::Input:: *)
(*SetFloatingPointFormat[binary64]*)
(*SetRoundingMode[NearestTiesToEven]*)


(* ::Input:: *)
(*machineEvaluate//ClearAll*)


(* ::Input:: *)
(*SetAttributes[machineEvaluate,HoldAll]*)


(* ::Input:: *)
(*machineEvaluate[a_*b_+c_+d__]:=machineEvaluate[CorrectlyRound[machineEvaluate[a]machineEvaluate[b]+machineEvaluate[c]]+d]*)
(*machineEvaluate[a_*b_+c_]:=CorrectlyRound[machineEvaluate[a]machineEvaluate[b]+machineEvaluate[c]]*)
(*machineEvaluate[a_ +Shortest[ b_]]:=CorrectlyRound[machineEvaluate[a]+machineEvaluate[b]]*)
(*machineEvaluate[a_*Shortest[b_]]:=CorrectlyRound[machineEvaluate[a]machineEvaluate[b]]*)
(*machineEvaluate[a_-b_]:=CorrectlyRound[machineEvaluate[a]-machineEvaluate[b]]*)
(*machineEvaluate[a_/b_]:=CorrectlyRound[machineEvaluate[a]/machineEvaluate[b]]*)
(*machineEvaluate[a_^2]:=CorrectlyRound[machineEvaluate[a ]machineEvaluate[a]]*)
(*machineEvaluate[a_^3]:=CorrectlyRound[machineEvaluate[a^2 ]machineEvaluate[a]]*)
(*machineEvaluate[a_^4]:=CorrectlyRound[machineEvaluate[a^2 ]machineEvaluate[a^2]]*)
(*machineEvaluate[a_?NumberQ]:=CorrectlyRound[a]*)
(*machineEvaluate[CorrectlyRound[a_]]:=CorrectlyRound[a]*)


(* ::Input:: *)
(*machineEvaluate[x^3*(-0.1666666666666666666666666517210414159141306579306540343219`30. +0.0083333331609395193764627166673865023798916269058631840465`30. x)]*)


(* ::Section:: *)
(*Sin*)


(* ::Input:: *)
(*Series[Sin[t],{t,0,5}]*)


(* ::Input:: *)
(*sinFn[t_]:=(Sin[t]-t)/t^3*)


(* ::Input:: *)
(*Limit[sinFn[t],t->0]*)


(* ::Input:: *)
(*Series[sinFn[t],{t,0,2}]*)


(* ::Subsection:: *)
(*Minimax  Polynomials Between Table Entries*)


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


(* ::Input:: *)
(*Map[Plot[Sin[x]-x-machineEvaluate[#[[2]][x]],{x,-#[[1]],#[[1]]},PlotRange->Full,WorkingPrecision->30]&,sinGeneralApproximations]*)


(* ::Input:: *)
(*Map[*)
(*FindMaximum[{Sin[x]-x-machineEvaluate[#[[2]][x]],x>=-#[[1]]&&x<=#[[1]]},{x,9#[[1]]/10},WorkingPrecision->30]&,*)
(*sinGeneralApproximations]*)


(* ::Input:: *)
(*N[Log[2,Max[Abs[Table[Sin[x]-x-machineEvaluate[sinGeneralApproximations[[2,2]][x]],{x,9 sinGeneralApproximations[[2,1]]/10,sinGeneralApproximations[[2,1]],sinGeneralApproximations[[2,1]]/1*^6}]]]],30]*)


(* ::Input:: *)
(*machineEvaluate[Evaluate[sinGeneralApproximations[[2,2]][x]]]*)


(* ::Input:: *)
(*CorrectlyRound[CoefficientList[sinGeneralApproximations[[2,2]][x],x]]*)


(* ::Input:: *)
(*HexLiteral[CorrectlyRound[CoefficientList[sinGeneralApproximations[[2,2]][x],x]]]*)


(* ::Input:: *)
(*N[Log[2,Max[Abs[Table[Sin[x]-x-machineEvaluate[sinGeneralApproximations[[4,2]][x]],{x,9 sinGeneralApproximations[[4,1]]/10,sinGeneralApproximations[[4,1]],sinGeneralApproximations[[4,1]]/1*^6}]]]],30]*)


(* ::Input:: *)
(*machineEvaluate[Evaluate[sinGeneralApproximations[[4,2]][x]]]*)


(* ::Input:: *)
(*CorrectlyRound[CoefficientList[sinGeneralApproximations[[4,2]][x],x]]*)


(* ::Input:: *)
(*HexLiteral[CorrectlyRound[CoefficientList[sinGeneralApproximations[[4,2]][x],x]]]*)


(* ::Input:: *)
(*Map[N[Sin[x]-x-machineEvaluate[#[[2]][x]]/.x->#[[1]],20]&,sinGeneralApproximations]*)


(* ::Input:: *)
(*Map[Log[2,N[Abs[Sin[x]-x-machineEvaluate[#[[2]][x]]/.x->#[[1]]],20]]&,sinGeneralApproximations]*)


(* ::Text:: *)
(*A radius of 1/256 gives us 71 bits, which is 18 bits more than the mantissa.  A radius of 1/1024 gives us 84 bits, which covers the 9 extra bits to reach the vicinity of zero.*)


(* ::Subsection:: *)
(*Minimax  Polynomials  Near  Zero*)


(* ::Subsubsection:: *)
(*Degree  1*)


(* ::Input:: *)
(*sinGeneralApproximation0Results=Map[{#,GeneralMiniMaxApproximation[*)
(*{t^2,If[t==0,-1/6,sinFn[t]],If[t==0,1,Sin[t]/t^3]},*)
(*{t,{0,#},1,0},*)
(*x,WorkingPrecision->30]}&,{1/512,1/1024}]*)


(* ::Input:: *)
(*sinGeneralApproximation0=Map[{#[[1]],Function[u, Evaluate[ u^3 (#[[2,2,1]]/.x->u^2)]]}&,sinGeneralApproximation0Results]*)


(* ::Input:: *)
(*Map[Plot[Sin[x]-x-machineEvaluate[#[[2]][x]],{x,-#[[1]],#[[1]]},WorkingPrecision->30,PlotRange->Full]&,sinGeneralApproximation0]*)


(* ::Input:: *)
(*Map[Plot[1-(x+machineEvaluate[#[[2]][x]])/Sin[x],{x,-#[[1]],#[[1]]},WorkingPrecision->30,PlotRange->Full]&,sinGeneralApproximation0]*)


(* ::Input:: *)
(*Map[FindMaximum[{1-(x+machineEvaluate[#[[2]][x]])/Sin[x],x>=-#[[1]]&&x<=#[[1]]},{x,9#[[1]]/10},WorkingPrecision->30]&,sinGeneralApproximation0]*)


(* ::Input:: *)
(*N[Log[2,Max[Abs[Table[1-(x+machineEvaluate[sinGeneralApproximation0[[2,2]][x]])/Sin[x],{x,9 sinGeneralApproximation0[[2,1]]/10,sinGeneralApproximation0[[2,1]],sinGeneralApproximation0[[2,1]]/1*^6}]]]],30]*)


(* ::Input:: *)
(*HexLiteral[CorrectlyRound[CoefficientList[sinGeneralApproximations[[2,2]][x],x]]]*)


(* ::Input:: *)
(*Map[N[Log[2,Abs[1-(x+machineEvaluate[#[[2]][x]])/Sin[x]/.x->#[[1]]]],30]&,sinGeneralApproximation0]*)


(* ::Subsubsection:: *)
(*Degree  2*)


(* ::Input:: *)
(*sinGeneralApproximation0Results2=Map[{#,GeneralMiniMaxApproximation[*)
(*{t^2,If[t==0,-1/6,sinFn[t]],If[t==0,1,Sin[t]/t^3]},*)
(*{t,{0,#},2,0},*)
(*x,WorkingPrecision->30]}&,{1/512}]*)


(* ::Input:: *)
(*sinGeneralApproximation02=Map[{#[[1]],Function[u, Evaluate[ u^3 (#[[2,2,1]]/.x->u^2)]]}&,sinGeneralApproximation0Results2]*)


(* ::Input:: *)
(*Map[Plot[Sin[x]-x-machineEvaluate[#[[2]][x]],{x,-#[[1]],#[[1]]},WorkingPrecision->30,PlotRange->Full]&,sinGeneralApproximation02]*)


(* ::Input:: *)
(*Map[Plot[1-(x+machineEvaluate[#[[2]][x]])/Sin[x],{x,-#[[1]],#[[1]]},WorkingPrecision->30,PlotRange->Full]&,sinGeneralApproximation02]*)


(* ::Input:: *)
(*Map[FindMaximum[{1-(x+machineEvaluate[#[[2]][x]])/Sin[x],x>=-#[[1]]&&x<=#[[1]]},{x,9#[[1]]/10},WorkingPrecision->30]&,sinGeneralApproximation02]*)


(* ::Input:: *)
(*N[Log[2,Max[Abs[Table[1-(x+machineEvaluate[sinGeneralApproximation02[[1,2]][x]])/Sin[x],{x,9 sinGeneralApproximation02[[1,1]]/10,sinGeneralApproximation02[[1,1]],sinGeneralApproximation02[[1,1]]/1*^6}]]]],30]*)


(* ::Input:: *)
(*Map[N[Log[2,Abs[1-(x+machineEvaluate[#[[2]][x]])/Sin[x]/.x->#[[1]]]],30]&,sinGeneralApproximation02]*)


(* ::Section:: *)
(*Cos*)


(* ::Input:: *)
(*Series[Cos[t],{t,0,5}]*)


(* ::Subsection:: *)
(*Minimax  Polynomials Between Table Entries*)


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
(*x,WorkingPrecision->40]},*)
(*{n,7,12}]*)


(* ::Input:: *)
(*cosGeneralApproximations1=Map[*)
(*{#[[1]],Function[u, Evaluate[ u^2(#[[2,2,1]]/.x->u^2)]]}&,*)
(*cosGeneralApproximationResults1]*)


(* ::Input:: *)
(*Map[Plot[Cos[x]-1-machineEvaluate[#[[2]][x]],{x,-#[[1]],#[[1]]},PlotRange->Full,WorkingPrecision->30]&,cosGeneralApproximations1]*)


(* ::Output:: *)
(*{Graphics[Annotation[{{{{}, {}, Annotation[{Directive[Opacity[1.], RGBColor[0.368417, 0.506779, 0.709798], AbsoluteThickness[2]], Line[CompressedData["*)
(*1:eJwB4QQe+yFib1JlAgAAAE0AAAACAAAALj8W6v//f7/aNY8C0kaDvG/5I4j5*)
(*+n+/JwdSW7cBg7ywszEm8/V/vz/FqTkpvYK8MihNYubrf79UvGAVqjSCvDcR*)
(*hNrM13+/QiGnOIgngbxA4/HKma9/v4xb+S9IOn68UofNqzNff7+1BK6KY4p2*)
(*vHfPhG1nvn6/3w/mSOfvYbzU1AvqwWF9v7ilPZZTbGo8HcBvFDccfL8Y5eOf*)
(*W4t6PLVrxAUP3Xq/IYBpSEEfgTxgjDHw2YJ5v0bbm6DZAoM8+JJ7iL8/eL9L*)
(*FLtaSzODPKMO3hmY4Xa/BKuGhME3gjycSjFy04l1v4m+QPnMeoA8gmxheClJ*)
(*dL8qrmdXkvF8PHsDqndy7XK/AJL4kNhAeDxhgM8k1qhxvya/yj5h6HM8lb3l*)
(*mJxqcL8YiOlygNJvPLrfKAysIm6/cochRn0RaDwjEEBCVJ5rv3kAI8VA4GE8*)
(*siqIauLjaL9wAJissb9YPN/FsSA2Nma/sWE403lYUDzlLJUyv7Zjv7p4mySG*)
(*9EQ8EX6pNi4BYb+3X/uBA+c3PCw27yyl81y/o3yrKBC0KTyDRO3QuXhXvx6A*)
(*yH9yqxY8FVSukFkXUr8s03MM9D8APLL2xQ/IJEq/yC3GPyfo4Tul2+WN6YQ+*)
(*vxenm/xitrA7XVEctNtjJL/TJUBiISxKO7kKe9XCdSM/BKjVVmF8RTs14gEe*)
(*gVg/P8dzB8l3mbI7I/C7+TrCST9cO3407OHgO0mj2c1OWFI/3y/JIg0tATxF*)
(*TRID9bVXP+yHWrcalBc8j1/XgDC3XD+wTs3Q7OUoPMdOHQ1QEmE/Jqn7Nm5E*)
(*ODztIRV+0ppjP8yICV2uhUQ8dXQrYY8WZj9zb9wqTQNQPNfcEFJmyGg/YpAz*)
(*mUFeWDxgeTznB0xrP2aqJF53KWE8wys3isMFbj+QVE3JlsRnPCYJvOikSHA/*)
(*rvpMHfwCbzwdvGtFBYhxP9gMkPhQe3M8APoCqXLicj+kLIeghxp4PPdRvV7F*)
(*JXQ/UnrSKW95fDzbNF8bJYR1PwLAgIdNcoA8cFcQESLcdj9gFCAM2jGCPBmU*)
(*5FgEHXg/WLwxIJwmgzyuW6Cn83h5P22dFe7WCYM8Vz1/SMi9ej/Ip/suG2GB*)
(*PLFebSI6/Hs/Irh484iIezz4CkMDuVV9P3MNjhCqo2s8UtE7Nh2Yfj9DVx+/*)
(*iUdYvAo7C8G8nX4/A39eVbv0WbzBpNpLXKN+P6kM3f5nolu8MHh5YZuufj+A*)
(*bhFkHAtfvA4ft4wZxX4/8nn6u936YrzJbDLjFfJ+P3QqucyTKWq8QAgpkA5M*)
(*fz94rB4twc10vPhx+BquUX8/BBbp4+lOdbyv28elTVd/P3uMj7mf0XW8Hq9m*)
(*u4xifz8+DZKTm9l2vPxVpOYKeX8/fRjDJ+PweLy3ox89B6Z/PwAUzOHIRn28*)
(*bg3vx6arfz+nLYwxItV9vCZ3vlJGsX8/ePnvYflkfryUSl1ohbx/PwD6rtWv*)
(*hH+8cvGakwPTfz8murTJSuiAvCpbah6j2H8/sqwv26AygbzhxDmpQt5/P/5z*)
(*onyCfYG8UJjYvoHpfz/IhuE7NhSCvAgCqEkh738/TOl19CNggry/a3fUwPR/*)
(*P12RsgrTrIK8dtVGX2D6fz8kNN1NwvmCvC4/Fur//38/2jWPAtJGg7xwFEvw*)
(**)
(*"]]}, "Charting`Private`Tag#1"]}}, {}}, <|"HighlightElements" -> <|"Label" -> {"XYLabel"}, "Ball" -> {"InterpolatedBall"}|>, "LayoutOptions" -> <|"PanelPlotLayout" -> <||>, "PlotRange" -> {{Rational[-1, 128], Rational[1, 128]}, {-3.3439643412375455`*^-17, 3.330733010085182*^-17}}, "Frame" -> {{False, False}, {False, False}}, "AxesOrigin" -> {0, 1.8489679631747698`*^-16}, "ImageSize" -> {360, 360/GoldenRatio}, "Axes" -> {True, True}, "LabelStyle" -> {}, "AspectRatio" -> GoldenRatio^(-1), "DefaultStyle" -> {Directive[Opacity[1.], RGBColor[0.368417, 0.506779, 0.709798], AbsoluteThickness[2]]}, "HighlightLabelingFunctions" -> <|"CoordinatesToolOptions" -> ({Identity[Part[#, 1]], Identity[Part[#, 2]]}& ), "ScalingFunctions" -> {{Identity, Identity}, {Identity, Identity}}|>, "Primitives" -> {}, "GCFlag" -> False|>, "Meta" -> <|"DefaultHighlight" -> {"Dynamic", None}, "Index" -> {}, "Function" -> Plot, "GroupHighlight" -> False|>|>, "DynamicHighlight"], AspectRatio -> GoldenRatio^(-1), Axes -> {True, True}, AxesLabel -> {None, None}, AxesOrigin -> {0, 1.8489679631747698`*^-16}, DisplayFunction -> Identity, Frame -> {{False, False}, {False, False}}, FrameLabel -> {{None, None}, {None, None}}, FrameTicks -> {{Automatic, Automatic}, {Automatic, Automatic}}, GridLines -> {None, None}, GridLinesStyle -> Directive[GrayLevel[0.5, 0.4]], ImagePadding -> All, Method -> {"DefaultBoundaryStyle" -> Automatic, "DefaultGraphicsInteraction" -> {"Version" -> 1.2, "TrackMousePosition" -> {True, False}, "Effects" -> {"Highlight" -> {"ratio" -> 2}, "HighlightPoint" -> {"ratio" -> 2}, "Droplines" -> {"freeformCursorMode" -> True, "placement" -> {"x" -> "All", "y" -> "None"}}}}, "DefaultMeshStyle" -> AbsolutePointSize[6], "ScalingFunctions" -> None, "CoordinatesToolOptions" -> {"DisplayFunction" -> ({(Identity[#]& )[Part[#, 1]], (Identity[#]& )[Part[#, 2]]}& ), "CopiedValueFunction" -> ({(Identity[#]& )[Part[#, 1]], (Identity[#]& )[Part[#, 2]]}& )}}, PlotRange -> {{Rational[-1, 128], Rational[1, 128]}, {-3.3439643412375455`*^-17, 3.330733010085182*^-17}}, PlotRangeClipping -> True, PlotRangePadding -> {{Scaled[0.02], Scaled[0.02]}, {Scaled[0.05], Scaled[0.05]}}, Ticks -> {Automatic, Automatic}],Graphics[Annotation[{{{{}, {}, Annotation[{Directive[Opacity[1.], RGBColor[0.368417, 0.506779, 0.709798], AbsoluteThickness[2]], Line[CompressedData["*)
(*1:eJwB4QQe+yFib1JlAgAAAE0AAAACAAAALj8W6v//b7+z55pCaEUjvG/5I4j5*)
(*+m+/j6P6DMUFI7ywszEm8/Vvv78EMFQWuiK8MihNYubrb7/GbwKGeDIivDcR*)
(*hNrM12+//5tPs+cjIbxA4/HKma9vv4aw0iFiNB68UofNqzNfb7+TenlI35AW*)
(*vHfPhG1nvm6/BXzRouTzAbzU1AvqwWFtv23naKVrZAo8HcBvFDccbL9Q20ko*)
(*AYUaPLVrxAUP3Wq/d0DbgGwfITxgjDHw2YJpv51Q099TACM8+JJ7iL8/aL8s*)
(*9/vbYjYjPKMO3hmY4Wa/z/6HyH44IjycSjFy04llv5fV61avfCA8gmxheClJ*)
(*ZL9oJiIDB/McPHsDqndy7WK/oA87vA0+GDxhgM8k1qhhv+O/K4Ra6RM8lb3l*)
(*mJxqYL9YW00kqNgPPLrfKAysIl6/2Z1E8P0NCDwjEEBCVJ5bv0UdYlVW3wE8*)
(*siqIauLjWL+JnwD+yMP4O9/FsSA2Nla/ZQvMLIZa8DvlLJUyv7ZTvzK3Fq10*)
(*9OQ7EX6pNi4BUb+1jqy9/+TXOyw27yyl80y/IzX7uYqlyTuDRO3QuXhHv+Ah*)
(*UR/Hx7Y7FVSukFkXQr+7NXsaq1ygO7L2xQ/IJDq/jF93RQ4Ggjul2+WN6YQu*)
(*v3xoM6j/UlA7XVEctNtjFL/uFHvaKP/vOrkKe9XCdRM/uXFXzWX46To14gEe*)
(*gVgvP2DLcXAPxVI7I/C7+TrCOT9m6UrtRuOAO0mj2c1OWEI/+zmNZHg2oTtF*)
(*TRID9bVHP1jLHnwai7c7j1/XgDC3TD9LAEqinNrIO8dOHQ1QElE/g7A9SnU2*)
(*2DvtIRV+0ppTPyKVSJw4iOQ7dXQrYY8WVj+gGy+tEAbwO9fcEFJmyFg/A7s7*)
(*4HZY+DtgeTznB0xbP0LLXH//JgE8wys3isMFXj+7qC9spcEHPCYJvOikSGA/*)
(*CJ7L40oJDzwdvGtFBYhhP2AWQ9xcdxM8APoCqXLiYj9qbdEyaR4YPPdRvV7F*)
(*JWQ/LfKbI216HDzbNF8bJYRlP6frgsHUciA8cFcQESLcZj8WMBurfi4iPBmU*)
(*5FgEHWg/rfjY1RAnIzyuW6Cn83hpP/qXEw7RCyM8Vz1/SMi9aj+I56Yai2Ah*)
(*PLFebSI6/Gs/kqPycn+DGzz4CkMDuVVtP7/Q8kCsrws8UtE7Nh2Ybj9CMwTK*)
(*kl74uwo7C8G8nW4/bsgePDj7+bvBpNpLXKNuP14sRj/6kPu7MHh5YZuubj9v*)
(*LkQ1ufb+uw4ft4wZxW4/G++w8+n0ArzJbDLjFfJuP5LxilMrMwq8QAgpkA5M*)
(*bz8jdjkTacwUvPhx+BquUW8/4nJtIbhKFbyv28elTVdvP5HkdTEw1BW8Hq9m*)
(*u4xibz9E5Hw6KdIWvPxVpOYKeW8/hwsPIT3sGLy3ox89B6ZvP+zMxYg8Tx28*)
(*bg3vx6arbz/6hqdZi9AdvCZ3vlJGsW8/V06b0dlrHryUSl1ohbxvP3EfQVjU*)
(*fx+8cvGakwPTbz/8NtQ/DOogvCpbah6j2G8/X989tzc0IbzhxDmpQt5vP/qn*)
(*Ewi5fCG8UJjYvoHpbz88wDzgsxcivAgCqEkh728/byppSp1gIry/a3fUwPRv*)
(*P2kDlNWlrCK8dtVGX2D6bz8cxbHtvvUivC4/Fur//28/s+eaQmhFI7wXk0F1*)
(**)
(*"]]}, "Charting`Private`Tag#1"]}}, {}}, <|"HighlightElements" -> <|"Label" -> {"XYLabel"}, "Ball" -> {"InterpolatedBall"}|>, "LayoutOptions" -> <|"PanelPlotLayout" -> <||>, "PlotRange" -> {{Rational[-1, 256], Rational[1, 256]}, {-5.22344812098875*^-19, 5.207543915783045*^-19}}, "Frame" -> {{False, False}, {False, False}}, "AxesOrigin" -> {0, 2.2146431015717308`*^-16}, "ImageSize" -> {360, 360/GoldenRatio}, "Axes" -> {True, True}, "LabelStyle" -> {}, "AspectRatio" -> GoldenRatio^(-1), "DefaultStyle" -> {Directive[Opacity[1.], RGBColor[0.368417, 0.506779, 0.709798], AbsoluteThickness[2]]}, "HighlightLabelingFunctions" -> <|"CoordinatesToolOptions" -> ({Identity[Part[#, 1]], Identity[Part[#, 2]]}& ), "ScalingFunctions" -> {{Identity, Identity}, {Identity, Identity}}|>, "Primitives" -> {}, "GCFlag" -> False|>, "Meta" -> <|"DefaultHighlight" -> {"Dynamic", None}, "Index" -> {}, "Function" -> Plot, "GroupHighlight" -> False|>|>, "DynamicHighlight"], AspectRatio -> GoldenRatio^(-1), Axes -> {True, True}, AxesLabel -> {None, None}, AxesOrigin -> {0, 2.2146431015717308`*^-16}, DisplayFunction -> Identity, Frame -> {{False, False}, {False, False}}, FrameLabel -> {{None, None}, {None, None}}, FrameTicks -> {{Automatic, Automatic}, {Automatic, Automatic}}, GridLines -> {None, None}, GridLinesStyle -> Directive[GrayLevel[0.5, 0.4]], ImagePadding -> All, Method -> {"DefaultBoundaryStyle" -> Automatic, "DefaultGraphicsInteraction" -> {"Version" -> 1.2, "TrackMousePosition" -> {True, False}, "Effects" -> {"Highlight" -> {"ratio" -> 2}, "HighlightPoint" -> {"ratio" -> 2}, "Droplines" -> {"freeformCursorMode" -> True, "placement" -> {"x" -> "All", "y" -> "None"}}}}, "DefaultMeshStyle" -> AbsolutePointSize[6], "ScalingFunctions" -> None, "CoordinatesToolOptions" -> {"DisplayFunction" -> ({(Identity[#]& )[Part[#, 1]], (Identity[#]& )[Part[#, 2]]}& ), "CopiedValueFunction" -> ({(Identity[#]& )[Part[#, 1]], (Identity[#]& )[Part[#, 2]]}& )}}, PlotRange -> {{Rational[-1, 256], Rational[1, 256]}, {-5.22344812098875*^-19, 5.207543915783045*^-19}}, PlotRangeClipping -> True, PlotRangePadding -> {{Scaled[0.02], Scaled[0.02]}, {Scaled[0.05], Scaled[0.05]}}, Ticks -> {Automatic, Automatic}],Graphics[Annotation[{{{{}, {}, Annotation[{Directive[Opacity[1.], RGBColor[0.368417, 0.506779, 0.709798], AbsoluteThickness[2]], Line[CompressedData["*)
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


(* ::Input:: *)
(*Map[*)
(*FindMaximum[{Cos[x]-1-machineEvaluate[#[[2]][x]],x>=-#[[1]]&&x<=#[[1]]},{x,9#[[1]]/10},WorkingPrecision->30]&,*)
(*cosGeneralApproximations1]*)


(* ::Input:: *)
(*N[Log[2,Max[Abs[Table[Cos[x]-1-machineEvaluate[cosGeneralApproximations1[[4,2]][x]],{x,9 cosGeneralApproximations1[[4,1]]/10,cosGeneralApproximations1[[4,1]],cosGeneralApproximations1[[4,1]]/1*^6}]]]],30]*)


(* ::Input:: *)
(*machineEvaluate[Evaluate[cosGeneralApproximations1[[4,2]][x]]]*)


(* ::Input:: *)
(*CorrectlyRound[CoefficientList[cosGeneralApproximations1[[4,2]][x],x]]*)


(* ::Input:: *)
(*HexLiteral[CorrectlyRound[CoefficientList[cosGeneralApproximations1[[4,2]][x],x]]]*)


(* ::Input:: *)
(*N[Log[2,Max[Abs[Table[Cos[x]-1-machineEvaluate[cosGeneralApproximations1[[5,2]][x]],{x,9 cosGeneralApproximations1[[5,1]]/10,cosGeneralApproximations1[[5,1]],cosGeneralApproximations1[[5,1]]/1*^6}]]]],30]*)


(* ::Input:: *)
(*machineEvaluate[Evaluate[cosGeneralApproximations1[[5,2]][x]]]*)


(* ::Input:: *)
(*CorrectlyRound[CoefficientList[cosGeneralApproximations1[[5,2]][x],x]]*)


(* ::Input:: *)
(*HexLiteral[CorrectlyRound[CoefficientList[cosGeneralApproximations1[[5,2]][x],x]]]*)


(* ::Input:: *)
(*N[Log[2,Max[Abs[Table[Cos[x]-1-machineEvaluate[cosGeneralApproximations1[[6,2]][x]],{x,9 cosGeneralApproximations1[[6,1]]/10,cosGeneralApproximations1[[6,1]],cosGeneralApproximations1[[6,1]]/1*^6}]]]],30]*)


(* ::Input:: *)
(*machineEvaluate[Evaluate[cosGeneralApproximations1[[6,2]][x]]]*)


(* ::Input:: *)
(*CorrectlyRound[CoefficientList[cosGeneralApproximations1[[6,2]][x],x]]*)


(* ::Input:: *)
(*HexLiteral[CorrectlyRound[CoefficientList[cosGeneralApproximations1[[6,2]][x],x]]]*)


(* ::Input:: *)
(*Map[Log[2,N[Abs[(Cos[x]-1-machineEvaluate[#[[2]][x]]/.x->#[[1]])],20]]&,cosGeneralApproximations1]*)


(* ::Text:: *)
(*A radius of 1/1024 gives us 72 bits, which covers the 18 bits we want to have in the accurate values.  For the binade with 1 leading zeroes, we must use the polynomial for radius 1/2048 works.  For binades with 2 leading zeroes or more, we must use the polynomial for radius 1/4096.*)


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
(*{n,7,10}]*)


(* ::Input:: *)
(*cosGeneralApproximations2=Map[*)
(*{#[[1]],Function[u, Evaluate[u^4(HornerForm[#[[2,2,1]]]/.x->u^2)]]}&,*)
(*cosGeneralApproximationResults2]*)


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
(*"]]}, "Charting`Private`Tag#1"]}}, {}}, <|"HighlightElements" -> <|"Label" -> {"XYLabel"}, "Ball" -> {"InterpolatedBall"}|>, "LayoutOptions" -> <|"PanelPlotLayout" -> <||>, "PlotRange" -> {{Rational[-1, 256], Rational[1, 256]}, {-1.0280650589930651`*^-25, 1.0363719845969535`*^-25}}, "Frame" -> {{False, False}, {False, False}}, "AxesOrigin" -> {0, 2.220446048107562*^-16}, "ImageSize" -> {360, 360/GoldenRatio}, "Axes" -> {True, True}, "LabelStyle" -> {}, "AspectRatio" -> GoldenRatio^(-1), "DefaultStyle" -> {Directive[Opacity[1.], RGBColor[0.368417, 0.506779, 0.709798], AbsoluteThickness[2]]}, "HighlightLabelingFunctions" -> <|"CoordinatesToolOptions" -> ({Identity[Part[#, 1]], Identity[Part[#, 2]]}& ), "ScalingFunctions" -> {{Identity, Identity}, {Identity, Identity}}|>, "Primitives" -> {}, "GCFlag" -> False|>, "Meta" -> <|"DefaultHighlight" -> {"Dynamic", None}, "Index" -> {}, "Function" -> Plot, "GroupHighlight" -> False|>|>, "DynamicHighlight"], AspectRatio -> GoldenRatio^(-1), Axes -> {True, True}, AxesLabel -> {None, None}, AxesOrigin -> {0, 2.220446048107562*^-16}, DisplayFunction -> Identity, Frame -> {{False, False}, {False, False}}, FrameLabel -> {{None, None}, {None, None}}, FrameTicks -> {{Automatic, Automatic}, {Automatic, Automatic}}, GridLines -> {None, None}, GridLinesStyle -> Directive[GrayLevel[0.5, 0.4]], ImagePadding -> All, Method -> {"DefaultBoundaryStyle" -> Automatic, "DefaultGraphicsInteraction" -> {"Version" -> 1.2, "TrackMousePosition" -> {True, False}, "Effects" -> {"Highlight" -> {"ratio" -> 2}, "HighlightPoint" -> {"ratio" -> 2}, "Droplines" -> {"freeformCursorMode" -> True, "placement" -> {"x" -> "All", "y" -> "None"}}}}, "DefaultMeshStyle" -> AbsolutePointSize[6], "ScalingFunctions" -> None, "CoordinatesToolOptions" -> {"DisplayFunction" -> ({(Identity[#]& )[Part[#, 1]], (Identity[#]& )[Part[#, 2]]}& ), "CopiedValueFunction" -> ({(Identity[#]& )[Part[#, 1]], (Identity[#]& )[Part[#, 2]]}& )}}, PlotRange -> {{Rational[-1, 256], Rational[1, 256]}, {-1.0280650589930651`*^-25, 1.0363719845969535`*^-25}}, PlotRangeClipping -> True, PlotRangePadding -> {{Scaled[0.02], Scaled[0.02]}, {Scaled[0.05], Scaled[0.05]}}, Ticks -> {Automatic, Automatic}],Graphics[Annotation[{{{{}, {}, Annotation[{Directive[Opacity[1.], RGBColor[0.368417, 0.506779, 0.709798], AbsoluteThickness[2]], Line[CompressedData["*)
(*1:eJwB4QQe+yFib1JlAgAAAE0AAAACAAAALj8W6v//X7+VI9Y6ZSE8Om/5I4j5*)
(*+l+/C7y6AeGZPjqwszEm8/Vfv3vAUB17wEA6MihNYubrX7+xDqzH39s5OjcR*)
(*hNrM11+/ZESJsq47NzpA4/HKma9fv/FshV/9gjU6UofNqzNfX7/9V/5SyxMk*)
(*OnfPhG1nvl6/kYOtpxsu+jnU1AvqwWFdvyw0+a6/Aja6HcBvFDccXL/aJVU9*)
(*KWE8urVrxAUP3Vq/tgovJ5auP7pgjDHw2YJZvxj/4AXsHz26+JJ7iL8/WL/p*)
(*JDiyna06uqMO3hmY4Va/1YJp+id4N7qcSjFy04lVv4IUQf45tzG6gmxheClJ*)
(*VL9baFwBJhMtunsDqndy7VK/ivgxIkfLI7phgM8k1qhRvy3fVcLPjB26lb3l*)
(*mJxqUL9gIiVZTCoTurrfKAysIk6/r0KMefk0CrojEEBCVJ5Lv7XvkUi5SwC6*)
(*siqIauLjSL9haHf0QhXxud/FsSA2Nka/hTd2Xayw4LnlLJUyv7ZDv3vFvaH9*)
(*fdK5EX6pNi4BQb/u4kvUYu2+uSw27yyl8zy/H+lJurN2prmDRO3QuXg3v5Cp*)
(*xAhgNo25FVSukFkXMr8kIosXTpxqubL2xQ/IJCq/FZhMGc82Qbml2+WN6YQe*)
(*vyUVljzA6w+5XVEctNtjBL8CQ5tmYxmRuLkKe9XCdQM/BMaIWXm6grg14gEe*)
(*gVgfP1NUjtpu8ta4I/C7+TrCKT+zAHkL9ctIuUmj2c1OWDI/zliqz+jKYrlF*)
(*TRID9bU3P93qXKBiIZO5j1/XgDC3PD/gWhHyXXOlucdOHQ1QEkE/1mDKdz23*)
(*wLntIRV+0ppDPxoHzF1NfdG5dXQrYY8WRj8gJ++zabXiudfcEFJmyEg/STEk*)
(*eqFN8blgeTznB0xLP0YCMWHdePy5wys3isMFTj/VwYdKeIUKuiYJvOikSFA/*)
(*HEgSlbG5E7odvGtFBYhRP3pdPorx5Bu6APoCqXLiUj9m4fJKGQckuvdRvV7F*)
(*JVQ/wkO7APSHKrrbNF8bJYRVP30H5VDzETO6cFcQESLcVj8LgscsHB83uhmU*)
(*5FgEHVg/AZFT6N5cO7quW6Cn83hZP/M8WMhqTEC6Vz1/SMi9Wj+wq4VDa2hA*)
(*urFebSI6/Fs/VHGzTX3kPrr4CkMDuVVdP83+XJmTyTK6UtE7Nh2YXj8iUZ2j*)
(*R7L6uQo7C8G8nV4/qJHAMV4X97nBpNpLXKNeP77nS2XU+fk5MHh5YZuuXj8B*)
(*xNpZIY/iuQ4ft4wZxV4/ZOFrdqje87nJbDLjFfJeP96ONIMqvhM6QAgpkA5M*)
(*Xz8kjZO4AqokOvhx+BquUV8/BAouQ482Ijqv28elTVdfPycHvbX1bCk6Hq9m*)
(*u4xiXz/Rjtoh8y8rOvxVpOYKeV8/Q04+R3oQNDq3ox89B6ZfP0LCrsY9MjM6*)
(*bg3vx6arXz/HfnxFXuc5OiZ3vlJGsV8/2459I4WZNTqUSl1ohbxfPxuHjSYu*)
(*+TU6cvGakwPTXz+/ImI1zSk8Oipbah6j2F8/qg2XyWNSODrhxDmpQt5fPwg9*)
(*bPS4xjw6UJjYvoHpXz/EreOQPSs/OggCqEkh718/1eRWXNSIOjq/a3fUwPRf*)
(*P2HOCgajWUA6dtVGX2D6Xz88TSeEgZY8Oi4/Fur//18/lSPWOmUhPDojYkRJ*)
(**)
(*"]]}, "Charting`Private`Tag#1"]}}, {}}, <|"HighlightElements" -> <|"Label" -> {"XYLabel"}, "Ball" -> {"InterpolatedBall"}|>, "LayoutOptions" -> <|"PanelPlotLayout" -> <||>, "PlotRange" -> {{Rational[-1, 512], Rational[1, 512]}, {-4.141932917378821*^-28, 4.228768671840614*^-28}}, "Frame" -> {{False, False}, {False, False}}, "AxesOrigin" -> {0, 2.220446049245711*^-16}, "ImageSize" -> {360, 360/GoldenRatio}, "Axes" -> {True, True}, "LabelStyle" -> {}, "AspectRatio" -> GoldenRatio^(-1), "DefaultStyle" -> {Directive[Opacity[1.], RGBColor[0.368417, 0.506779, 0.709798], AbsoluteThickness[2]]}, "HighlightLabelingFunctions" -> <|"CoordinatesToolOptions" -> ({Identity[Part[#, 1]], Identity[Part[#, 2]]}& ), "ScalingFunctions" -> {{Identity, Identity}, {Identity, Identity}}|>, "Primitives" -> {}, "GCFlag" -> False|>, "Meta" -> <|"DefaultHighlight" -> {"Dynamic", None}, "Index" -> {}, "Function" -> Plot, "GroupHighlight" -> False|>|>, "DynamicHighlight"], AspectRatio -> GoldenRatio^(-1), Axes -> {True, True}, AxesLabel -> {None, None}, AxesOrigin -> {0, 2.220446049245711*^-16}, DisplayFunction -> Identity, Frame -> {{False, False}, {False, False}}, FrameLabel -> {{None, None}, {None, None}}, FrameTicks -> {{Automatic, Automatic}, {Automatic, Automatic}}, GridLines -> {None, None}, GridLinesStyle -> Directive[GrayLevel[0.5, 0.4]], ImagePadding -> All, Method -> {"DefaultBoundaryStyle" -> Automatic, "DefaultGraphicsInteraction" -> {"Version" -> 1.2, "TrackMousePosition" -> {True, False}, "Effects" -> {"Highlight" -> {"ratio" -> 2}, "HighlightPoint" -> {"ratio" -> 2}, "Droplines" -> {"freeformCursorMode" -> True, "placement" -> {"x" -> "All", "y" -> "None"}}}}, "DefaultMeshStyle" -> AbsolutePointSize[6], "ScalingFunctions" -> None, "CoordinatesToolOptions" -> {"DisplayFunction" -> ({(Identity[#]& )[Part[#, 1]], (Identity[#]& )[Part[#, 2]]}& ), "CopiedValueFunction" -> ({(Identity[#]& )[Part[#, 1]], (Identity[#]& )[Part[#, 2]]}& )}}, PlotRange -> {{Rational[-1, 512], Rational[1, 512]}, {-4.141932917378821*^-28, 4.228768671840614*^-28}}, PlotRangeClipping -> True, PlotRangePadding -> {{Scaled[0.02], Scaled[0.02]}, {Scaled[0.05], Scaled[0.05]}}, Ticks -> {Automatic, Automatic}],Graphics[Annotation[{{{{}, {}, Annotation[{Directive[Opacity[1.], RGBColor[0.368417, 0.506779, 0.709798], AbsoluteThickness[2]], Line[CompressedData["*)
(*1:eJwB4QQe+yFib1JlAgAAAE0AAAACAAAALj8W6v//T79IxKY+IpvUOW/5I4j5*)
(*+k+/akGyzqbf1jmwszEm8/VPvwlWDRyUnoS5MihNYubrT7+KuEb8RyZ8OTcR*)
(*hNrM10+/IfhmUXuVyDlA4/HKma9Pv75fuxLXhNM5UofNqzNfT78Rwng19IOx*)
(*uXfPhG1nvk6/vBymxpqOrjnU1AvqwWFNv2A0JwdMzry5HcBvFDccTL+QsC/H*)
(*mGayubVrxAUP3Uq/FFsmt7rLu7lgjDHw2YJJvxv9Yvuj98S5+JJ7iL8/SL+J*)
(*feduB5SxuaMO3hmY4Ua/3Y4X+JqquLmcSjFy04lFvzOa4st2wb+5gmxheClJ*)
(*RL+f1FARcY+tuXsDqndy7UK/zt5V3P3PoblhgM8k1qhBv4D3jWDhIqG5lb3l*)
(*mJxqQL+HHjbNTgR7ubrfKAysIj6/dpeHI82rZLkjEEBCVJ47v1L85L3/MIy5*)
(*siqIauLjOL/eAmOpBnhEud/FsSA2Nja/8ZS8sgDEcbnlLJUyv7YzvzqcDbzR*)
(*zTI5EX6pNi4BMb92RpJB9Rk5uSw27yyl8yy/SxPS1j8bQ7mDRO3QuXgnvx2Y*)
(*I5LQvCQ5FVSukFkXIr8ylkMLG88KubL2xQ/IJBq/QOcKYFbS+jil2+WN6YQO*)
(*v6MgKgLHGMC4XVEctNtj9L5/NOrRMM9JuLkKe9XCdfM+jvOPQ0lPL7g14gEe*)
(*gVgPPyVCReiVp5M4I/C7+TrCGT8fi66Y0yjhOEmj2c1OWCI/vZDlxr4F+DhF*)
(*TRID9bUnPzBHgUbENzK5j1/XgDC3LD/RsF+eTa8/ucdOHQ1QEjE/82NcpP7g*)
(*MTntIRV+0pozP8IjOfXvcze5dXQrYY8WNj9U1frjJr1CudfcEFJmyDg/bOYM*)
(*YANsaLlgeTznB0w7P/X/BhPVG2e5wys3isMFPj9IrGJqpBOYuSYJvOikSEA/*)
(*bEm1RbrldLkdvGtFBYhBP9+E5orfX5y5APoCqXLiQj/NK1L/wcewufdRvV7F*)
(*JUQ/UPKXZ9TXtLnbNF8bJYRFPxMfLfz+Joq5cFcQESLcRj8mZn8ljDGzuRmU*)
(*5FgEHUg/hcCwVuddvrmuW6Cn83hJPzqpjHjo6bG5Vz1/SMi9Sj/i7eyJHsa2*)
(*ubFebSI6/Es/qP8oRihvuLn4CkMDuVVNP2mgvLezksS5UtE7Nh2YTj+aHb++*)
(*1BjIOQo7C8G8nU4/qmMPZhK2tDnBpNpLXKNOP/Ifu2K4LLM5MHh5YZuuTj+r*)
(*Mxr0vrLIuQ4ft4wZxU4/hnZzC6SZyTnJbDLjFfJOP7DCSKyOQaA5QAgpkA5M*)
(*Tz9e9ue6HLXAOfhx+BquUU8/USU6dfyb0Tmv28elTVdPP0Gq/oWIIZC5Hq9m*)
(*u4xiTz+qW2VJkqWhOfxVpOYKeU8/DK1FuE8UqLm3ox89B6ZPP/o19IFQUMC5*)
(*bg3vx6arTz8rVlmG/XC8OSZ3vlJGsU8/o2inWWZNqrmUSl1ohbxPP98AnPJt*)
(*hMo5cvGakwPTTz+OGCYe/fyhOSpbah6j2E8/lNUw3qecwTnhxDmpQt5PP26n*)
(*HkT8g6i5UJjYvoHpTz9NWf5x1Y7BOQgCqEkh708//42NxMcT1jm/a3fUwPRP*)
(*P/gTaxidH7W5dtVGX2D6Tz+Qm14cWeK3OS4/Fur//08/SMSmPiKb1DlwAmHN*)
(**)
(*"]]}, "Charting`Private`Tag#1"]}}, {}}, <|"HighlightElements" -> <|"Label" -> {"XYLabel"}, "Ball" -> {"InterpolatedBall"}|>, "LayoutOptions" -> <|"PanelPlotLayout" -> <||>, "PlotRange" -> {{Rational[-1, 1024], Rational[1, 1024]}, {-2.4354331366140135`*^-30, 4.511029897537687*^-30}}, "Frame" -> {{False, False}, {False, False}}, "AxesOrigin" -> {0, 2.22044604925029*^-16}, "ImageSize" -> {360, 360/GoldenRatio}, "Axes" -> {True, True}, "LabelStyle" -> {}, "AspectRatio" -> GoldenRatio^(-1), "DefaultStyle" -> {Directive[Opacity[1.], RGBColor[0.368417, 0.506779, 0.709798], AbsoluteThickness[2]]}, "HighlightLabelingFunctions" -> <|"CoordinatesToolOptions" -> ({Identity[Part[#, 1]], Identity[Part[#, 2]]}& ), "ScalingFunctions" -> {{Identity, Identity}, {Identity, Identity}}|>, "Primitives" -> {}, "GCFlag" -> False|>, "Meta" -> <|"DefaultHighlight" -> {"Dynamic", None}, "Index" -> {}, "Function" -> Plot, "GroupHighlight" -> False|>|>, "DynamicHighlight"], AspectRatio -> GoldenRatio^(-1), Axes -> {True, True}, AxesLabel -> {None, None}, AxesOrigin -> {0, 2.22044604925029*^-16}, DisplayFunction -> Identity, Frame -> {{False, False}, {False, False}}, FrameLabel -> {{None, None}, {None, None}}, FrameTicks -> {{Automatic, Automatic}, {Automatic, Automatic}}, GridLines -> {None, None}, GridLinesStyle -> Directive[GrayLevel[0.5, 0.4]], ImagePadding -> All, Method -> {"DefaultBoundaryStyle" -> Automatic, "DefaultGraphicsInteraction" -> {"Version" -> 1.2, "TrackMousePosition" -> {True, False}, "Effects" -> {"Highlight" -> {"ratio" -> 2}, "HighlightPoint" -> {"ratio" -> 2}, "Droplines" -> {"freeformCursorMode" -> True, "placement" -> {"x" -> "All", "y" -> "None"}}}}, "DefaultMeshStyle" -> AbsolutePointSize[6], "ScalingFunctions" -> None, "CoordinatesToolOptions" -> {"DisplayFunction" -> ({(Identity[#]& )[Part[#, 1]], (Identity[#]& )[Part[#, 2]]}& ), "CopiedValueFunction" -> ({(Identity[#]& )[Part[#, 1]], (Identity[#]& )[Part[#, 2]]}& )}}, PlotRange -> {{Rational[-1, 1024], Rational[1, 1024]}, {-2.4354331366140135`*^-30, 4.511029897537687*^-30}}, PlotRangeClipping -> True, PlotRangePadding -> {{Scaled[0.02], Scaled[0.02]}, {Scaled[0.05], Scaled[0.05]}}, Ticks -> {Automatic, Automatic}]}*)


(* ::Input:: *)
(*Map[*)
(*FindMaximum[{Cos[x]-1+x^2/2-machineEvaluate[#[[2]][x]],x>=-#[[1]]&&x<=#[[1]]},{x,9#[[1]]/10},WorkingPrecision->30]&,*)
(*cosGeneralApproximations2]*)


(* ::Input:: *)
(*N[Log[2,Max[Abs[Table[Cos[x]-1+x^2/2-machineEvaluate[cosGeneralApproximations2[[4,2]][x]],{x,9 cosGeneralApproximations2[[4,1]]/10,cosGeneralApproximations2[[4,1]],cosGeneralApproximations2[[4,1]]/1*^6}]]]],30]*)


(* ::Input:: *)
(*CorrectlyRound[CoefficientList[cosGeneralApproximations2[[4,2]][x],x]]*)


(* ::Input:: *)
(*HexLiteral[CorrectlyRound[CoefficientList[cosGeneralApproximations2[[4,2]][x],x]]]*)


(* ::Input:: *)
(*Map[Log[2,N[Abs[(Cos[x]-1+x^2/2-machineEvaluate[#[[2]][x]]/.x->#[[1]])],20]]&,cosGeneralApproximations2]*)


(* ::Input:: *)
(*Map[N[Cos[x]-1-machineEvaluate[#[[2]][x]]/.x->#[[1]],20]&,cosGeneralApproximations2]*)


(* ::Input:: *)
(*Map[Log[2,N[Abs[(Cos[x]-1+x^2/2-machineEvaluate[#[[2]][x]]/.x->#[[1]])],20]]&,cosGeneralApproximations2]*)


(* ::Text:: *)
(*A radius of 1/128 gives us 77 bits, which is compatible with what we have for the Sin function.*)


(* ::Section:: *)
(*Cutoff Points*)


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
