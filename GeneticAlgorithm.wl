(* ::Package:: *)

BeginPackage["GeneticAlgorithm`"]


GAMaximize::usage = "GAMaximize[func, vars] maximizes the given function func numerically with the Genetic Algorithm. Variables should be written in form of {\!\(\*SubscriptBox[\(x\), \(1\)]\)\[Element]Reals, \!\(\*SubscriptBox[\(x\), \(2\)]\)\[Element]Integers}. 
References: 
[1] K.Deep,M.Thakur,A new crossover operator for real coded genetic algorithms,Applied Mathematics and Computation.188 (2007) 895\[Dash]911. https://doi.org/10.1016/j.amc.2006.10.047.
[2] K.Deep,M.Thakur,A new mutation operator for real coded genetic algorithms,Applied Mathematics and Computation.193 (2007) 211\[Dash]230. https://doi.org/10.1016/j.amc.2007.03.046.
[3] K. Deep, K.P. Singh, M.L. Kansal, C. Mohan, A real coded genetic algorithm for solving integer and mixed integer optimization problems, Applied Mathematics and Computation. 212 (2009) 505\[Dash]518. https://doi.org/10.1016/j.amc.2009.02.044.
[4] R. Hinterding, Gaussian mutation and self-adaption for numeric genetic algorithms, in: Proceedings of 1995 IEEE International Conference on Evolutionary Computation, IEEE, Perth, WA, Australia, 1995: p. 384. https://doi.org/10.1109/ICEC.1995.489178.
Possible options includes:
PrecisionGoal: integer, default value 5
AccuracyGoal: integer, default value 5
MaxIterations: integer, default value 10000
PopulationSize: integer, default value 10 * degree of freedem
Bounds: default value None
ShrinkingFactor: real, default value 1, must be smaller than 1
Scale: real, default value: 1
TournamentSize: integer, default value 3
MutationProbability: real, default value 0.01
CrossoverProbability: real, default value 0.8
MaxUnchangedGeneration: integer, default value MaxIterations/4
bReal: real, default value 0.15
bInt: real, default value 0.35
a: real, default value 0
pReal: real, default value 10
pInt: real, default value 4
"
PopulationSize::usage = "An option for GAMaximize for specifying the population size of a generation in the GA. Default 30."
ShrinkingFactor::usage = "An option for GAMaximize for the shrinking factor of the standard deviation during the Gaussian mutation. Default 1."
TournamentSize::usage = "An option for GAMaximize for specifying the tournament size during the tournament selection in the GA. Default 3."
MutationProbability::usage = "An option for GAMaximize for specifying the mutation probability in the GA. Default 0.005."
CrossoverProbability::usage = "An option for GAMaximize for specifying the crossover probability in the GA. Default 0.8."
MaxUnchangedGeneration::usage = "An option for GAMaximize for specifying the maximum allowed generation number for the known maximum value. Default half of MaxIterations."

GAMinimize::usage = "GAMinimize[func, vars] minimizes the given function func numerically with the Genetic Algorithm. Variables should be written in form of {\!\(\*SubscriptBox[\(x\), \(1\)]\)\[Element]Reals, \!\(\*SubscriptBox[\(x\), \(2\)]\)\[Element]Integers}. 
References: 
[1] K.Deep,M.Thakur,A new crossover operator for real coded genetic algorithms,Applied Mathematics and Computation.188 (2007) 895\[Dash]911. https://doi.org/10.1016/j.amc.2006.10.047.
[2] K.Deep,M.Thakur,A new mutation operator for real coded genetic algorithms,Applied Mathematics and Computation.193 (2007) 211\[Dash]230. https://doi.org/10.1016/j.amc.2007.03.046.
[3] K. Deep, K.P. Singh, M.L. Kansal, C. Mohan, A real coded genetic algorithm for solving integer and mixed integer optimization problems, Applied Mathematics and Computation. 212 (2009) 505\[Dash]518. https://doi.org/10.1016/j.amc.2009.02.044.
[4] R. Hinterding, Gaussian mutation and self-adaption for numeric genetic algorithms, in: Proceedings of 1995 IEEE International Conference on Evolutionary Computation, IEEE, Perth, WA, Australia, 1995: p. 384. https://doi.org/10.1109/ICEC.1995.489178.
Possible options includes:
PrecisionGoal: integer, default value 5
AccuracyGoal: integer, default value 5
MaxIterations: integer, default value 10000
PopulationSize: integer, default value 10 * degree of freedem
Bounds: default value None
ShrinkingFactor: real, default value 1, must be smaller than 1
Scale: real, default value: 1
TournamentSize: integer, default value 3
MutationProbability: real, default value 0.01
CrossoverProbability: real, default value 0.8
MaxUnchangedGeneration: integer, default value MaxIterations/4
bReal: real, default value 0.15
bInt: real, default value 0.35
a: real, default value 0
pReal: real, default value 10
pInt: real, default value 4
"
PopulationSize::usage = "An option for GAMaximize and GAMinimize for specifying the population size of a generation in the GA. Default 30."
ShrinkingFactor::usage = "An option for GAMaximize and GAMinimize for the shrinking factor of the standard deviation during the Gaussian mutation. Default 1."
TournamentSize::usage = "An option for GAMaximize and GAMinimize for specifying the tournament size during the tournament selection in the GA. Default 3."
MutationProbability::usage = "An option for GAMaximize and GAMinimize for specifying the mutation probability in the GA. Default 0.005."
CrossoverProbability::usage = "An option for GAMaximize and GAMinimize for specifying the crossover probability in the GA. Default 0.8."
MaxUnchangedGeneration::usage = "An option for GAMaximize and GAMinimize for specifying the maximum allowed generation number for the known maximum/minimum value. Default half of MaxIterations."


Begin["Private`"]


$bReal=0.15;
$bInt=0.35;
$a=0.;
$pc=0.8;
$pm=0.005;
$pReal=10;
$pInt=4;
$k=3;(* tournament size *)
$sigma=1;(* initial standard deviation for Gaussian mutation *)
$shrink=1;(* shrink factor for standard deviation in Gaussian mutation *)


(* Laplace crossover, reference: K.Deep,M.Thakur,A new crossover operator for real coded genetic algorithms,Applied Mathematics and Computation.188 (2007) 895\[Dash]911. https://doi.org/10.1016/j.amc.2006.10.047. *)
cLaplaceCrossover=Compile[{
{x,_Real,2},
{types,_Integer,1},(* 1 for integer, 0 for real *)
{pc,_Real},(* crossover probability *)
{a,_Real},(* coefficient a for laplace crossover *)
{bReal,_Real},(* coefficient bReal for laplace crossover *)
{bInt,_Real}},(* coefficient bInt for laplace crossover *)
Module[{dof=Length[x[[1]]],popsize=Length[x],roulette,matingPool=Most@{x[[1]]},matingPoolSize,partedMatingPool,len,u,r,bli,beta,population=x,x1,x2,y},
(* generate random numbers and add chromosomes into mating pool *)
roulette=RandomReal[{0,1},popsize];
Table[If[roulette[[i]]<=pc,AppendTo[matingPool,x[[i]]]],{i,popsize}];
matingPoolSize=Length[matingPool];
(* if there are odd numbers of chromosomes, get rid of the last one and add it to the population (it has offspring with itself) *)
If[OddQ[matingPoolSize],AppendTo[population,matingPool[[-1]]];matingPool=Most[matingPool]];
partedMatingPool=Partition[matingPool,2];
len=Length[partedMatingPool];
(* offspring generation *)
u=RandomReal[{0,1},{len,dof}];
r=RandomReal[{0,1},{len,dof}];
x1=partedMatingPool[[All,1]];
x2=partedMatingPool[[All,2]];
bli=If[#>0,bInt,bReal]&/@types;
beta=a+Map[If[#<=1/2,-1,1]&,r,{2}] Table[bli,len] Map[Log,u,{2}];
Join[population,x1+beta Map[Abs,x1-x2,{2}],x2+beta Map[Abs,x1-x2,{2}]]
],CompilationTarget->"C",Parallelization->True,RuntimeAttributes->{Listable},RuntimeOptions->"Speed"];


(* Power mutation, reference: K.Deep,M.Thakur,A new mutation operator for real coded genetic algorithms,Applied Mathematics and Computation.193 (2007) 211\[Dash]230. https://doi.org/10.1016/j.amc.2007.03.046. *)
cPowerMutation=Compile[{
{x,_Real},
{lower,_Real},
{upper,_Real},
{pm,_Real}(* mutation probability *),
{p,_Real}},(* coefficient p for power mutation *)
Module[{pl,s1,s,r=RandomReal[],t},
pl=RandomReal[{0,1}];
If[pl<pm(* roulette, if win, then mutate, otherwise return the original value *),
s1=RandomReal[];
s=s1^p;
t=(x-lower)/(upper-lower);
x+If[t<r,-s(x-lower),s(upper-x)],x]],CompilationTarget->"C",Parallelization->True,RuntimeAttributes->{Listable},RuntimeOptions->"Speed"];


cGaussianMutation=Compile[{
{x,_Real},
{pm,_Real},(* mutation probability *)
{sigma0,_Real},(* initial standard deviation for Gaussian mutation *)
{shrink,_Real},
{gen,_Integer},
{max,_Integer}},
Module[{pl,sigma=sigma0 (1-shrink gen/max)},
pl=RandomReal[{0,1}];
If[pl<pm(* roulette, if win, then mutate, otherwise return the original value *),
x+RandomVariate[NormalDistribution[0,sigma]],x]],CompilationTarget->"C",Parallelization->True,RuntimeAttributes->{Listable},RuntimeOptions->"Speed"];


cMutation=Compile[{{x,_Real,2},{types,_Integer,1},{bounded,_Integer,1},{lowers,_Real,1},{uppers,_Real,1},{pm,_Real},{pReal,_Real},{pInt,_Real},{sigma0,_Real},{shrink,_Real},{gen,_Integer},{max,_Integer}},
Module[{x0,type,boundQ,dof=Length[x[[1]]],popsize=Length[x],y},
Table[
x0=x[[i,j]];
type=types[[j]];
boundQ=bounded[[j]];
If[boundQ==1,
(* apply power mutation for bounded variables *)
cPowerMutation[x0,lowers[[j]],uppers[[j]],pm,If[type==1,pInt,pReal]],
(* applu gaussian mutation for unbounded variables *)
cGaussianMutation[x0,pm,sigma0,shrink,gen,max]]
,{i,1,popsize},{j,1,dof}]],CompilationTarget->"C",Parallelization->True,RuntimeAttributes->{Listable},RuntimeOptions->"Speed",CompilationOptions->{"InlineExternalDefinitions"->True}];


cTruncation=Compile[{{x,_Real,2},{types,_Integer,1}},
Module[{dof=Length[x[[1]]],popsize=Length[x],x0,type,r},
Table[x0=x[[i,j]];
type=types[[j]];If[type==1,r=RandomReal[];If[r<=0.5,Floor[x0],Ceiling[x0]],x0],{i,popsize},{j,dof}]],CompilationTarget->"C",Parallelization->True,RuntimeAttributes->{Listable},RuntimeOptions->"Speed"];


cApplyBounds=Compile[{{x,_Real,2},{bounded,_Integer,1},{lowers,_Real,1},{uppers,_Real,1}},
Module[{dof=Length[x[[1]]],popsize=Length[x],x0,boundQ},
Table[x0=x[[i,j]];
boundQ=bounded[[j]];If[boundQ==1,Max[Min[x0,uppers[[j]]],lowers[[j]]],x0],{i,popsize},{j,dof}]],CompilationTarget->"C",Parallelization->True,RuntimeAttributes->{Listable},RuntimeOptions->"Speed"];


cGAOperator=Compile[{
{x,_Real,2},
{types,_Integer,1},(* 1 for integer, 0 for real *)
{bounded,_Integer,1},(* 1 for bounded, 0 for unbounded *)
{lowers,_Real,1},
{uppers,_Real,1},
{pc,_Real},(* crossover probability *)
{pm,_Real},(* mutation probability *)
{aC,_Real},(* coefficient a for laplace crossover *)
{bCReal,_Real},(* coefficient bReal for laplace crossover *)
{bCInt,_Real},(* coefficient bInt for laplace crossover *)
{pReal,_Real},(* coefficient pReal for power mutation *)
{pInt,_Real},(* coefficient pInt for power mutation *)
{sigma0,_Real},(* initial standard deviation for Gaussian mutation *)
{shrink,_Real}(* shrinkage coefficient of standard deviation for Gaussian mutation *),
{gen,_Integer}(* current generation number *),
{max,_Integer}(* maximum generation number *)
},
Module[{dof=Length[x[[1]]](* degrees of freedom *),popsize=Length[x](* population size *),population,transposed},
(* crossover *)
population=cLaplaceCrossover[x,types,pc,aC,bCReal,bCInt];
(* mutation *)
population=cMutation[population,types,bounded,lowers,uppers,pm,pReal,pInt,sigma0,shrink,gen,max];
(* truncation *)
population=cTruncation[population,types];
cApplyBounds[population,bounded,lowers,uppers]],CompilationTarget->"C",Parallelization->True,RuntimeAttributes->{Listable},RuntimeOptions->"Speed",CompilationOptions->{"InlineExternalDefinitions"->True}];


(* Tournament selection *)
cTournamentSelection=Compile[{{fitnesses,_Real,1},{x,_Real,2},{populationSize,_Integer},{tournamentSize,_Integer}},
Module[{selection=Most@{x[[1]]},random,indices=Most@{0},length=Length[x],winner,fitnessMax,index},
Table[indices=Most@{0};
While[Length[indices]<tournamentSize,
random=RandomInteger[{1,length}];
If[!MemberQ[indices,random],AppendTo[indices,random]]];
winner=x[[indices[[1]]]];
fitnessMax=fitnesses[[indices[[1]]]];
Table[index=indices[[j]];
If[fitnesses[[index]]>fitnessMax,
fitnessMax=fitnesses[[index]];
winner=x[[index]]],{j,2,tournamentSize}];
AppendTo[selection,winner],{i,populationSize}];
selection],CompilationTarget->"C",Parallelization->True,RuntimeAttributes->{Listable},RuntimeOptions->"Speed"];


Options[GAMinimize]=Options[GAMaximize]={PrecisionGoal->5,AccuracyGoal->5,MaxIterations->100,PopulationSize->Automatic,Bounds->None,ShrinkingFactor->$shrink,Scale->$sigma,TournamentSize->$k,MutationProbability->$pm,CrossoverProbability->$pc,
"bReal"->$bReal,"bInt"->$bInt,"a"->$a,"pReal"->$pReal,"pInt"->$pInt,MaxUnchangedGeneration->Automatic};


GAMaximize[function_,variables_List,OptionsPattern[]]:=
Module[{generation=0,
popsize,
toursize=OptionValue[TournamentSize],
convergeQ=False,
vars=variables[[All,1]],
dof,
types=If[#===Integers,1,0]&/@variables[[All,2]](* 1 for integer, 0 for Real *),
bounded,lowers,uppers,
population,
max,
mean,
knownMax,
lastMaxGen=1,
sol,
knownBest,
lastMean,
aC=N[OptionValue["a"]],
bCInt=N[OptionValue["bInt"]],
bCReal=N[OptionValue["bReal"]],
pReal=N[OptionValue["pReal"]],
pInt=N[OptionValue["pInt"]],
pc=N[OptionValue[CrossoverProbability]],
pm=N[OptionValue[MutationProbability]],
sigma0=N[OptionValue[Scale]],
shrink=N[OptionValue[ShrinkingFactor]],
maxGen=OptionValue[MaxIterations],
fitnesses,
continueQ=True,
acc=OptionValue[AccuracyGoal],
pre=OptionValue[PrecisionGoal],
id,
maxUnchangedGen=If[OptionValue[MaxUnchangedGeneration]===Automatic,Max[100,IntegerPart[OptionValue[MaxIterations]/4]],OptionValue[MaxUnchangedGeneration]]},
dof=Length[vars];
popsize=If[OptionValue[PopulationSize]===Automatic, Max[10 dof, 50], OptionValue[PopulationSize]];
bounded=If[OptionValue[Bounds]===None,ConstantArray[0,dof],If[#===None,0,1]&/@OptionValue[Bounds]];
lowers=If[OptionValue[Bounds]===None,ConstantArray[0.,dof],If[#===None,0.,N[#[[1]]]]&/@OptionValue[Bounds]];
uppers=If[OptionValue[Bounds]===None,ConstantArray[0.,dof],If[#===None,0.,N[#[[2]]]]&/@OptionValue[Bounds]];
(* initialization *)
population=Table[MapThread[If[#1==0,RandomReal[{-100,100}],RandomReal[{#2,#3}]]&,{bounded,lowers,uppers}],{popsize}];
population=cTruncation[population,types];
population=cApplyBounds[population,bounded,lowers,uppers];
max=mean=Apply[function,First[population]];
knownMax=max-100;
sol=First[population];
Monitor[
While[continueQ,
generation++;
If[generation>1,population=cTournamentSelection[fitnesses,population,popsize,toursize]];
population=cGAOperator[population,types,bounded,lowers,uppers,pc,pm,aC,bCReal,bCInt,pReal,pInt,sigma0,shrink,generation,maxGen];
fitnesses=Apply[function,#]&/@population;
max=Max[fitnesses];
{lastMean,mean}={mean,Mean[fitnesses]};
convergeQ=
(If[knownMax>$MachineEpsilon,
Abs[max-knownMax]<=10^(-acc)&&Abs[1-(max/(knownMax+$MachineEpsilon))]< 10^(-pre),
Abs[max-knownMax]<=10^(-acc)]&&
If[Abs[lastMean]>$MachineEpsilon,
(Abs[mean-lastMean]<=10^(-acc)&&Abs[1-(mean/(lastMean+$MachineEpsilon))]< 10^(-pre)),
Abs[mean-lastMean]<=10^(-acc)])||(generation-lastMaxGen>maxUnchangedGen);
If[max>knownMax,
lastMaxGen=generation;
knownMax=max;
sol=population[[Ordering[fitnesses,-1][[1]]]]];
(*Print[knownMax, " ", lastMaxGen, " ", maxUnchangedGen];*)
continueQ=((!convergeQ)&&(generation<maxGen))],
(*Print["Generation ",generation,", Known maximum = ",knownMax]],*)
Row[{"Current generation = " <> ToString[generation]<> "\nKnown maximum = "<>ToString[knownMax] <> "\nCurrent mean = "<>ToString[mean] <> "\nGeneration number for the current maximum  = "<>ToString[lastMaxGen]}]];
If[generation==maxGen&&(!convergeQ),Print["Max iteration steps of " <>ToString[maxGen]<> " reached but the results do not converge!"]];
id=Ordering[fitnesses,-1][[1]];
{knownMax,Thread[vars->sol]}]


GAMinimize[function_,variables_List,OptionsPattern[]]:=
Module[{generation=0,
popsize,
toursize=OptionValue[TournamentSize],
convergeQ=False,
vars=variables[[All,1]],
dof,
types=If[#===Integers,1,0]&/@variables[[All,2]](* 1 for integer, 0 for Real *),
bounded,lowers,uppers,
population,
min,
mean,
knownMin,
lastMaxGen=1,
sol,
knownBest,
lastMean,
aC=N[OptionValue["a"]],
bCInt=N[OptionValue["bInt"]],
bCReal=N[OptionValue["bReal"]],
pReal=N[OptionValue["pReal"]],
pInt=N[OptionValue["pInt"]],
pc=N[OptionValue[CrossoverProbability]],
pm=N[OptionValue[MutationProbability]],
sigma0=N[OptionValue[Scale]],
shrink=N[OptionValue[ShrinkingFactor]],
maxGen=OptionValue[MaxIterations],
fitnesses,
continueQ=True,
acc=OptionValue[AccuracyGoal],
pre=OptionValue[PrecisionGoal],
id,
maxUnchangedGen=If[OptionValue[MaxUnchangedGeneration]===Automatic,Max[100,IntegerPart[OptionValue[MaxIterations]/4]],OptionValue[MaxUnchangedGeneration]]},
dof=Length[vars];
popsize=If[OptionValue[PopulationSize]===Automatic, Max[10 dof, 50], OptionValue[PopulationSize]];
bounded=If[OptionValue[Bounds]===None,ConstantArray[0,dof],If[#===None,0,1]&/@OptionValue[Bounds]];
lowers=If[OptionValue[Bounds]===None,ConstantArray[0.,dof],If[#===None,0.,N[#[[1]]]]&/@OptionValue[Bounds]];
uppers=If[OptionValue[Bounds]===None,ConstantArray[0.,dof],If[#===None,0.,N[#[[2]]]]&/@OptionValue[Bounds]];
(* initialization *)
population=Table[MapThread[If[#1==0,RandomReal[{-100,100}],RandomReal[{#2,#3}]]&,{bounded,lowers,uppers}],{popsize}];
population=cTruncation[population,types];
population=cApplyBounds[population,bounded,lowers,uppers];
min=mean=Apply[function,First[population]];
knownMin=min+100;
sol=First[population];
Monitor[
While[continueQ,
generation++;
If[generation>1,population=cTournamentSelection[-fitnesses,population,popsize,toursize]];
population=cGAOperator[population,types,bounded,lowers,uppers,pc,pm,aC,bCReal,bCInt,pReal,pInt,sigma0,shrink,generation,maxGen];
fitnesses=Apply[function,#]&/@population;
min=Min[fitnesses];
{lastMean,mean}={mean,Mean[fitnesses]};
convergeQ=
(If[knownMin>$MachineEpsilon,
Abs[min-knownMin]<=10^(-acc)&&Abs[1-(min/(knownMin+$MachineEpsilon))]< 10^(-pre),
Abs[min-knownMin]<=10^(-acc)]&&
If[Abs[lastMean]>$MachineEpsilon,
(Abs[mean-lastMean]<=10^(-acc)&&Abs[1-(mean/(lastMean+$MachineEpsilon))]< 10^(-pre)),
Abs[mean-lastMean]<=10^(-acc)])||(generation-lastMaxGen>maxUnchangedGen);
If[min<knownMin,
lastMaxGen=generation;
knownMin=min;
sol=population[[Ordering[fitnesses,1][[1]]]]];
(*Print[knownMax, " ", lastMaxGen, " ", maxUnchangedGen];*)
continueQ=((!convergeQ)&&(generation<maxGen))],
(*Print["Generation ",generation,", Known maximum = ",knownMax]],*)
Row[{"Current generation = " <> ToString[generation]<> "\nKnown minimum = "<>ToString[knownMin] <> "\nCurrent mean = "<>ToString[mean] <> "\nGeneration number for the current minimum  = "<>ToString[lastMaxGen]}]];
If[generation==maxGen&&(!convergeQ),Print["Max iteration steps of " <>ToString[maxGen]<> " reached but the results do not converge!"]];
id=Ordering[fitnesses,-1][[1]];
{knownMin,Thread[vars->sol]}]


End[]


EndPackage[]
