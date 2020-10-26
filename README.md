# Mathematica-GA

A Mathematica implementation of genetic algorithm, which seems to be faster than Mathematica built-in `NMaximize` or `NMinimize` with the option `"DifferentialEvolution"`.

Two main functions in this packages are `GAMaximize` and `GAMinimize`.

`GAMaximize[func, vars]`/`GAMinimize[func, vars]`  maximizes/minimizes the given function func numerically with the Genetic Algorithm. 

Variables should be written in form of `{Element[x, Reals],Element[y, Integers]}`. 

# Example
```
func[x_, y_] := (x + 0.1)^2 + (y + 0.1)^2
GAMinimize[func, {Element[x, Reals], Element[y, Integers]}, PopulationSize -> 100]
```

# References: 

[1] K.Deep,M.Thakur,A new crossover operator for real coded genetic algorithms,Applied Mathematics and Computation.188 (2007) 895\[Dash]911. https://doi.org/10.1016/j.amc.2006.10.047.

[2] K.Deep,M.Thakur,A new mutation operator for real coded genetic algorithms,Applied Mathematics and Computation.193 (2007) 211\[Dash]230. https://doi.org/10.1016/j.amc.2007.03.046.

[3] K. Deep, K.P. Singh, M.L. Kansal, C. Mohan, A real coded genetic algorithm for solving integer and mixed integer optimization problems, Applied Mathematics and Computation. 212 (2009) 505\[Dash]518. https://doi.org/10.1016/j.amc.2009.02.044.

[4] R. Hinterding, Gaussian mutation and self-adaption for numeric genetic algorithms, in: Proceedings of 1995 IEEE International Conference on Evolutionary Computation, IEEE, Perth, WA, Australia, 1995: p. 384. https://doi.org/10.1109/ICEC.1995.489178.


# Options

`PrecisionGoal`: integer, default value 5

`AccuracyGoal`: integer, default value 5

`MaxIterations`: integer, default value 10000

`PopulationSize`: integer, default value 10 * number of variables

`Bounds`: default value None

`ShrinkingFactor`: real, default value 1, must be smaller than 1

`Scale`: real, default value: 1

`TournamentSize`: integer, default value 3

`MutationProbability`: real, default value 0.01

`CrossoverProbability`: real, default value 0.8

`MaxUnchangedGeneration`: integer, default value MaxIterations/4

`bReal`: real, default value 0.15

`bInt`: real, default value 0.35

`a`: real, default value 0

`pReal`: real, default value 10

`pInt`: real, default value 4
