# Hochbaum's Pseudoflow (HPF) Algorithm for (2-Piecewise Linear) Fully Parametric Minimum Cut — v3
This package provides a fully parametric implementation of Hochbaum's Pseudoflow algorithm for minimum cut on directed graphs. In the parametric minimum cut problem, the capacity of source-adjacent arcs is monotone non-decreasing in the parameter `lambda` whereas the capacity of sink-adjacent arcs is monotone non-increasing in `lambda`. This solver requires that the capacities of source and sink adjacent arcs are 2-piecewise linear in `lambda`: `capacity = max(0, constant + multiplier * lambda)`.

This fully parametric solver finds the optimal minimum cut for all `lambda` values in a given range. The solution for all lambda values is represented with `O(n)` intervals for the parameter lambda. In each interval, the optimal minimum cut remains the same.

Version 3 is a faster rewrite of [pseudoflow-parametric-cut-v2](https://github.com/hochbaumGroup/pseudoflow-parametric-cut-v2). It finds the same breakpoints and is **8–40× faster** on the benchmark graphs below, with the gap growing with graph size.

A simple parametric minimum cut solver that provides the optimal minimum cut for a given list of arc capacities is available [here](https://riot.ieor.berkeley.edu/Applications/Pseudoflow/parametric.html), and a non-parametric maximum flow version of pseudoflow is available [here](https://riot.ieor.berkeley.edu/Applications/Pseudoflow/maxflow.html).

This implementation uses a variant of the fully parametric HPF algorithm as described in:
>    DS Hochbaum (2008), The Pseudoflow algorithm: A new algorithm for the maximum flow problem. Operations Research, 58(4):992-1009.

## What is new in v3
The algorithm is the same as in v2: the range `[lower bound, upper bound]` is split recursively, deciding for each interval whether it contains 0, 1, or more breakpoints. What changes is the cost of each step.

- **Recursive contraction.** Minimum cuts are nested in `lambda`, so nodes already in the source set at the lower end of an interval, or in the sink set at the upper end, are merged into a super-source and a super-sink. Each subproblem is built from its parent's already-contracted graph, so subproblems shrink as the recursion goes deeper. Contracted arcs keep their `constant` and `multiplier`, and parallel arcs are merged correctly under the `max(0, ...)` floor.
- **Cheaper intersection of cut functions.** The next `lambda` to test is computed from the contracted arc list of the current level, restricted to arcs adjacent to nodes whose side differs between the two cuts. v2 scans every arc of the original graph at every recursive call.
- **Specialized intersection for symmetric instances.** When every parametric arc is a source arc `max(0, c + m*lambda)` or a sink arc `max(0, -c - m*lambda)` of the same node, as in densest-subgraph problems, the instance is detected automatically and a faster intersection is used.
- **Cache-friendly data layout.** Nodes and arcs are addressed by `uint` indices instead of pointers.
- **Automatic rounding.** The `max(0, ...)` floor is enabled automatically when some arc capacity becomes negative within the `lambda` range; the input file no longer has a rounding flag.
- **Compact output.** One line per breakpoint, listing the nodes that join the source set at that `lambda`.
- **Pure C.** The MATLAB and Python wrappers of v2 are not included. The default build prints nothing to stdout and compiles without warnings.

Information from earlier solves is reused only structurally, through contraction: each subproblem contains only the nodes not yet decided by the cuts at its interval endpoints. Pseudoflow state (flows and labels) is not carried between solves, and *free runs* (see pg.15 of the paper above) are not used. These remain possible further improvements.

## Instructions for C
From the repository root, compile the `hpf` executable with `make` (or run `make` in `src/pseudoflow/c`).

To execute the solver, use:
```bash
src/pseudoflow/c/hpf input-file.txt output-file.txt
```

The input file should contain the graph structure and is assumed to have the following format:
```
    c <comment lines>
    p <# nodes> <# arcs> <lower bound> <upper bound>
    n <source node> s
    n <sink node> t
    a <from-node> <to-node> <constant capacity> <lambda multiplier>
```
where the `a` line is repeated for each arc. The file should satisfy the following conditions:
- Nodes are labeled `0 .. <# nodes> - 1`.
- `<lambda multiplier>` is non-negative if `<from-node> == <source node>` and `<to-node> != <sink-node>`.
- `<lambda multiplier>` is non-positive if `<from-node> != <source node>` and `<to-node> == <sink-node>`.
- `<lambda multiplier>` is zero if `<from-node> != <source node>` and `<to-node> != <sink-node>`.

Fields after `<upper bound>` on the `p` line are ignored, so v2 input files are read without changes. Note that the rounding flag of v2 is ignored as well: negative capacities are always rounded to zero.

The solver will generate the following output file:
```
t <time (in sec) read data> <time (in sec) initialize> <time (in sec) solve>
s <# arc scans> <# mergers> <# pushes> <# relabels> <# gaps>
p <number of breakpoints = k>
l <lambda> <node> <node> ...
```
The `l` line appears once for each breakpoint, in increasing order of `lambda`, and lists the nodes that are added into the source set at that `lambda` value.

See `src/pseudoflow/c/example` for a small example. To get progress and tracing output, build with `make CFLAGS='-c -fpic -Wall -std=gnu99 -O4 -DVERBOSE'` in `src/pseudoflow/c`.

## Samples
The `samples` directory contains four densest-subgraph instances with their output. In each instance, every vertex `v` of the input graph has a source arc of capacity `max(0, lambda - deg(v))` and a sink arc of capacity `max(0, deg(v) - lambda)`, every edge has unit-capacity arcs, and `lambda` ranges over `[0, 1000]`.

| Instance | Wall time | Peak memory |
|---|---|---|
| `close-cliques` | 0.13 s | 23 MB |
| `ego-facebook` | 0.14 s | 26 MB |
| `email-Enron` | 0.55 s | 71 MB |
| `soc-Slashdot0811` | 1.38 s | 175 MB |

Times are for one core of an Intel Xeon Gold 6330 (Savio cluster, UC Berkeley), with gcc 8.5.

## Benchmarks
Solve time (third field of the `t` line) of v2 and v3 on densest-subgraph instances built from SNAP graphs. The full table is in `benchmarks/v2_vs_v3.tsv`, and `benchmarks/v3_all_runs.tsv` has v3 results on 36 graphs, up to `wikidata` (2113 s, 6231 breakpoints).

| Graph | v2 (s) | v3 (s) | Speedup |
|---|---|---|---|
| ego-twitter | 65.29 | 1.61 | 40.6× |
| com-dblp | 142.78 | 3.73 | 38.3× |
| com-amazon | 169.18 | 5.33 | 31.8× |
| soc-Epinions1 | 10.89 | 0.42 | 25.7× |
| soc-Slashdot0811 | 11.80 | 0.47 | 24.9× |
| soc-Slashdot0902 | 14.14 | 0.59 | 24.1× |
| email-Enron | 3.75 | 0.21 | 17.8× |
| ego-facebook | 0.43 | 0.05 | 8.1× |

The breakpoints found by v3 (both the `lambda` values and the node sets) are identical to those of v2 with its correctness fixes applied. Without those fixes, v2 reports one breakpoint fewer on each of these instances, at the low end of the `lambda` range.

## Citation
When using this code, please cite:
> DS Hochbaum (2008), The Pseudoflow algorithm: A new algorithm for the maximum flow problem. Operations Research, 58(4):992-1009.

> DS Hochbaum and B Chandran (2009), A Computational Study of the Pseudoflow and Push-relabel Algorithms for the Maximum Flow Problem. Operations Research, 57(2):358-376.

> Q Spaen, B Fishbain and DS Hochbaum, Hochbaum's Pseudo-flow C Implementation, http://riot.ieor.berkeley.edu/riot/Applications/Pseudoflow/maxflow.html

## License
See `LICENSE.md`. Free for educational, research, and not-for-profit purposes; contact UC Berkeley's Office of Technology Licensing for commercial use.
