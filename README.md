# BTTP

This code solves the bipartite traveling tournament problem (BTTP).

Our source code is written by c++. When the number of teams in each league satisfies $n = 16$, it can genetate a good schedule for BTTP.

① The coordinates (latitude, longitude) of home venues of 32 teams (including 16 teams in the Western Conference and 16 teams in the Eastern Conference) are provided in locations.txt.

② The 3-cycle construction with the swapping heuristics can be found in BTTP_3_cycle.cpp, and the 3-path construction with the swapping heuristics can be found in BTTP_3_path.cpp.

-- The input has already been given, and the output includes the best found result generated by our algorithm with the a running time.

-- The 16 teams in the Western Conference are labeled as {1,2,...,16}, and the 16 teams in the Eastern Conference are labeled as {17,18,...,32}.

-- The output schedule is displayed by a $2n\times 2n$ matrix, where the $i$-th row indicates team $t_i$, the $j$-th column indicates the $j$-th day in the schedule, the item $t_{i,j}$ (resp., $-t_{i,j}$) on the $i$-th row and $j$-th column indicates that team $t_i$ plays an away (resp., home) game with team $t_{i,j}$.

③ The brute-and-force algorithm for calculating the independent lower bound can be found in GetILB.cpp.

-- The input has already been given, and the output is the value of the independent lower bound. It may take half an hour to calculate the independent lower bound.
