SoPlex version 2.0.0 [mode: optimized] [precision: 8 byte] [rational: gmp] [githash: 568f354]
Copyright (c) 1996-2014 Konrad-Zuse-Zentrum fuer Informationstechnik Berlin (ZIB)

int:simplifier = 0
int:readmode = 1
int:solvemode = 2
real:feastol = 1e-20
real:opttol = 1e-20

Reading (rational) LP file <output/lp_file_gene_deletion_Formate_CO_1.lp> . . .
Reading took 0.46 seconds.

LP has 4206 rows 5427 columns and 89797 nonzeros.

Added 1 slack columns to transform rows to equality form.
Initial floating-point solve . . .
Equilibrium scaling LP

type |   time |   iters | facts |  shift   |    value
  L  |    0.0 |       0 |     1 | 0.00e+00 | 0.00000000e+00
  L  |    0.0 |     200 |     2 | 0.00e+00 | 0.00000000e+00
  L  |    0.0 |     400 |     4 | 0.00e+00 | 0.00000000e+00
  L  |    0.1 |     600 |     6 | 0.00e+00 | 0.00000000e+00
  L  |    0.1 |     800 |     7 | 0.00e+00 | 1.70256423e+00
  L  |    0.1 |    1000 |     8 | 0.00e+00 | 1.75354660e+00
  L  |    0.1 |    1200 |     9 | 0.00e+00 | 1.79193723e+00
  L  |    0.1 |    1400 |    10 | 0.00e+00 | 2.47190620e+00
  L  |    0.2 |    1600 |    11 | 0.00e+00 | 2.47190643e+00
  L  |    0.2 |    1800 |    12 | 0.00e+00 | 2.47190643e+00
  L  |    0.3 |    2000 |    13 | 0.00e+00 | 2.47190643e+00
  L  |    0.3 |    2200 |    14 | 0.00e+00 | 2.47244041e+00
  L  |    0.3 |    2400 |    16 | 6.15e-05 | 8.26743527e+00
  L  |    0.4 |    2600 |    17 | 1.37e-04 | 8.27588902e+00
  L  |    0.4 |    2800 |    19 | 1.37e-04 | 8.27588902e+00
  L  |    0.5 |    3000 |    21 | 1.37e-04 | 8.27588902e+00
  L  |    0.6 |    3200 |    22 | 1.98e-03 | 8.43785206e+00
  L  |    0.6 |    3200 |    22 | 1.98e-03 | 8.43785206e+00
  L  |    0.6 |    3400 |    23 | 3.34e-03 | 8.67960127e+00
  L  |    0.6 |    3400 |    23 | 3.34e-03 | 8.67960127e+00
  L  |    0.8 |    3600 |    24 | 3.70e-03 | 8.68314815e+00
  E  |    1.1 |    3800 |    67 | 3.26e-07 | 8.68425080e+00
  E  |    1.1 |    3800 |    67 | 3.26e-07 | 8.68425080e+00
  E  |    1.1 |    3800 |    67 | 3.26e-07 | 8.68425080e+00
  E  |    1.2 |    3800 |    69 | 3.26e-07 | 8.68425080e+00
  E  |    1.2 |    3800 |    69 | 3.26e-07 | 8.68425080e+00

Floating-point optimal.
Max. bound violation = 1.7579765700197535707e-07
Max. row violation = 2.7868536021459301056e-07
Max. dual violation = 9.5875343149468946006e-07
Refined floating-point solve . . .

type |   time |   iters | facts |  shift   |    value
  E  |    0.0 |       0 |     1 | 4.31e+00 | 6.47349698e-05
  L  |    0.1 |     200 |     2 | 3.65e-11 | 1.16484120e+10
  L  |    0.2 |     400 |     3 | 3.65e-11 | 1.18934069e+10
  L  |    0.3 |     600 |     4 | 3.65e-11 | 1.23695942e+10
  L  |    0.5 |     800 |     8 | 5.39e-02 | 1.54006201e+13
  L  |    0.5 |     800 |     8 | 5.39e-02 | 1.54006201e+13
  L  |    0.5 |     800 |     8 | 5.39e-02 | 1.54006201e+13
  L  |    0.5 |     800 |     8 | 5.39e-02 | 1.54006201e+13
  L  |    0.5 |     800 |     8 | 5.39e-02 | 1.54006201e+13
  L  |    0.5 |     800 |     8 | 5.39e-02 | 1.54006201e+13
  L  |    0.5 |     800 |     8 | 5.39e-02 | 1.54006201e+13
  L  |    0.5 |     800 |     8 | 5.39e-02 | 1.54006201e+13
  L  |    0.5 |     800 |     8 | 5.39e-02 | 1.54006201e+13
  L  |    0.8 |    1000 |    16 | 0.00e+00 | 1.15607409e+10
  L  |    0.9 |    1200 |    17 | 0.00e+00 | 1.18538974e+10
  L  |    0.9 |    1400 |    18 | 0.00e+00 | 1.21909404e+10
  L  |    1.1 |    1600 |    21 | 0.00e+00 | 1.86806150e+11

Numerical troubles during floating-point solve.
Increasing Markowitz threshold.
Solving from scratch.
Equilibrium scaling LP

type |   time |   iters | facts |  shift   |    value
  L  |    0.0 |       0 |     3 | 2.03e-02 | 1.61883075e-02
type |   time |   iters | facts |  shift   |    value
  L  |    0.0 |       0 |     3 | 2.03e-02 | 1.61883075e-02
  L  |    0.0 |     200 |     4 | 2.03e-02 | -7.01864988e-01
  L  |    0.0 |     400 |     5 | 2.03e-02 | -5.61144391e-03
  L  |    0.0 |     600 |     6 | 2.03e-02 | -4.08300452e-03
  L  |    0.1 |     800 |     7 | 2.03e-02 | -3.98821723e-03
  L  |    0.1 |    1000 |     8 | 2.03e-02 | -3.92592204e-03
  L  |    0.1 |    1200 |    10 | 2.03e-02 | -3.68489840e-03
  L  |    0.2 |    1400 |    11 | 2.03e-02 | 2.43262964e-02
  L  |    0.2 |    1600 |    12 | 2.03e-02 | 8.04042963e-02
  L  |    0.3 |    1800 |    13 | 2.04e-02 | 1.43146283e-01
  L  |    0.3 |    2000 |    14 | 2.04e-02 | 3.94188353e-01
  L  |    0.4 |    2200 |    15 | 2.06e-02 | 2.59463616e+00
  L  |    0.4 |    2400 |    16 | 2.10e-02 | 2.32293176e+00
  L  |    0.6 |    2600 |    18 | 1.07e-01 | 5.47538553e+07
  L  |    0.6 |    2600 |    18 | 1.07e-01 | 5.47538553e+07
  L  |    0.6 |    2600 |    18 | 1.07e-01 | 5.47538553e+07
  L  |    0.6 |    2600 |    18 | 1.07e-01 | 5.47538553e+07
  L  |    0.6 |    2600 |    18 | 1.07e-01 | 5.47538553e+07
  L  |    0.6 |    2600 |    18 | 1.07e-01 | 5.47538553e+07
  L  |    0.6 |    2600 |    18 | 1.07e-01 | 5.47538553e+07
  L  |    0.6 |    2600 |    18 | 1.07e-01 | 5.47538553e+07
  L  |    0.6 |    2600 |    18 | 1.07e-01 | 5.47538553e+07
  L  |    0.6 |    2600 |    18 | 1.07e-01 | 5.47538553e+07
  L  |    0.6 |    2600 |    18 | 1.07e-01 | 5.47538553e+07
  L  |    0.6 |    2600 |    18 | 1.07e-01 | 5.47538553e+07
  L  |    0.8 |    2800 |    20 | 4.16e+00 | 2.41957918e+09
  L  |    0.8 |    2800 |    20 | 4.16e+00 | 2.41957918e+09
  L  |    0.8 |    2800 |    20 | 4.16e+00 | 2.41957918e+09
  L  |    0.8 |    2800 |    20 | 4.16e+00 | 2.41957918e+09
  L  |    0.8 |    2800 |    20 | 4.16e+00 | 2.41957918e+09
  L  |    0.8 |    2800 |    20 | 4.16e+00 | 2.41957918e+09
  L  |    0.8 |    2800 |    20 | 4.16e+00 | 2.41957918e+09
  L  |    0.8 |    2800 |    20 | 4.16e+00 | 2.41957918e+09
  L  |    0.8 |    2800 |    20 | 4.16e+00 | 2.41957918e+09
  L  |    0.8 |    2800 |    20 | 4.16e+00 | 2.41957918e+09
  L  |    0.8 |    2800 |    20 | 4.16e+00 | 2.41957918e+09
  L  |    0.8 |    2800 |    20 | 4.16e+00 | 2.41957918e+09
  L  |    0.8 |    2800 |    20 | 4.16e+00 | 2.41957918e+09
  L  |    0.8 |    2800 |    20 | 4.16e+00 | 2.41957918e+09
  L  |    0.8 |    2800 |    20 | 4.16e+00 | 2.41957918e+09
  L  |    0.8 |    2800 |    20 | 4.16e+00 | 2.41957918e+09
  L  |    0.8 |    2800 |    20 | 4.16e+00 | 2.41957918e+09
  L  |    0.8 |    2800 |    20 | 4.16e+00 | 2.41957918e+09
  L  |    0.8 |    2800 |    20 | 4.16e+00 | 2.41957918e+09
  L  |    0.8 |    2800 |    20 | 4.16e+00 | 2.41957918e+09
  L  |    0.8 |    2800 |    20 | 4.16e+00 | 2.41957918e+09
  L  |    0.8 |    2800 |    20 | 4.16e+00 | 2.41957918e+09
  L  |    1.0 |    3000 |    26 | 1.89e+01 | 8.71453808e+09
  L  |    1.1 |    3200 |    30 | 2.42e+01 | 1.08916126e+10
  L  |    1.2 |    3400 |    31 | 2.47e+01 | 1.09059637e+10
  L  |    1.2 |    3400 |    31 | 2.47e+01 | 1.09059637e+10
  L  |    1.3 |    3600 |    33 | 2.47e+01 | 1.12650668e+10
  L  |    1.4 |    3800 |    34 | 2.47e+01 | 1.14432225e+10
  L  |    1.5 |    4000 |    35 | 2.47e+01 | 1.17456042e+10
  L  |    1.8 |    4200 |    45 | 2.86e+01 | 1.47961070e+11
  L  |    1.8 |    4200 |    45 | 2.86e+01 | 1.47961070e+11
  L  |    1.8 |    4200 |    45 | 2.86e+01 | 1.47961070e+11
  L  |    1.8 |    4200 |    45 | 2.86e+01 | 1.47961070e+11
  L  |    1.9 |    4400 |    49 | 2.03e-02 | -1.32235310e+01
  L  |    1.9 |    4600 |    50 | 2.03e-02 | -1.17735945e+01
  L  |    1.9 |    4800 |    51 | 2.03e-02 | -8.28540093e-01
  L  |    1.9 |    5000 |    52 | 2.03e-02 | -7.73761221e-01
  L  |    1.9 |    5200 |    53 | 2.03e-02 | -6.12239746e-01
  L  |    2.0 |    5400 |    54 | 2.03e-02 | -5.99610701e-01
  L  |    2.0 |    5600 |    56 | 2.03e-02 | -5.53683506e-01
  L  |    2.0 |    5800 |    57 | 2.03e-02 | -6.99699280e-01
type |   time |   iters | facts |  shift   |    value
  L  |    2.1 |    6000 |    58 | 2.04e-02 | -6.96045062e-01
  L  |    2.1 |    6200 |    59 | 2.04e-02 | -6.61097177e-01
  L  |    2.1 |    6400 |    60 | 2.04e-02 | -2.06236206e-01
  L  |    2.2 |    6600 |    63 | 2.31e-02 | 6.07512531e+07
  L  |    2.3 |    6800 |    64 | 2.31e-02 | 1.06713828e+10
  L  |    2.4 |    7000 |    65 | 1.79e-01 | 1.07313684e+10
  L  |    2.5 |    7200 |    66 | 2.29e-01 | 1.09758690e+10
  L  |    2.6 |    7400 |    67 | 2.29e-01 | 1.09852526e+10
  L  |    2.6 |    7400 |    67 | 2.29e-01 | 1.09852526e+10
  L  |    2.7 |    7600 |    69 | 2.97e-01 | 1.11372365e+10
  L  |    2.7 |    7600 |    69 | 2.97e-01 | 1.11372365e+10
  L  |    2.7 |    7600 |    69 | 2.97e-01 | 1.11372365e+10
  L  |    2.7 |    7600 |    69 | 2.97e-01 | 1.11372365e+10
  L  |    2.7 |    7600 |    69 | 2.97e-01 | 1.11372365e+10
  L  |    2.8 |    7800 |    73 | 5.49e+00 | 1.23303664e+10
  L  |    2.9 |    8000 |    74 | 5.49e+00 | 1.25565868e+10
  L  |    3.0 |    8200 |    77 | 5.55e+00 | 1.29364738e+10
  L  |    3.3 |    8400 |    81 | 4.46e+01 | 4.77397129e+10
  L  |    3.3 |    8400 |    81 | 4.46e+01 | 4.77397129e+10
  L  |    3.6 |    8600 |    90 | 4.59e+01 | 1.12668690e+11
  L  |    3.9 |    8800 |    98 | 8.45e+01 | 3.02163040e+11

Switching scaling.

type |   time |   iters | facts |  shift   |    value
  L  |    0.0 |       0 |     1 | 1.83e+00 | 1.43625283e-02
  L  |    0.0 |     200 |     2 | 1.83e+00 | -2.04836265e+01
  L  |    0.1 |     400 |     3 | 1.83e+00 | -2.55897218e+02
  L  |    0.1 |     600 |     4 | 1.83e+00 | -1.84682598e+01
  L  |    0.1 |     800 |     5 | 1.83e+00 | -1.70198481e+01
  L  |    0.1 |    1000 |     6 | 1.83e+00 | -1.66878732e+01
  L  |    0.2 |    1200 |     7 | 1.83e+00 | -1.63470100e+01
  L  |    0.2 |    1400 |     8 | 1.83e+00 | -1.32452097e+01
  L  |    0.2 |    1600 |     9 | 1.83e+00 | -1.28559598e+01
  L  |    0.3 |    1800 |    10 | 1.83e+00 | -1.28547163e+01
  L  |    0.3 |    2000 |    12 | 1.83e+00 | -1.28522881e+01
  L  |    0.4 |    2200 |    13 | 1.84e+00 | -1.28469403e+01
  L  |    0.4 |    2400 |    14 | 1.84e+00 | -1.28460752e+01
  L  |    0.4 |    2600 |    15 | 1.84e+00 | -1.28461133e+01
  L  |    0.5 |    2800 |    16 | 1.84e+00 | 1.03465604e+10
  L  |    0.5 |    2800 |    16 | 1.84e+00 | 1.03465604e+10
  L  |    0.6 |    3000 |    17 | 7.65e+00 | 1.08607062e+10
  L  |    0.6 |    3000 |    17 | 7.65e+00 | 1.08607062e+10
  L  |    0.7 |    3200 |    18 | 7.65e+00 | 1.15163761e+10
  L  |    0.7 |    3400 |    19 | 7.66e+00 | 1.17989476e+10
  L  |    0.8 |    3600 |    20 | 7.67e+00 | 1.22273634e+10
  L  |    0.8 |    3600 |    20 | 7.67e+00 | 1.22273634e+10
  L  |    0.8 |    3600 |    20 | 7.67e+00 | 1.22273634e+10
  L  |    0.8 |    3600 |    20 | 7.67e+00 | 1.22273634e+10
  L  |    0.8 |    3600 |    20 | 7.67e+00 | 1.22273634e+10
  L  |    0.8 |    3600 |    20 | 7.67e+00 | 1.22273634e+10
  L  |    0.8 |    3600 |    20 | 7.67e+00 | 1.22273634e+10
  L  |    0.8 |    3600 |    20 | 7.67e+00 | 1.22273634e+10
  L  |    0.8 |    3600 |    20 | 7.67e+00 | 1.22273634e+10
  L  |    0.8 |    3600 |    20 | 7.67e+00 | 1.22273634e+10
  L  |    0.8 |    3600 |    20 | 7.67e+00 | 1.22273634e+10
  L  |    0.8 |    3600 |    20 | 7.67e+00 | 1.22273634e+10
  L  |    1.0 |    3800 |    21 | 7.99e+00 | 4.98514396e+10
  L  |    1.2 |    4000 |    22 | 8.03e+00 | 1.85602717e+11
  L  |    1.5 |    4200 |    28 | 8.06e+00 | 1.85923688e+11

Floating-point infeasible.
Setting up LP to test for feasibility.
Initial floating-point solve . . .
Equilibrium scaling LP

type |   time |   iters | facts |  shift   |    value
  L  |    0.0 |       0 |     1 | 0.00e+00 | -1.00000000e+00

Floating-point optimal.
Max. bound violation = 0
Max. row violation = 0
Max. dual violation = 15.439681955265582403
Refined floating-point solve . . .

type |   time |   iters | facts |  shift   |    value
  E  |    0.0 |       0 |     1 | 0.00e+00 | 0.00000000e+00

Floating-point optimal.
Max. bound violation = 0
Max. row violation = 0
Max. dual violation = 2.9029940296466697725e-07
Refined floating-point solve . . .

type |   time |   iters | facts |  shift   |    value
  E  |    0.0 |       0 |     1 | 0.00e+00 | 0.00000000e+00

Floating-point optimal.
Max. bound violation = 0
Max. row violation = 0
Max. dual violation = 1.0128207880890229749e-09
Refined floating-point solve . . .

type |   time |   iters | facts |  shift   |    value
  E  |    0.0 |       0 |     1 | 0.00e+00 | 0.00000000e+00

Floating-point optimal.
Max. bound violation = 0
Max. row violation = 0
Max. dual violation = 1.1027586355023549162e-20
Refined floating-point solve . . .

type |   time |   iters | facts |  shift   |    value
  E  |    0.0 |       0 |     1 | 0.00e+00 | 0.00000000e+00

Floating-point optimal.
Max. bound violation = 0
Max. row violation = 0
Max. dual violation = 3.2776763210371532444e-30
Tolerances reached.
Primal infeasible.  Dual Farkas ray available.

SoPlex status       : problem is solved [infeasible]
Solving time (sec)  : 9.94
Iterations          : 18799
Objective value     : 1.00000000e+100

No primal solution available.No primal solution available.
No dual solution available.

Statistics
==========

SoPlex status       : problem is solved [infeasible]
Original problem    : 
  Columns           : 5427
              boxed : 36
        lower bound : 5427
        upper bound : 36
               free : 0
  Rows              : 4206
             ranged : 4205
                lhs : 4206
                rhs : 4205
               free : 0
  Nonzeros          : 89797
         per column : 16.55
            per row : 21.35
           sparsity : 0.00
    min. abs. value : 0.00
    max. abs. value : 3440491270.24
Objective sense     : minimize
Solution (rational) : 
  Objective value   : 1.0000000000000000159e+100
Violation (rational): 
  Max/sum bound     : - / -
  Max/sum row       : - / -
  Max/sum redcost   : - / -
  Max/sum dual      : - / -
Total time          : 10.39
  Reading           : 0.45
  Solving           : 9.94
  Preprocessing     : 0.01 (0.10% of solving time)
  Simplex           : 8.77 (88.23% of solving time)
  Synchronization   : 0.08 (0.80% of solving time)
  Transformation    : 0.03 (0.30% of solving time)
  Other             : 1.05 (10.56% of solving time)
Refinements         : 5
  Stalling          : 3
Iterations          : 18799
  From scratch      : 17150 (91.23%)
  From basis        : 1649 (8.77%)
  Primal            : 201 (1.07%)
  Dual              : 18598 (98.93%)
  Bound flips       : 0
LU factorizations   : 265
  Factor. frequency : 70.94 iterations per factorization
  Factor. time      : 1.88
LU solves           : 50193
  Solve frequency   : 2.67 solves per iteration
  Solve time        : 2.71
