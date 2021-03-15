Complete the table below with your results, and then provide your interpretation at the end.

Note that:

- When calculating the parallel speed-up S, use the time output by the code, which corresponds
  to the parallel calculation and does not include reading in the file or performing the serial check.

- Take as the serial execution time the time output by the code when run with a single process
  (hence the speed-up for 1 process must be 1.0, as already filled in the table).


No. Process:                        Mean time (average of 3 runs)           Parallel speed-up, S:
===========                         ============================:           ====================
1                                   0.002342                                1.0
2                                   0.001464                                1.59972
4                                   0.001762                                1.32917

Architecture that the timing runs were performed on:
Mac OS big Sur with a 64-bit x86 CPU. 

A brief interpretation of these results (2-3 sentences should be enough):
Running with two processes shows a significant increase in performance, likely the overhead of communication is offset due to the number of data points. 
Using four processes still shows a relevant increase in speed over serial performance, but has suffered in comparison to the two process run; likely due
to the increase in communication needed. If ran on a network of machines there will also liklely be an even larger slow down when running with a high
number of processes. 