# task-graph-simulation

We provide here additional material for the paper "A general approach to scheduling and checkpointing workflows" by Li Han, Valentin Le Fèvre, Louis-Claude Canon, Yves Robert and Frédéric Vivien.

The code available is the code used for all simulations presented in the paper. The Makefile should compile it into an executable named "simu", so compile by using

make

You can run simulations using

./simu input.csv ckpt_strategy horizon

Input.csv is the input file needed to run the simulations. It has to be formatted as a csv file with the following structure:
one line "parent,child,weight,files" at the beginning
for each dependence between two tasks, a line of the form:
Task1,Task2,Val,{'file1': val1, 'file2': val2, ..., 'filen': valn}
where Task1 is the name of the task that sends data to Task2, Val is a (now) useless value, and then the description of the files with the time needed to transfer them.
Then a line "task_id,weight,proc,ckpt_s1,ckpt_s2,ckpt_all,ckpt_none" followed by a line for each task:
Task,W,P,c1,c2,c3,c4
where Task is the name of Task (can be anything except pure integers), W is the duration of the task, P is the processor on which it is mapped, c[1-4] is a boolean (1 or 0) indicating whether the task is checkpointed for each of the 4 strategies available. c3 should be 1 and c4 should be 0 to respect "ckpt_all" and "ckpt_none" strategies.
Finally, one line for each processor of the form:
i,task1,task2,...,taskn
where i is the number of the processor, and task1 to taskn the names of the tasks mapped on it. The order is important as it represents the list scheduling.
Input.csv must have a special format too: [Name]_[NbTasks]_numP[NbProcs]_lam[ErrorRate]_[anything].csv

ckpt_strategy is the choice of the strategy: ckpt_s1, ckpt_s2, ckpt_all or ckpt_none (as defined in the input file). every is also an acceptable keyword and runs the 4 strategies, while almost runs everything except ckpt_none (for time concerns).

horizon is the maximum time of the simulator. Choosing a big horizon increases the time spent generating the errors, while choosing a too small horizon will prevent the simulator from finishing the execution and will return -1 as makespan.
