# Vouw
*Pattern Mining in Reduce-Fold Cellular Automata*

Vouw is a command-line interface for printing an analysing Reduce-Fold Cellular Automata (RFCA) and is an academic implementation of the VOUW algorithm. For the visualisation of RFCA, please see the user-friendly RFExplore:
https://github.com/mickymuis/rfexplore

## Builing Vouw

Vouw is coded in C99 and should run on anywhere. I may have used some Linux specific functions calls, so if you run in to any trouble please notify me. 
Compiling is easy, just run in your terminal:
```
git clone https://github.com/mickymuis/vouw.git
cd vouw
./configure
cd build
make
```
The `vouw` executable is now located in `build`.

## Printing automata

Let's print our first automaton! For this we pick automaton 2.2.6 and call the `print` module:
```
./vouw print -m 2 -b 2 -r 6 -i 00110 -f 5
```
If everything is well, this should produce the following output:
```
1 1 0 1 0 0 0 1 1 0 
 0 1 1 1 0 0 1 0 1 
  1 0 0 1 0 1 1 1 
   1 0 1 1 1 0 0 
    1 1 0 0 1 0 
     0 1 0 1 1 
      1 1 1 0 
       0 0 1 
        0 1 
         1 
```

The arguments are the same for every module (the first argument after `./vouw`). Type `./vouw` without arguments to see decriptions for each parameter.

## Analysing with VOUW

We can use the VOUW algorithm to compress an automaton and view the compressed output, code table and compression ratio:
```
mik@trommel:~/code/vouw/build$ ./vouw encode -m 2 -b 2 -r 6 -i 00110 -f 5

```
Which appearently compresses to the following with a ratio of 81%
```
C . B . . A . B . B 
 . . A . B . B . B 
  . B . B . . A C 
   C . . A . . A 
    C . . A . B 
     . B . . A 
      C C . B 
       C . B 
        . B 
         C 

```
Now the real power lies in the compression of an automaton with the code table *of another automaton*. This can be used, for example, as measure of distance. Let's compress a bigger version of 2.2.6 with the code table of 2.2.7:
```
./vouw encode -m 2 -b 2 -r 6 -i 00110 -f 20 using -r 7

```
Notice the extra `using` argument. The ratio is now almost 68%, indicating that 2.2.6 and 2.2.7 are quite similar. If we now were to compress the same automaton with the code table generated from 2.2.8, we would obtain a ratio of 108%, indicating that they are less similar.
