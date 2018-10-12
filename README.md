# Parallel_MPI_2Dsorting
Thread level implementation of 2D sorting in MPI

The sorting criteria will change slightly and the treatment of holes will change.

You have to do 2D sorting using multiple nodes — one process per node, each process should use openMP to maximally utilize all allocated cores on that node. The 2D grid of now contains only keys (and not string value). However, key is now a two-tuple <int, float> comprising an int and a floating point number.

The goal is to sort them so that each row is sorted from left to right and each column is sorted from top to bottom. The algorithm you should follow is:

Repeat until each row and column are sorted but no more than 4 times:

    Sort each row in increasing order
    
    Sort each column in increasing order

The file format is explicitly sparse. This time there is a single input file listing a sequence of <r, c, i, f> 4-tuples, in no particular order. r is the row number and c is the column number for key <i, f>. An <r, c> pair is never repeated. There is no compression of holes this time, i.e., a hole remains a hole.
(To be clear, the format is native binary as before.)

Sorting order is as follows:

    Row order is given by:
           min(<i, f>,  <j, g>)  is min(f, g)
                      if (f == g) the one in the smaller column is smaller.
    Column order is given by:
           min(<i, f>,  <j, g>)  is min(i, j)
                      if (i == j) the one in the smaller row is smaller.

The maximum number of rows and columns will be specified on the command line as follows (to be included along with mpirun)

sort2d2 <inputfile> <outputfile> <maxrows> <maxcolumns>

The output is in a single file in the same format as the input, except the keys are listed finally in row major order.
