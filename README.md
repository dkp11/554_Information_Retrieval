*********************************************************************************************
**************		   554 - Information Retrieval			*********************
**************			   Homework 2				*********************
**************		PageRank algorithm Implementation		*********************
*********************************************************************************************
Github link:
https://github.com/dkp11/554_Information_Retrieval/
********************************************************************************************
Programming Language: 	R
Script File Name: 	PageRank.R
No of Implementation:	2

*********************************************************************************************
* Implementation
*********************************************************************************************
The general implementation rule is same or common across methods. The implementation follows
in line with the algorithm.
1. Take the input adjacency matrix
2. Create a graph using igraph library
3. Get the adjacency matrix from the graph
4. Do necessary preprocessing to the adjacency matrix
5. Create a unit vector and then a bookmark vector
6. Calculate the page rank for all the pages using the pageRank formula
7. Get the eigen vectors
8. Sort the vectors descendingly to get top 100 nodes
9. Get their corresponding values
*********************************************************************************************

*********************************************************************************************
Common steps:
=============
The steps listed here are common to all methods and hence should be followed before 
going on to individual methods:
1. Open PageRank.R script file in a R environment. I used R version 3.1.2
2. Select lines 1 through 24 and run line by line for mentioned behaviour
3. These lines does the following:
	- load the 'igraph' library
	- load the input data from the "AM.txt" file to a matrix (the file name should be AM.tx)
	- Create a graph out of the matrix (inputGraph)
	- Creates an adjacency matrix out of the graph (G)
	- declares dampening factor constant (d = 0.85)
4. From here on, we can go for any of the methods to calculate the PageRank.
*********************************************************************************************

*********************************************************************************************
Method 1:
=========
Uses the power method or the Eigen vector function in R. A much simplified implementation. 
Bigger memory and much time is a constraint in case of larger matrices.

1. Choose which method to use - Eigen or power. Default is eigen.
2. Run from lines 29 - 69, line by line.
2. Will take some time (for larger graphs alone) but will give the result vector.
*********************************************************************************************

*********************************************************************************************
Method 2:
=========
This method uses inbuilt library and produces results instanteously. This is used for verification purposes alone.

1. Run line 79.
*********************************************************************************************

*********************************************************************************************
Important Note:
===============
The program was tested in the following environment:
64-bit Windows Server environment with 64/128 GB (RAM) virtual memory.
Testing the same in a 64-bit Windows 8.1 with 8 GB (RAM) failed due to insufficient memory.

Also running times of the methods vary depending on few mins to several hours based on 
- method selected
- resource available
- input file size

Also to note, if a subset of nodes were selected, then the page ranking will be different 
and a different top 100 values and nodes were returned. 
*********************************************************************************************

*********************************************************************************************
Results:
========
Input 1:- Six node graph
========================

The ordering was same in all methods. But the numerical values are different.

Nodes in prder:
---------------
	5 	    2 		3 	  4 	     1 		6

Values in the above order of nodes:
-----------------------------------
Method 1 - Eigen:
-----------------
 0.51474883   0.50288776   0.47359977   0.47359977   0.17342812   0.05887897

Method 1 - Power:
-----------------
 0.0012636881 0.0012345696 0.0011626687 0.0011626687 0.0004257592 0.0001445456

Method 2:
---------
 0.03694043   0.23496947   0.21050244   0.21050244   0.22279983   0.08428539 
 
Input 2:- 81433 node graph
==========================

The results obtained in both the methods were same for the order of nodes but different in values:

Top 100 nodes based on the PageRank value (Reads from left to right, top to bottom: 1 to 100)
---------------------------------------------------------------------------------------------
  [1]  1804 30383  4510  4470  7633  1065  1579    91  4468  4291  4621  8263  4550
 [14]  4521  4091  3276  7893  4826  4469   190   176  4263  3300  3673   710 13217
 [27]  4513 11299 24274   255   440   371 21365  2524  8825    96  9024  4467  1052
 [40] 14007 11102 19328  3365  4508  4575  2110  4586  3129  4772 12764  1395  8082
 [53] 18815 12132  8510   342  2878  4528  5274 13997 15640 12402 10549   659  4671
 [66]  4594  4477 12784  4504 34798  8478 12327 15913   192  9095  3112  4805 20970
 [79]   628   257  9083 13684   239 37834 14473   902   793   702  4544 11742  9529
 [92]  3206  4541  1906  9792  3321  7717 21343   158  4733

Top 100 values based on the PageRank algorithm (Reads from left to right, top to bottom: 1 to 100)
--------------------------------------------------------------------------------------------------
 [1] 0.0096407001 0.0093931740 0.0039245391 0.0028897442 0.0027942227 0.0026862274
 [7] 0.0026237740 0.0024911912 0.0024360264 0.0024221044 0.0018796827 0.0018204389
 [13] 0.0017502665 0.0017013465 0.0016549865 0.0015114594 0.0014668878 0.0014646926
 [19] 0.0014603717 0.0014272081 0.0013640730 0.0013493886 0.0013378711 0.0012967479
 [25] 0.0012172393 0.0011875079 0.0011743047 0.0011695446 0.0011645818 0.0011584006
 [31] 0.0011513301 0.0011136162 0.0010691942 0.0010406083 0.0009849123 0.0009818481
 [37] 0.0009797909 0.0009622440 0.0009466164 0.0009290173 0.0009202798 0.0009112039
 [43] 0.0008936824 0.0008702739 0.0008671125 0.0008560494 0.0008474022 0.0008173247
 [49] 0.0008139462 0.0008063718 0.0008032386 0.0008019976 0.0007995251 0.0007924632
 [55] 0.0007828019 0.0007804547 0.0007740613 0.0007512787 0.0007443707 0.0007435256
 [61] 0.0007320001 0.0007229492 0.0007204874 0.0007199263 0.0007144972 0.0007136729
 [67] 0.0007116203 0.0007056514 0.0006981857 0.0006823621 0.0006806075 0.0006799051
 [73] 0.0006593956 0.0006542951 0.0006471924 0.0006450648 0.0006444335 0.0006401155
 [79] 0.0006331471 0.0006210566 0.0006167713 0.0006061622 0.0006027354 0.0006014844
 [85] 0.0005901081 0.0005883785 0.0005878007 0.0005848247 0.0005846955 0.0005725401
 [91] 0.0005678998 0.0005620070 0.0005543582 0.0005527643 0.0005511434 0.0005366670
 [97] 0.0005362828 0.0005318937 0.0005298201 0.0005279586
*********************************************************************************************

*********************************************************************************************
**************			  End of Homework 2			*********************
**************			   Happy Grading!!!			*********************
*********************************************************************************************
