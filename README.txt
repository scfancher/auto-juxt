DphaseplotLL.cpp is our code for calculating the relative error in the center cell of a cluster of cells communicating via autocrine or juxtacrine signaling. There are 8 parameters listed at the top of the code that define the space over which the code looks.
1) sn is the starting number of cells.
2) N is the largest number of cells. Together, sn and N define the range of population sizes over which the code calculates.
3) nn is the number of population sizes the code looks at. These sizes will start at sn and grow exponentially until they reach the maximum of N so as to make the range of population sizes linear on a logarithmic scale (if the exponential growth is too small, linear growth will be used to insure population sizes are integers).
4) sl is the starting lambda/a value.
5) fl is the largest lambda/a value. Together, sl and fl define the range of lambda/a values over which the code calculates.
6) nl is the number of lambda/a values the code looks at. Similarly to nn, these values will start at sl and grow exponentially until they reach a maximum of fl (again, if exponential growth is too small, linear growth will be used to insure lambda/a values that are not too similar).
7) D defines the diffusion rate of the cells in the monty carlo portion of the autocrine system. This dictates how far each cell moves during each step of the monty carlo simulation.
8) bl defines the number of consecutive moves that must be rejected due to their causing an increase in the relative error before the code stops the monty carlo loop. Increasing bl increases the probability that the final configuration exists at a minimum.

The output file, NLdata.txt, has 5 columns of data.
1) The first column gives the lambda/a value for that row.
2) The second column gives the population size for that row.
3) The third column gives the average nearest-neighbor distance of cells in the optimized autocrine system.
4) The fourth column gives the relative error of the center cell of the juxtacrine system.
5) The fifth column gives the relative error of the center cell of the optimized autocrine system.