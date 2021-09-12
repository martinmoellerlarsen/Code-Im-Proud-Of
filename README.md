# Code-Im-Proud-Of
Eden Growth Code for PhD fellow application
The following code is a model of bacterial growth on a 3D lattice called Eden Growth. The method places a new bacteria at the edge of the colony at each step
thereby slowly growing the colony outwards from an initial seed in the middle of the lattice. In this version we simulate a colony where two species of bacteria 
is competing for space within the same colony. The model randomly selects a site on the surface of the colony at each step and the from the set of empty sites adjecent to the surface sites places a new bacteria of the same type as its parent site. This process is then iterated over a specified number of times

The lattice is a 3D cube with sites containing a 0 are defined as empty, sites containg a 1 is occupied by a bacteria of type 1 and sites occupied by a 2 is
a bacteria of type 2

Realistic colonies that corresponds in size to the bacterial colonies i grew experimentally are of an order of magnitude 10^6 bacteria
This may take many hours to simulate, the code is set to a small colony of approx 100 bacteria for testing and can be either printed directly to the terminal
or saved as a .raw file to visualize using appropriate software

Included in the folder is a small gif ColonyGrowth.gif that shows a colony growing up to 2*10^6 bacteria
Also an image of the final colony ColonyImage.png is included
