# include "armadillo.hpp"
# include <vector>
# include <algorithm>
# include <random>

arma::ucube start_lattice(int size, std::mt19937 &randomElement);
// initialize the lattice. with the bottom of the lattice filled with "bacteria"

arma::ucube add_sphere(arma::ucube lattice, int x, int y, int z, int r);
// Add a sphere of radius r and center (x,y,z) at where growth is impossible.

std::vector<int> aggregate(arma::ucube lattice);
// Vector containg the idices of the initial aggregate of bectrias of type 1 and 2

std::vector<int> surface(arma::ucube lattice, std::vector<int> aggregate);
// vector containing indexes of the initial surface

int random_surface_element(std::vector<int> surface, std::mt19937 &randomElement);
// Draw a random element from the surface

int random_growth_element(arma::ucube lattice, int element, std::mt19937 &randomElement);
// Select a random element for the model to grow into, given a random surface element

arma::ucube Growth_step(arma::ucube lattice, int growth_element);
// Grow into the element chosen by random growth element.

std::vector<int> surface_update(arma::ucube lattice, std::vector<int> surface, int growth_element);

void Eden_Growth(arma::ucube &Colony, int n_steps, std::mt19937 &randomElement); 
// Preform "Eden steps" until the minimum of the colony suface is greater than the value "cutoff".

arma::ucube Prevented_Growth(arma::ucube lattice);
// A lattice containing only the cube where growth is impossible, for visualization purposes.

arma::ucube Colony_type1(arma::ucube lattice);
// lattice containing only the colny of type 1

arma::ucube Colony_type2(arma::ucube lattice);
// lattice containing only the colny of type 2

arma::ucube All_Colonies(arma::ucube lattice);

void print(std::vector<int> const &input);
// Print a std vector


void Save_lattice_as_pgm(arma::ucube lattice, int step_number);

int main(int argc, char const *argv[]){
	std::mt19937 randomElement(522130432);
	arma::ucube Colony;

	//Chose a lattice size and number of steps
	//Simulations of colonies of size order 10^6 similar size to the in vitro colonies
	//May many hours to run
	//Avoid setting steps to large compared to the lattice size as the model may grow out of bounds
	int lattice_size = 20;
	int Growth_steps = 100;


	Colony = start_lattice(lattice_size, randomElement);
	Eden_Growth(Colony, Growth_steps, randomElement);

	// Save the colony as a raw file, that can be visualise in fx. Paraview or similar software
	Colony.save("TwoTypes.raw", arma::raw_binary);

	// Prints the lattice directly to the terminal
	Colony.print();

	}

arma::ucube start_lattice(int size, std::mt19937 &randomElement){
	// Initiate an empty lattice of the specified size
	arma::ucube lattice(size, size, size, arma::fill::zeros);

	
	// Add seed of a colony with radius r, where bacteria of type 1 and 2 are placed randomly within a sphere of radius 3
	int r = 3;
	int x = lattice.n_cols/2;
	int y = lattice.n_rows/2;
	int z = lattice.n_slices/2;
	std::uniform_int_distribution<int> randN(1,2);

	for (int i = r-x-1; i < r+x+1; ++i)
	{
		for (int j = r-y-1; j < r+y+1; ++j)
		{
			for (int k = r-z-1; k < r+z+1; ++k)
			{
				if (sqrt(pow(i-x,2)+pow(j-y,2)+pow(k-z,2)) <= r)
				{
					int random = randN(randomElement);
					lattice(i,j,k) = random;
				}
			}
		}
	}




	return lattice;
}


// Function that finds all sites in the lattice occupied by bacteria
std::vector<int> aggregate(arma::ucube lattice){
	arma::uvec agg1 = find(lattice == 1);
	arma::uvec agg2 = find(lattice == 2);
	arma::uvec agg = join_cols(agg1,agg2);
	std::vector<int> aggregate = arma::conv_to< std::vector<int> >::from(agg);
	return aggregate;
}

// Fuction that finds all surface sites of the colony ie. all the sites that can grow in the following step
std::vector<int> surface(arma::ucube lattice, std::vector<int> aggregate){
	std::vector<int> s;
	for (int i = 0; i < aggregate.size(); ++i)
	{
		if (lattice(aggregate.at(i)+1)==0)
		{
			s.push_back(aggregate.at(i));
		}
		else if (lattice(aggregate.at(i)-1)==0)
		{
			s.push_back(aggregate.at(i));
		}
		else if (lattice(aggregate.at(i)+lattice.n_rows)==0)
		{
			s.push_back(aggregate.at(i));
		}
		else if (lattice(aggregate.at(i)-lattice.n_rows)==0)
		{
			s.push_back(aggregate.at(i));
		}
		else if (lattice(aggregate.at(i)+lattice.n_rows*lattice.n_cols)==0)
		{
			s.push_back(aggregate.at(i));
		}
		else if (lattice(aggregate.at(i)-lattice.n_rows*lattice.n_cols)==0)
		{
			s.push_back(aggregate.at(i));
		}
	}
	std::random_shuffle(s.begin(), s.end());

	return s;
}


// Selects a random element from a vector containing the indicies of all surface elements
int random_surface_element(std::vector<int> surface, std::mt19937 &randomElement){
	std::uniform_int_distribution<int> randN(0,surface.size()-1);
	int random = randN(randomElement);
	int j = surface.at(random);
	return j;
}


// Given a surface element select randomly one of the possible neighbour sites to grow to
int random_growth_element(arma::ucube lattice, int element, std::mt19937 &randomElement){
	std::vector<int> O;
	if(lattice(element+1)==0){
		O.push_back(element+1);
	}
	if(lattice(element-1)==0){
		O.push_back(element-1);
	}
	if (lattice(element+lattice.n_rows)==0)
	{
		O.push_back(element+lattice.n_rows);
	}
	if (lattice(element-lattice.n_rows)==0)
	{
		O.push_back(element-lattice.n_rows);
	}
	if (lattice(element+lattice.n_rows*lattice.n_cols)==0)
	{
		O.push_back(element+lattice.n_rows*lattice.n_cols);
	}
	if (lattice(element-lattice.n_rows*lattice.n_cols)==0)
	{
		O.push_back(element-lattice.n_rows*lattice.n_cols);
	}
	std::uniform_int_distribution<int> rand_N(0,O.size()-1);
	int random = rand_N(randomElement);
	int j = O.at(random);
	//if(lattice(element)==1){
		//lattice(j) = 1;
	//}
	//else if(lattice(element)==2){
		//lattice(j) = 2;
	//}
	return j;
}

// Fill the selected site with a bacetria if the parent bacteria is type 1 fill with 1 and vice versa
arma::ucube Growth_step(arma::ucube lattice, int element, int growth_element){
	if(lattice(element)==1){
		lattice(growth_element) = 1;
	}
	else if(lattice(element)==2){
		lattice(growth_element) = 2;
	}
	return lattice;
}

// Update the surface vector by removing the surface site selected in the previous steps if it no longer has any empty neighbours
// Add the newly grown to site if it has any empty neighbours
// More efficient than redefining the aggregate and surface each step
std::vector<int> surface_update(arma::ucube lattice, std::vector<int> surface, int growth_element){
	std::vector<int> s = surface;
	if (lattice(growth_element+1)==0)
	{
		s.push_back(growth_element);
	}
	else if (lattice(growth_element-1)==0)
	{
		s.push_back(growth_element);
	}
	else if (lattice(growth_element+lattice.n_rows)==0)
	{
		s.push_back(growth_element);
	}
	else if (lattice(growth_element-lattice.n_rows)==0)
	{
		s.push_back(growth_element);
	}
	else if (lattice(growth_element+lattice.n_rows*lattice.n_cols)==0)
	{
		s.push_back(growth_element);
	}
	else if (lattice(growth_element-lattice.n_rows*lattice.n_cols)==0)
	{
		s.push_back(growth_element);
	}

	std::vector<int> n;
	n.push_back(growth_element+1);
	n.push_back(growth_element-1);
	n.push_back(growth_element+lattice.n_rows);
	n.push_back(growth_element-lattice.n_rows);
	n.push_back(growth_element+lattice.n_rows*lattice.n_cols);
	n.push_back(growth_element-lattice.n_rows*lattice.n_cols);

	std::vector<int>::iterator it;
	for (int i = 0; i < n.size(); ++i)
	{
		it = std::find(s.begin(), s.end(), n.at(i));
			if (it != s.end())
		{
			int index = std::distance(s.begin(), it);
			s.erase(s.begin()+index);
		}
	}

	for (int i = 0; i < n.size(); ++i)
	{
		if(lattice(n.at(i))==1 || lattice(n.at(i))==2){
			if (lattice(n.at(i)+1)==0)
			{
				s.push_back(n.at(i));
			}
			else if (lattice(n.at(i)-1)==0)
			{
				s.push_back(n.at(i));
			}
			else if (lattice(n.at(i)+lattice.n_rows)==0)
			{
				s.push_back(n.at(i));
			}
			else if (lattice(n.at(i)-lattice.n_rows)==0)
			{
				s.push_back(n.at(i));
			}
			else if (lattice(n.at(i)+lattice.n_rows*lattice.n_cols)==0)
			{
				s.push_back(n.at(i));
			}
			else if (lattice(n.at(i)-lattice.n_rows*lattice.n_cols)==0)
			{
				s.push_back(n.at(i));
			}

		}
		
	}
	std::random_shuffle(s.begin(), s.end());
	return s;
}


// Inplements the Eden Growth algorithm
void Eden_Growth(arma::ucube &Colony, int n_steps, std::mt19937 &randomElement){
	int counter = 0;
	std::vector<int> a = aggregate(Colony);
	std::vector<int> s = surface(Colony, a);

	std::vector<double> time;
	double t = 0.0;
	time.push_back(t);
	int int_time = 0;

	while(counter < n_steps){
		int element = random_surface_element(s, randomElement);
		int element2 = random_growth_element(Colony, element, randomElement);
		Colony = Growth_step(Colony, element, element2);
		s = surface_update(Colony,s, element2);
		t+= 1/(double)(s.size());
		//std::cout << t << "\n";
		//if(t - time.back() >= 1){
			//int_time++;
			//time.push_back(t);
			//Save_lattice_as_pgm(Colony, int_time);
			//std::cout << t << "\n";
		//}
		counter++;
		
		
	}
}


// Finds all sites of type 1
arma::ucube Colony_type1(arma::ucube lattice){
	lattice.elem( find(lattice != 1) ).zeros();
	return lattice;
}

// Finds all sites of type 2
arma::ucube Colony_type2(arma::ucube lattice){
	lattice.elem( find(lattice != 2) ).zeros();
	return lattice;
}

void print(std::vector<int> const &input){
	for (int i = 0; i < input.size(); ++i)
	{
		std::cout << input.at(i) << " ";
	}
	std::cout << "\n";
}

// Function that saves an instance of the lattice at a given time
void Save_lattice_as_pgm(arma::ucube lattice, int time_stamp){
	std::string sn = std::to_string(time_stamp);
	std::string time_string = "t" + sn;

	for (int i = 0; i < lattice.n_slices; ++i)
	{
		std::string str = std::to_string(i);
		lattice.slice(i).save(time_string + "_"+ str+ ".pgm", arma::pgm_binary);
	}
}
