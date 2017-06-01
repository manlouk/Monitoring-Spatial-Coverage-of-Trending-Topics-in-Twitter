#include<string>
#include<vector>
#include "Set.cpp"
using namespace std;

class Cluster {
		
	private:

		int population;
		Set cluster_elements;
		void increase_population(int);
		
	public:
		Cluster();		
		
		Cluster(Set);
		
		Cluster(Set,int);
		
		int get_cluster_population();
		
		Set get_cluster_elements();

		void set_cluster_elements(Set);
		
		double Jaccard_Index(Set);
	
		double overlap_measure(Set);
		
		void add_cluster_element(Set); 
		
		void merge(Cluster);
		
		void sorting();
		
		char get_first_char();
		
		
		
};				
