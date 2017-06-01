#include<algorithm>
#include<vector>
#include "cluster.h"

using namespace std;

Cluster::Cluster(Set cluster_elements,int population){
	this->cluster_elements=cluster_elements;
	this->population=population;
}
Cluster::Cluster(){
	this->population=0;
}
Cluster::Cluster(Set cluster_elements){
	this->population=1;
	this->cluster_elements=cluster_elements;
}
int Cluster::get_cluster_population(){
	return population;
}


Set Cluster::get_cluster_elements(){
	return cluster_elements;
}
void Cluster::sorting(){
	sort(this->cluster_elements.begin(),this->cluster_elements.end());
}
double Cluster::Jaccard_Index(Set element){
	
	Set A=this->get_cluster_elements();
	
	Set num=set_intersection(A,element);
		
	Set den=set_union(A,element);
	
	return ((num.size())/(double)(den.size()));	
}
double Cluster::overlap_measure(Set element){
	
	Set A=this->get_cluster_elements();
	
	int num=set_intersection(A,element).size();
	int den=min(A.size(),element.size());
	return num/((double)den);
}
void Cluster::increase_population(int N){
	population+=N;
}



void Cluster::set_cluster_elements(Set x){
	this->cluster_elements=x;
}
void Cluster::add_cluster_element(Set element){
	
		Set new_cluster=set_union(this->get_cluster_elements(),element);
		this->set_cluster_elements(new_cluster);		
		this->increase_population(1);	
}
void Cluster::merge(Cluster B)
{
	Set new_cluster=set_union(this->get_cluster_elements(),B.get_cluster_elements());
	this->set_cluster_elements(new_cluster);
	this->increase_population(B.get_cluster_population());
	
}
char Cluster::get_first_char(){
	char c;	
	Set hashtags = this->get_cluster_elements();
		for(int i=0;i<hashtags[0].size();i++){
			if((hashtags[0].at(i)>='a')&&(hashtags[0].at(i)<='z')){	
				c=hashtags[0].at(i);
				break;
			}
	}
	
	return c;
}
	
