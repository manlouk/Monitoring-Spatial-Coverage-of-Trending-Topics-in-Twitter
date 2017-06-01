//Title: structures.h
//Description: Basic data structures and function for manipulating spatial grid partitioning.
//Platform(s): gcc 4.8.3.
//Date: 29/12/2007
//Revision: 17/3/2015

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <list>
#include <map>
#include <queue>
#include <sys/timeb.h>
#include <ctime>
#include <set>
#include "cluster.cpp"
using namespace std;

typedef vector<string> Set;
//The workspace for all spatial objects (universe) in case this is not deduced from input

#define X_MIN 505000
#define Y_MIN 155000
#define X_MAX 555000
#define Y_MAX 205000


//---------------------------------------------------------------------------
//Get system time in milliseconds
inline unsigned int get_time() {
	timeb t;
	ftime(&t);
	return t.time*1000+t.millitm;
}
vector<string> gethashtags(const string& message){
	int i=0;
	int k,s;
	bool flag;
	vector<string> hashtags;
	while(i<message.size()){
		flag=0;
		if(message[i]=='#'){
			flag=1;
			 k=i;
			 s=i;
			while(((message[s]!=',')&&(message[s]!='\n'))&&(s<message.size())){s++;}
			string hash(message,k,s-k);
			//cout<<hash<<endl;
			hashtags.push_back(hash);	
			
		}
		 i++;
	
	}
	return hashtags;
}

//vector<unsigned int>changed_cells;
//typedef pair<Set,set<unsigned int> > Topics_Cells;
typedef struct{
	Set name;
	set<pair<unsigned int,unsigned int> > cells;
}Topics_Cells;
//topic with dense cell zones in grid
typedef struct Topic_Zones{
	Set topic_name;
	//string label;
	vector<set<unsigned int> > cluster_cells;
	unsigned int freq;
} topic_zones;


/*void printTopicCells(vector<Topics_Cells> A)
{
	//cout<<"*******START WINDOW*********"<<endl;
	for(int i=0;i<A.size();i++){
		Set topic=A[i].first;
		set<unsigned int> v_cells=A[i].second;
		set<unsigned int>::iterator it1,it2;
		
		int k=0;
		for(it1=v_cells.begin();it1!=v_cells.end();++it1){
			for(it2=v_cells.begin();it2!=v_cells.end();++it2){
				if(isNeighboorCell(*it1,*it2,GranX))
					k++;
			}
		}
		cout<<"{";
		for(int c=0;c<topic.size();c++){
			cout<<topic[c]<<" ";
		}
		cout<<"};{";
		for(it1=v_cells.begin();it1!=v_cells.end();++it1){
			
			cout<<*it1<<" ";
			//k++;
		}
		cout<<" };"<<k<<endl;
	}
	//cout<<"*********END OF WINDOW********"<<endl;
}
*/
//A rectangular box defined by its lower left and upper right corner.
typedef struct boxInfo 
{
	double x_min;
	double y_min;
	double x_max;
	double y_max;
        //Constructors for boxes
	boxInfo(const double x1, const double y1, const double x2, const double y2) 
	{
            x_min = x1;
            y_min = y1;
            x_max = x2;
            y_max = y2;
	}
        boxInfo() {};
} box;


//Positional data for a given object
typedef struct Obj_tuple {
        unsigned int ts;	//Location timestamp
	unsigned int id;	//Object id
	Set tags;	
	bool fresh;		//Indicator for newly received location
	double x;               //May be used for storing LONgitude coordinates
	double y;               //May be used for storing LATitude coordinates
} TObject;


class Compare
{
public:
    bool operator() (Cluster c1, Cluster c2)
    {
        return ((c1.get_cluster_elements()).size()<(c2.get_cluster_elements()).size());
    }
};

//Collection of object locations indexed by their ID
typedef map<unsigned int, TObject* > TObjChain;       //CAUTION: pointer reference to object location
typedef vector<Cluster>Clusters_per_cell;
typedef struct Cluster_struct{
	Clusters_per_cell clusters_set;//Contents of cluster
	unsigned int ts;//cluster timestamp...latest time a tweet inserted in cluster
} Clusters_cell;
//---------------------------------------------------------------------------
typedef vector<Clusters_cell> P_cell;//struct that holds Clusters_cell with different timestamp 
//Check whether an object's location is contained within the given rectangle.
struct myclass {
  bool operator() (Cluster c1 ,Cluster c2) { return (c1.get_cluster_population()>c2.get_cluster_population());}
} clusters_compare;
bool pointInRect(double x, double y, box rect)
{
	return ( ((x >= rect.x_min) && (x <= rect.x_max) && (y >= rect.y_min) && (y <= rect.y_max)) ? true : false);
}

//Check if two segments overlap along the same axis, i.e. they have a common interval.
bool segmentOverlap(double a_min, double a_max, double b_min, double b_max)
{
	return ( ((a_max < b_min) || (b_max < a_min)) ? false : true);
}

//Intersection exists if only these rectangles have overlapping extents over both axes.
bool rectIntersect(box rect1, box rect2)
{
	return (segmentOverlap(rect1.x_min, rect1.x_max, rect2.x_min, rect2.x_max) && segmentOverlap(rect1.y_min, rect1.y_max, rect2.y_min, rect2.y_max));
}

//Check whether the first rectangle is fully contained within the second one.
bool rectContain(box rect1, box rect2)
{
	return ((rect1.x_min >= rect2.x_min) && ((rect1.y_min >= rect2.y_min) && (rect1.x_max <= rect2.x_max) && (rect1.y_max <= rect2.y_max)) ? true : false);
}
