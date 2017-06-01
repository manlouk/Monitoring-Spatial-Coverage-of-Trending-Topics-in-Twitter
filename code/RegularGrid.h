//Title: RegularGrid.h
//Description: Regular grid partitioning for manipulating point data streams.
//Platform(s): gcc 4.8.3.
//Date: 29/12/2007
//Revision: 17/3/2015

#include "structures.h"

class RegularGrid {
public:
  RegularGrid(double &, double &, unsigned int &, unsigned int &);	//1st Constructor
  RegularGrid(double &, double &, double &);				//2nd Constructor
  RegularGrid(unsigned int &, unsigned int &);                          //3rd Constructor
  RegularGrid(unsigned int &, unsigned int &, box);                     //4th Constructor
  ~RegularGrid();							//Destructor	

  bool Allocate(unsigned int &);
  inline void UpdateObject(TObject &,unsigned long,double);
  void printGridState();
  void remove_past_clusters(unsigned int&,unsigned long);
  void remove_oldest_pane(unsigned long);
  void getCellObjects(unsigned int &);
  void rangeQuery(box &, vector<long> &);
  unsigned int ts;				//Current timestamp value (Heartbeat)
  unsigned int cell_cnt;
  void find_top_topics(double,int,double,unsigned long,long,map<char,vector<Topics_Cells> >&,map<char,vector<Topics_Cells> >,int ,bool);
 inline void printCellCenters();
void printTopicCells(map<char,vector<Topics_Cells> >&,unsigned int,double,ofstream &, int &);
void printBox();
vector<topic_zones> findDenseCells(map<char,vector<Topics_Cells> >);
void printTopicZones(vector<topic_zones>,unsigned int,ofstream & );
//void TopicMonitoring(vector<Topics_Cells>,vector<Topics_Cells>);	
//void topic_expansion(vector<Topics_Cells>,vector<Topics_Cells>);	
//Total number of cells

private:
  //Cell contents. Note that pointer references to object locations are used in cell listings (data type TObjChain)
  struct GridCell 
  {
	TObjChain ObjInfo;    //Object locations assigned into this cell
	P_cell ClustInfo;        
	box cellBox;          //Cell rectangle specified by its coordinates
	bool processed;
	vector<Cluster> old_batch;
	vector<Cluster> new_batch;
	vector<Cluster> past_window;
  };
  inline int number_of_clusters(multimap<char,Cluster>);
  vector<unsigned int> ObjAssignments;	//For each object, remember the cell it has been allocated to.
   inline void Clustering_cell(unsigned int &,TObject &,double,unsigned long);	
  inline unsigned int HashLocation(double &, double &);
  inline map<unsigned int, bool> HashBox(box &);
  inline box getCellBox(unsigned int &);
  void refineCandidates(const unsigned int &, box &, vector<long> &);
  multimap<char,Cluster>  find_topics(unsigned int &,double,int,double,unsigned long,map<char,vector<Topics_Cells> >&,map<char,vector<Topics_Cells> >,int,bool,int &);
  unsigned int obj_cnt;             //Number of distinct objects indexed in the grid
 //vector<Cluster> merge_clusters_from_panes(unsigned int &,int,double);
multimap<char,Cluster> merge_clusters_from_panes(unsigned int &);
 void merge_panes(multimap<char,Cluster>&,Clusters_cell,unsigned int *);
  double XMIN, YMIN, XMAX, YMAX;    //Space bounds (universe of discourse)  
  unsigned int GranX, GranY;        //granularity of each dimension for hashing
  double width, height;             //width, height of the 2D space.  
  double dx, dy;                    //the dimensions of each cell
  GridCell *cell;                   //Cell table
  bool isNeighboorCell(unsigned int ,unsigned int ,unsigned int );
 void merge_changed_clusters(set<unsigned int> &,set<unsigned int> );
//void smarty_merge(vector<Cluster>& ,vector<Cluster> ,vector<Cluster>,double );
 set<unsigned int> update(unsigned int cellID,set<unsigned int>);
bool belongs_to_cluster(unsigned int ,set<unsigned int> );
//inline bool tolerance(Cluster,vector<Topics_Cells>,double,double);
//inline bool popular_topic_condition(Cluster ,bool,double,vector<Topics_Cells>);
inline multimap<char,Cluster> popular_topics(multimap<char,Cluster>,bool,int,double ,map<char,vector<Topics_Cells> >);
};
