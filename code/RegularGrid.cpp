#include "RegularGrid.h"
#include <algorithm>
#include <ctime>
//Constructor #1, getting explicitly the granularity of x,y axes
RegularGrid::RegularGrid(double &pWidth, double &pHeight, unsigned int &pGranX, unsigned int &pGranY) 
{
	//In older implementations, dimensions were artificially expanded by 0.01% ...
	//...to avoid assignments at the boundaries. This practice is NOT used here.
	width = pWidth; 
	height = pHeight;
	GranX = pGranX; 
	GranY = pGranY;
	//Compute cell dimensions
	dx = width/(double)GranX; 
	dy = height/(double)GranY;
}

//Constructor #2, implements square cells, by getting their size
RegularGrid::RegularGrid(double &pWidth, double &pHeight, double &pSize) 
{
	width = pWidth;
	height = pHeight;
	dx = dy = pSize;
	//Compute granularity for axes
	GranX = (unsigned int)ceil(width/dx); 
	GranY = (unsigned int)ceil(height/dy);
}

//Constructor #3, assuming a fixed 2-dimensional Universe, but allowing user to specify the granularity of x,y axes
RegularGrid::RegularGrid(unsigned int &pGranX, unsigned int &pGranY) 
{
        XMIN = X_MIN; YMIN = Y_MIN; XMAX = X_MAX; YMAX = Y_MAX;     //Universe defined in structures.h
	width = XMAX - XMIN;    
	height = YMAX - YMIN;
	GranX = pGranX; 
	GranY = pGranY;
	//Compute cell dimensions
	dx = width/(double)GranX; 
	dy = height/(double)GranY;
}

//Constructor #4, allowing user to specify the granularity of x,y axes and a 2-dimensional Universe
RegularGrid::RegularGrid(unsigned int &pGranX, unsigned int &pGranY, box pBox) 
{
        XMIN = pBox.x_min; YMIN = pBox.y_min; XMAX = pBox.x_max; YMAX = pBox.y_max;   //Universe given as argument
	width = XMAX - XMIN;    
	height = YMAX - YMIN;
	GranX = pGranX; 
	GranY = pGranY;
	//Compute cell dimensions
	dx = width/(double)GranX; 
	dy = height/(double)GranY;
}


//Destructor
RegularGrid::~RegularGrid()
{
	//No need to maintain information about object allocations anymore
	ObjAssignments.clear();

	//Destroy contents for all cells
	for(unsigned int cid=0; cid<cell_cnt; cid++)
		cell[cid].ObjInfo.clear();

	delete[] cell;	 //deallocate Cell table
}

void RegularGrid::printBox(){
	cout<<XMIN<<" "<<XMAX<<" "<<YMIN<<" "<<YMAX<<endl;
}

//Allocate memory and enable the use of this structure.
//It must be called explicitly after the constructor.
bool RegularGrid::Allocate(unsigned int &pobj_cnt) 
{
	obj_cnt = pobj_cnt;
	cell_cnt = GranX*GranY;
	try 
	{
		cell = new GridCell[cell_cnt];		//Allocate Cell table

		//Initialize Cell table
		for (unsigned int i=0; i<cell_cnt; ++i) 
                {
                    cell[i].cellBox = getCellBox(i);     //Store cell coordinates
                    cell[i].processed = false;
                }

		//Initially, assume that all objects are allocated to the first cell
		ObjAssignments.assign (pobj_cnt, 0);
	}
	catch(...) { return false; }

	return true;
}


//Evaluates hash function (resulting into cell ranges between 0 and GranX*GranY-1 ).
//Returns the cell id, i.e. a pointer to the Cell table.
inline unsigned int RegularGrid::HashLocation(double &x, double &y) 
{
    return GranX*(unsigned int)((y-YMIN)/dy)+(unsigned int)((x-XMIN)/dx); 
}
//find center of cell
inline void RegularGrid::printCellCenters()
{
	double x=XMIN;
	double y=YMIN;
	unsigned int cellID=0;
	for(int i=0;i<GranX;i++){
		for(int j=0;j<GranY;j++){
			cout<<cellID<<";"<<x+dx/2<<";"<<y+dy/2<<endl;
			x+=dx;	
			cellID+=1;
		}
		x=XMIN;
		y+=dy;
	}
	
}

/************Cluster in cell with id cellID dummy method**********/
inline void RegularGrid::Clustering_cell(unsigned int &cellID,TObject &tweet,double threshold,unsigned long time)
{
	//P_cell Pc=cell[cellID].ClustInfo;//set of set of clusters
	Clusters_cell Cc,*C;
	bool exists=false;
	//P_cell fresh_batch;
	for(int i=0;i<cell[cellID].ClustInfo.size();i++){
		Clusters_cell Cc=cell[cellID].ClustInfo[i];
		if(Cc.ts==time){
			exists=true;
			break;
		} 
	}
	
	//int N=cell[cellID].ClustInfo.size()-1;
	//if((cell[cellID].ClustInfo[N]).ts==time)exists=true;
	if((cell[cellID].ClustInfo.empty())||(exists==false)){
		//cout<<"cluster new pane"<<endl;
		//P_cell cell[cellID].ClustInfo=cell[cellID].ClustInfo;
		Set new_cluster=tweet.tags;
		Cluster n_C=Cluster(new_cluster);
		n_C.sorting();
		Cc.clusters_set.push_back(n_C);	
		Cc.ts=time;	
		cell[cellID].ClustInfo.push_back(Cc);
		//fresh_batch.push_back(Cc);	
			
				
	}
	else{

	C=&cell[cellID].ClustInfo[cell[cellID].ClustInfo.size()-1];
	double max=-1;
	int max_index=-1;
	int max_i;
	//cout<<1<<endl;
	//cout<<C->clusters_set.size()<<endl;
	for(int i=0;i<(C->clusters_set).size();i++){
		
		double dist=C->clusters_set[i].Jaccard_Index(tweet.tags);//compare clusters with tweet 
		//cout<<dist<<endl;
		if(dist>max){
			max=dist;
			max_i=i;
		}
		
	}
	if(max>=threshold)
		max_index=max_i;
	//cout<<1<<endl;
	if(max_index!=-1){
		
		C->clusters_set[max_index].add_cluster_element(tweet.tags);
		C->clusters_set[max_index].sorting();
		C->ts=time;
		//fresh_batch.push_back(*C);
	}
	else{
		Set new_cluster=tweet.tags;
		Cluster n_C=Cluster(new_cluster);
		//n_C.increase_population(1);
		n_C.sorting();
		C->clusters_set.push_back(n_C);	
		C->ts=time;
		//fresh_batch.push_back(*C);
		}
	
	
    	}
  		
	//cell[cellID].ClustInfo=Pc;		
}

//Updates the new position (x,y) of the specified moving object id (oid)
inline void RegularGrid::UpdateObject(TObject &cur_obj,unsigned long time,double threshold) 
{
	unsigned int i;

	
	
	i = HashLocation(cur_obj.x, cur_obj.y);
		if((i>=0)&&(i<GranX*GranY)){
		//cell[i].processed=true;
		Clustering_cell(i,cur_obj,threshold,time);
		
		
	}	
	
	
		
		
	
}

bool RegularGrid::isNeighboorCell(unsigned int cellA,unsigned int cellB,unsigned int GranX)
{//checks if cell A and cell B are neighboors 
	return ((cellA+1==cellB)||(cellA-1==cellB)||(cellA+GranX==cellB)||(cellA-GranX==cellB)||(cellA+GranX+1==cellB)||(cellA-(GranX+1)==cellB)||(cellA-(GranX-1)==cellB)||(cellA+GranX-1==cellB));
	
}

/** dummy merge method***/
void RegularGrid::merge_panes(multimap<char,Cluster> &merged_clusters,Clusters_cell p2,unsigned int *window_population)//merge 2 clusters sets from different panes
{	
	Clusters_per_cell set_c2=p2.clusters_set;
	Clusters_per_cell set_c;
	if(merged_clusters.empty()){
		//cout<<"merged_cluster is empty"<<endl;
		for(int i=0;i<set_c2.size();i++){
			if(set_c2[i].get_cluster_population()>=0){	
				*window_population+=set_c2[i].get_cluster_population();
				merged_clusters.insert(make_pair(set_c2[i].get_first_char(),set_c2[i]));
			}	
		}	
	}
	else{
		for(int i=0;i<set_c2.size();i++)
		{
			bool is_merged=false;
			Cluster c2=set_c2[i];
			char c = c2.get_first_char();
			//cout<<c<<endl;
			if(c2.get_cluster_population()>=0){
				//*window_population+=c2.get_cluster_population();
				double max_sim=-1;
				multimap<char,Cluster>::iterator max_clust;
				std::pair <std::multimap<char,Cluster>::iterator, std::multimap<char,Cluster>::iterator> ret;
				if(merged_clusters.find(c)!=merged_clusters.end()){
					 ret = merged_clusters.equal_range(c);
				 	for (std::multimap<char,Cluster>::iterator it=ret.first; it!=ret.second; ++it){
			//for(int j=0;j<merged_clusters[c].size();j++){
				/*need to be fix */
			/*	if(c2.Similarity(merged_clusters[j].get_cluster_elements())>=min(c2.get_cluster_elements().size(),merged_clusters[j].get_cluster_elements().size())/2.0){
					merged_clusters[j].merge(c2);
					is_merged=true;
					break;
					//merged_clusters[j]=c;
				}
			
			}*/
			
				double dist=c2.Jaccard_Index((it->second).get_cluster_elements());
				if(dist>max_sim){
					max_sim=dist;
					max_clust=it;
				}
			}
			if(max_sim>=0.5){
				(max_clust->second).merge(c2);
				is_merged=true;
			}
				
			else
					merged_clusters.insert(make_pair(c2.get_first_char(),c2));
				
		}
		else 
			merged_clusters.insert(make_pair(c2.get_first_char(),c2));
	}
	}
		
}
	//cout<<window_population<<endl;	
}
multimap<char,Cluster> RegularGrid::merge_clusters_from_panes(unsigned int &cellID){//merge clusters(similar) from diferenet pane of current window in same cell 
	
	P_cell P;
	P_cell Pc=cell[cellID].ClustInfo;
	multimap<char,Cluster> merged_clusters;
	unsigned int window_population=0;
	for(int i=0;i<Pc.size();i++){
		//int s=clock(); 
		merge_panes(merged_clusters,Pc[i],&window_population);
		//cout<<merged_clusters.size()<<endl;
		//int e=clock();
		//cout<<"simple merge "<<(double)((e-s)*1000)/CLOCKS_PER_SEC<<endl;	
		
	}
	//cout<<window_population<<endl;
	//cout<<merged_clusters.size()<<endl;
	return 	merged_clusters;
				

}
inline int RegularGrid::number_of_clusters(multimap<char,Cluster> m)
{	
	int s=0;
	for(char c='a';c<='z';c++){
		if(m.find(c)!=m.end()){
			s+=m.count(c);
		}
	}
return s;
}			
		
/*void RegularGrid::smarty_merge(vector<Cluster> &merged_clusters,vector<Cluster> fresh_batch,vector<Cluster> old_batch,double threshold){

	vector<int> marked_old(old_batch.size(),0);
	vector<int> marked_fresh(fresh_batch.size(),0);
	for(int i=0;i<merged_clusters.size();i++){
		//cout<<"old_batch "<<old_batch.size()<<endl;
		double max_dist=-1;
		int maxindex=-1;
		for(int j=0;j<old_batch.size();j++){
			if(marked_old[j]==0){
				double dist = merged_clusters[i].Jaccard_Index(old_batch[j].get_cluster_elements());
				if(dist>max_dist){
					max_dist=dist;
					maxindex=j;
				}
			}
		}
		if((max_dist>0)&&((is_subset(merged_clusters[i].get_cluster_elements(),old_batch[maxindex].get_cluster_elements()))||((old_batch[maxindex].get_cluster_elements(),is_subset(old_batch[maxindex].get_cluster_elements(),merged_clusters[i].get_cluster_elements()))))){	
			int N=merged_clusters[i].get_cluster_population()-old_batch[maxindex].get_cluster_population();
			marked_old[maxindex]=1;	
			if(N<=0){

					merged_clusters.erase(merged_clusters.begin()+i,merged_clusters.begin()+i+1);

					i--;

				}

				else{

					Cluster C= Cluster(merged_clusters[i].get_cluster_elements(),N);
					merged_clusters[i]=C;		
				}		
		
			}
		}
	for(int i=0;i<merged_clusters.size();i++){
		double max_dist=-1;
		int maxindex=-1;
		for(int j=0;j<fresh_batch.size();j++){
			if(marked_fresh[j]!=1){
			double dist = merged_clusters[i].Jaccard_Index(fresh_batch[j].get_cluster_elements());
			if(dist>max_dist){
				max_dist=dist;
				maxindex=j;
			}
		}
		}
		//cout<<"new batch"<<endl;
		if((max_dist>threshold)&&((is_subset(merged_clusters[i].get_cluster_elements(),fresh_batch[maxindex].get_cluster_elements()))||is_subset(fresh_batch[maxindex].get_cluster_elements(),merged_clusters[i].get_cluster_elements()))){
		int N=merged_clusters[i].get_cluster_population()+fresh_batch[maxindex].get_cluster_population();
		Cluster C= Cluster(set_union(merged_clusters[i].get_cluster_elements(),fresh_batch[maxindex].get_cluster_elements()),N);
		//fresh_batch.erase(fresh_batch.begin()+maxindex,fresh_batch.begin()+maxindex+1);				
		marked_fresh[maxindex]=1;		
		merged_clusters[i]=C;
		}
	
	}
	for(int i=0;i<fresh_batch.size();i++){
		if(marked_fresh[i]==0)
			merged_clusters.push_back(fresh_batch[i]);
	}
	//cout<<merged_clusters.size()<<endl;
}
 */
/*vector<Cluster> RegularGrid::merge_clusters_from_panes(unsigned int &cellID,int new_time_pane,double threshold){//merge clusters(similar) from diferenet pane of current window in same cell 
	
	//P_cell P;
	//P_cell Pc=cell[cellID].ClustInfo;
	vector<Cluster> fresh_batch;	
	//vector<Cluster> merged_clusters=cell[cellID].past_window;
	//cout<<cell[cellID].old_batch.size()<<endl;
	//cout<<merged_clusters.size()<<endl;	
	if(cell[cellID].old_batch.size()==0){
		cell[cellID].past_window.clear();
		for(int i=0;i<cell[cellID].ClustInfo.size();i++){
			merge_panes(cell[cellID].past_window,cell[cellID].ClustInfo[i],threshold);
		}
		//cout<<"dump way"<<endl;
	}
	else{
		for(int i=0;i<cell[cellID].ClustInfo.size();i++){
		
		if(cell[cellID].ClustInfo[i].ts>new_time_pane){
			//cout<<"new clusters "<<endl;
			for(int j=0;j<cell[cellID].ClustInfo[i].clusters_set.size();j++)
				fresh_batch.push_back(cell[cellID].ClustInfo[i].clusters_set[j]);
			}
		}
		
		
		smarty_merge(cell[cellID].past_window,fresh_batch,cell[cellID].old_batch,threshold);
		//cout<<"smart way"<<endl;
	}
	
	
	return 	cell[cellID].past_window;
				

}
*/
void RegularGrid::printTopicCells(map<char,vector<Topics_Cells> > &A,unsigned int t_front,double max_freq,ofstream &intensity,int &topic_tweets)
{
	//cout<<"*******START WINDOW*********"<<endl;
	//vector<bool> del(A.size(),false);
	for(char c='a';c<='z';c++){
		if(A.find(c)!=A.end()){	
		for(int i=0;i<A[c].size();i++){
		Set topic=A[c][i].name;
		set<pair<unsigned int,unsigned int> > v_cells=A[c][i].cells;
		set<pair<unsigned int,unsigned int> >::iterator it1,it2;
		
		/*int k=0;
		for(it1=v_cells.begin();it1!=v_cells.end();++it1){
			for(it2=v_cells.begin();it2!=v_cells.end();++it2){
				if(isNeighboorCell((*it1).first,(*it2).first,GranX))
					k++;
			}
		}*/
		intensity<<t_front<<";";
		unsigned int sum=0;
		for(it1=v_cells.begin();it1!=v_cells.end();++it1){
				intensity<<"'";
				for(int h=0;h<topic.size();h++){
					intensity<<topic[h]<<" ";
				}
				intensity<<"';";
				intensity<<(*it1).first<<";"<<(*it1).second<<endl;
				sum+=(*it1).second;
		}
		//intensity<< v_cells.size()<<";";
		topic_tweets+=sum;
		//cout<<sum<<endl;	
		if(sum<max_freq){
			A[c].erase(A[c].begin()+i,A[c].begin()+i+1);
			i--;
		}
		
		/*}
		cout<<"{";
		for(int c=0;c<topic.size();c++){
			cout<<topic[c]<<" ";
		}
		cout<<"};{";
		for(it1=v_cells.begin();it1!=v_cells.end();++it1){
			
			cout<<"("<<(*it1).first<<","<<(*it1).second<<")"<<" ";
			
		}
		*/
	}	
	}	
	}	
	
}
void RegularGrid::find_top_topics(double threshold,int window_tweets,double percent,unsigned long window_rear_t,long front_t,map<char,vector<Topics_Cells> > &v_topic_cells,map<char,vector<Topics_Cells>> past_topics,int n_time_pane,bool tolerance_condition)//rear time value of window -->range-front_t_pane
{	
	multimap<char,Cluster> m;
	int s1=0;
	for (unsigned int i=0; i<cell_cnt; ++i){

		//if(cell[i].processed){

			 m=find_topics(i,threshold,window_tweets,percent,window_rear_t,v_topic_cells,past_topics,n_time_pane,tolerance_condition,s1);
			 // cout<<"past "<<cell[i].past_window.size()<<endl;
			//cout<<"
		//	m.clear();
		//sort(m.begin(),m.end(),clusters_compare);

		/*if(i==1272){

			
			//cout<<"top-"<<k<<" topic is: "<<endl;

			//int size=min(m.size(),(unsigned long)k);// if the elements are less k  
			
				for(int j=0;j<m.size();j++){
									
					
						Set S=m[j].get_cluster_elements();
						
						cout<<front_t<<";"<<i<<";{";

						for(int c=0;c<S.size();c++){
		
							cout<<S[c];

							if(c<S.size()-1)

								cout<<",";

						}
	
						cout<<"}"<<";"<<m[j].get_cluster_population()<<endl;

					
				}
		}*/
	   //}
	/****print topics in cell cellID***/
	}
	cout<<s1<<endl;
	
}
/*inline bool RegularGrid::popular_topic_condition(Cluster c,bool tolerance_condition,double max_frequency,map<char,vector<Cluster> > past_topics)
{
	bool tol;
	//double per_km = (2500.0)/(50*50);
	if(past_topics.size()>0){
	if(tolerance_condition){
		tol= (((c.get_cluster_population())>=max_frequency)||(tolerance(c,past_topics,0.2,max_frequency)));
		
	}
	else{
		tol= ((c.get_cluster_population())>=max_frequency);
	}
	}
	else{
		tol= ((c.get_cluster_population())>=max_frequency);
	}
	return tol;
}*/
inline multimap<char,Cluster> RegularGrid::popular_topics(multimap<char,Cluster> topics,bool tolerance_condition,int window_tweets,double percent,map<char,vector<Topics_Cells> > past_topics)
{
	multimap<char,Cluster> clean_topics;
	/*for(int i=0;i<topics.size();i++){
		if(popular_topic_condition(topics[i],tolerance_condition,max_frequency,past_topics)){
			clean_topics.push_back(topics[i]);
		}
	}
	*/
	for(char c='a';c<='z';c++){
		if(topics.find(c)!=topics.end()){
		  std::pair <std::multimap<char,Cluster>::iterator, std::multimap<char,Cluster>::iterator> ret;
    		  ret = topics.equal_range(c);
		  for (std::multimap<char,Cluster>::iterator it=ret.first; it!=ret.second; ++it){
			Cluster Cl = it->second;
			//cout<<(float)Cl.get_cluster_population()/window_tweets<<endl;
			if(Cl.get_cluster_population()>=(percent/100.0)*window_tweets){
				clean_topics.insert(make_pair(c,Cl));	
	}
	}
	}	
	}
	return clean_topics;
}
		
			
multimap<char,Cluster> RegularGrid::find_topics(unsigned int &cellID,double threshold,int window_tweets,double percent,unsigned long time,map<char,vector<Topics_Cells> > &v_topic_cells,map<char,vector<Topics_Cells> > past_topics,int n_time_pane,bool tolerance_condition, int &s)
{	
	
	//remove_oldest_pane(time);//remove panes that doesn't belong to current window
	remove_past_clusters(cellID,time);
	//unsigned int pop_topics_s=get_time();
	//int s=clock();
	//vector<Cluster> m=merge_clusters_from_panes(cellID,n_time_pane,threshold);
	//int e=clock();
	multimap<char,Cluster> m=merge_clusters_from_panes(cellID);
	//unsigned int pop_topics_e=get_time();
	//cout<<pop_topics_e-pop_topics_s<<endl;
	
	s+=number_of_clusters(m);
	multimap<char,Cluster> clean_topics=popular_topics(m,tolerance_condition,window_tweets,percent,past_topics);
	//for(int i=0;i<m.size();i++){
		//if(popular_topic_condition(m[i],tolerance_condition,max_frequency,past_topics)){
			//clean_topics.push_back(m[i]);
		//unsigned int find_topics_s=get_time();
		//unsigned int pop_topics_s=get_time();
		//cout<<clean_topics.size()<<endl;
	//int s=clock();
	//cout<<clean_topics.size()<<endl;
		for(char c='a';c<='z';c++){
			 std::pair <std::multimap<char,Cluster>::iterator, std::multimap<char,Cluster>::iterator> ret;
			if(clean_topics.find(c)!=clean_topics.end()){
    		  	ret = clean_topics.equal_range(c);
		  for (std::multimap<char,Cluster>::iterator it=ret.first; it!=ret.second; ++it){
			//cout<<it->first<<endl;
			bool topic_exists=false;
			if(v_topic_cells.size()>0){
				int maxindex=-1;
				double max_sim=-1;
				for(int j=0;j<v_topic_cells[(it->first)].size();j++){
					double dist = ((it->second)).Jaccard_Index(v_topic_cells[(it->first)][j].name);
					if(dist>max_sim){
						max_sim=dist;
						maxindex=j;
					}
				}
				if(max_sim>=threshold){
					topic_exists=true;

					Set a=set_union(v_topic_cells[(it->first)][maxindex].name,((it->second)).get_cluster_elements());
					pair<unsigned int,unsigned int> tuple(cellID,((it->second)).get_cluster_population());
					v_topic_cells[(it->first)][maxindex].name=a;
					v_topic_cells[(it->first)][maxindex].cells.insert(tuple);
					
				}
			}
			if(topic_exists==false){
				pair<unsigned int,unsigned int> tuple(cellID,((it->second)).get_cluster_population());
				set<pair<unsigned int,unsigned int> > v_cellID;
				v_cellID.insert(tuple);
				Topics_Cells t;
				t.name=((it->second)).get_cluster_elements();
				t.cells.insert(tuple);
				v_topic_cells[(it->first)].push_back(t);					
			}
		//}
	}
	}
	}
	//cout<<(double)(clock()-s)/CLOCKS_PER_SEC<<endl;
	return clean_topics;
}	

/*inline bool RegularGrid::tolerance(Cluster topic,map<char,vector<Topics_Cells> > previous_topics,double sim_tolerance,double max_frequency){
	int maxindex=-1;
	double max_sim=-1;
	if(topic.get_cluster_population()>=(1.0-sim_tolerance)*max_frequency){
		
		for(int i=0;i<previous_topics.size();i++){
			double dist = topic.Jaccard_Index(previous_topics[i].name);
			if(dist>0.5){
				max_sim=dist;
				maxindex=i;
				break;
			}
		}
	if(max_sim>0.5)
		return true;
	else return false;
	}
	else 
		return false;

}*/
void RegularGrid::remove_oldest_pane(unsigned long t)
{
	for (unsigned int i=0; i<cell_cnt; i++){ 
		//if(cell[i].processed)
			remove_past_clusters(i,t);//remove clusters_cell  in cell i that have timestamp value less than t		
	}
}
void RegularGrid::remove_past_clusters(unsigned int &cellID,unsigned long time)
{	//remove clusters that cluster.ts<time
	
	vector<Cluster>past_batch;
	for(int i=0;i<cell[cellID].ClustInfo.size();i++){
		//cout<<"Pc["<<i<<"]"<<".ts= "<<Pc[i].ts<<"time="<<time<<endl;
		if(cell[cellID].ClustInfo[i].ts<time+1){
			//int j=i;
			for(int j=0;j<cell[cellID].ClustInfo[i].clusters_set.size();j++){	
				past_batch.push_back(cell[cellID].ClustInfo[i].clusters_set[j]);
			}
			cell[cellID].ClustInfo.erase(cell[cellID].ClustInfo.begin()+i,cell[cellID].ClustInfo.begin()+i+1);//erase oldest pane 		
			i--;
	
			}
	}
	cell[cellID].old_batch.swap(past_batch);
	past_batch.clear();
			
}

//Returns a list of cells overlapping with the given MBB
inline map<unsigned int, bool> RegularGrid::HashBox(box &pBox) 
{
    map<unsigned int, bool> cellList;               //List of cells returned with indicators for partial overlap
   
    unsigned int lc, uc, c, i, gdx;
    bool partialCover = false;
    box cellBox;

    //First, hashing box corners
    lc = HashLocation(pBox.x_min, pBox.y_min); 
    uc = HashLocation(pBox.x_max, pBox.y_max);
        
    //Range of cells affected in each row of the grid
    gdx = uc%GranX - lc%GranX;
    
    //Find all cells covered by this box after examining the cells of its corners
    for(c=lc; c<=uc-gdx; c+=GranX) 	
        for(i=c; i<=c+gdx; ++i) 	
        {	
            if (i >= cell_cnt)		//i is the cell currently being examined for overlap
                continue;
            
            //Check if cell box is fully contained within the query range
            partialCover = !rectContain(cell[i].cellBox, pBox);
                	
            //Insert this cell into the chain according to its total/partial overlap	
            cellList.insert(pair<unsigned int, bool>(i, partialCover));
        }

    return cellList; 
}


//Calculate cell box coordinates (may be also used for printing them into KML)
inline box RegularGrid::getCellBox(unsigned int &cellID)
{
    //Cell matrix indices along x and y axes
    unsigned int i = cellID % GranX;   
    unsigned int j = cellID / GranX;
    
    double x_min = XMIN + i * dx;
    double y_min = YMIN + j * dy;
    double x_max = XMIN + (i+1) * dx;
    double y_max = YMIN + (j+1) * dy;

    return box(x_min, y_min, x_max, y_max);
}


//Type all contents of each cell
void RegularGrid::printGridState()
{	//unsigned int i=42;
	for (unsigned int i=0; i<cell_cnt; ++i) 
		getCellObjects(i);				//Currently retained object locations
}


//Type object locations currently contained in the given cell
void RegularGrid::getCellObjects(unsigned int &cellID)
{
	//Object contents
	map < unsigned int, TObject * > objSequence = cell[cellID].ObjInfo;
	map < unsigned int, TObject * >::iterator objIter;
	P_cell clus = cell[cellID].ClustInfo;
	//Iterate through the sequence of objects in that cell
	/*for( objIter = objSequence.begin(); objIter != objSequence.end(); objIter++ ) 
	{
		TObject *o = objIter->second;
		cout << "CELL " << cellID << " OBJ <" << o->ts << " " << o->id << " " << o->x << " " << o->y << ">" << endl;
	}
	*/
	//cout<<clus.size()<<endl;
	for(int i=0;i<clus.size();i++){
		
	//	cout<<clus.size()<<endl;
		Clusters_cell C=clus[i];
		Clusters_per_cell A=C.clusters_set;
		for(int j=0;j<A.size();j++){
			//cout<<1<<endl;
			//cout<<"Cluster in cell "<<cellID<<" and timestamp "<<C.ts<<" has elements {";
				Set Cl=A[j].get_cluster_elements();
				int freq=A[j].get_cluster_population();
				for(int k=0;k<Cl.size();k++){
					if(k<=Cl.size()-1)					
						cout<<Cl[k]<<" , ";
				}
				
				//cout<<"} and frequency "<<freq<<endl;
		}
		
	}
		
	
}

bool RegularGrid::belongs_to_cluster(unsigned int cellID,set<unsigned int> C)
{	
	bool belong=false;
	for(set<unsigned int>::iterator it=C.begin();it!=C.end();++it){
		if(isNeighboorCell(cellID,*it,GranX)){
			belong=true;
			break;
		}
	}
	return belong;
}
set<unsigned int> RegularGrid::update(unsigned int cellID,set<unsigned int> cells)
{
		set<unsigned int> cluster;
		cells.insert(cellID);
		cluster=cells;
		return cluster;
}
void RegularGrid::merge_changed_clusters(set<unsigned int> &M,set<unsigned int> cells)
{
	for(set<unsigned int>::iterator t=cells.begin();t!=cells.end();++t){
		M.insert(*t);
	}
}
vector<topic_zones> RegularGrid::findDenseCells(map<char,vector<Topics_Cells> > T)
{
	vector<topic_zones> tz;
	topic_zones res;
	//int start_t=clock();
	for(char c='a';c<='z';c++){
		if(T.find(c)!=T.end()){
		for(int i=0;i<T[c].size();i++){
		//cout<<1<<endl;
			res.topic_name=T[c][i].name;//set of hashtags --> topic name 
			set<pair<unsigned int,unsigned int> > cells=T[c][i].cells;//set of cells id in topic T[i]
			for(set<pair<unsigned int,unsigned int> >::iterator it=cells.begin();it!=cells.end();++it){
			//cout<<1<<endl;
				if(res.cluster_cells.size()>=1){
					bool merged=false;
					set<unsigned int>  changed_clusters;//initialize here the set of changed clusters to empty
					set<unsigned int> Cl;				
					for(int j=0;j<res.cluster_cells.size();j++){
					//cout<<1<<endl;
					//set<unsigned int>  changed_clusters;
						if(belongs_to_cluster((*it).first,res.cluster_cells[j]))//check if cell id belongs to certain cluster
						{
							merged=true;
							//cout<<1<<endl;
							Cl=update((*it).first,res.cluster_cells[j]);//add cell id to cluster
							res.cluster_cells.erase(res.cluster_cells.begin()+j,res.cluster_cells.begin()+j+1);
							//merge clusters where *it cell id updated
							j--;
							merge_changed_clusters(changed_clusters,Cl);//update set of changed clusters  					
						}	
					}
						Cl.clear();
					//changed_clusters.clear();	
					//update set of topic clusters
					//res.cluster_cells.erase(res.cluster_cells.begin()+j,res.cluster_cells.begin()+j+1);
						
				
					if(merged==false){
						set<unsigned int> cluster_c;
							cluster_c.insert((*it).first);
							bool exists=false;
							for(int k=0;k<res.cluster_cells.size();k++){
								if(res.cluster_cells[k]==cluster_c){
									exists=true;
									break;
							}
						}				
							if(!exists) 
								res.cluster_cells.push_back(cluster_c);
						
				}
				else {
						res.cluster_cells.push_back(changed_clusters);
				}		
			}
			else{
						
						set<unsigned int> cluster_c;
						cluster_c.insert((*it).first);
						res.cluster_cells.push_back(cluster_c);
								
			}
		}
		res.freq=0;
		for(set<pair<unsigned int,unsigned int> >::iterator it=cells.begin();it!=cells.end();++it){
				res.freq+=(*it).second;		
		}
		tz.push_back(res);
		res.cluster_cells.clear();
	}
	}
	}
	/*int stop_t=clock();
	cout<<"time elapsed : "<<stop_t-start_t<<" ms"<<endl;
	*/	
	return tz;
}

		
void RegularGrid::printTopicZones(vector<topic_zones> T,unsigned int t_front,ofstream &results)
{
	for(int i=0;i<T.size();i++){
		Set name=T[i].topic_name;
		//print name of topic
		results<<t_front<<";";
		results<<"{"; 
		for(int n=0;n<name.size();n++)
			results<<name[n]<<" ";	
		results<<"};{";
		vector<set<unsigned int> > clusters = T[i].cluster_cells;
		int s=0;
		for(int j=0;j<clusters.size();j++){
			results<<"{";
			for(set<unsigned int>::iterator it=clusters[j].begin();it!=clusters[j].end();++it){
				results<<*it<<" ";
				s++;
			}
			results<<"} ";
		}
		results<<"};";
		results<<T[i].freq<<";";
		results<<(float)s/clusters.size()<<endl;
	}
}
//Perform a range query using a rectangle against the grid index
void RegularGrid::rangeQuery(box &qBox, vector<long> &qryResults) 
{
    map < unsigned int, TObject * > objSequence;
    map < unsigned int, TObject * >::iterator objIter;
    
    //Hash query rectangle in order to find the proper cell(s)
    map<unsigned int, bool> cellList = HashBox(qBox);
    
    //Identify object locations involved in the cell(s) just found
    map<unsigned int, bool>::iterator cellIter;
    for (cellIter=cellList.begin(); cellIter != cellList.end(); cellIter++ )
        if (cellIter->second)       //Query range partially overlapping this cell
            refineCandidates(cellIter->first, qBox, qryResults);           //Refinement is necessary
        else                        //Query range totally covering this cell
        {   //So, report all objects within this cell
            objSequence = cell[cellIter->first].ObjInfo;
            for( objIter = objSequence.begin(); objIter != objSequence.end(); objIter++ )
                qryResults.push_back(objIter->second->id);
        }
}

//Check each object indexed by the given cell for possible containment within the given query rectangle
void RegularGrid::refineCandidates(const unsigned int &cellID, box &qBox, vector<long> &qryResults)
{
    map < unsigned int, TObject * >::iterator objIter;
    
    //Object locations indexed in this cell
    map < unsigned int, TObject * > objSequence = cell[cellID].ObjInfo;
    
    //Iterate through this sequence (TOTALLY contained)
    for( objIter = objSequence.begin(); objIter != objSequence.end(); objIter++ ) 
    {
        if (pointInRect(objIter->second->x, objIter->second->y, qBox))      //Location within with query rectangle
            qryResults.push_back(objIter->second->id);                     //... so report qualifying object
    }
}
