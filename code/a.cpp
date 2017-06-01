#include "post_processing.cpp"
#include <chrono>
using namespace std;
using namespace chrono;
ifstream in;
ofstream out;
ofstream results;
ofstream monitor;
ofstream intensity;
bool first=true;
unsigned long long timestamp,id;
double x,y;
unsigned int gx;
unsigned   int obj_cnt=200000000;

int slide;
int range;
string tags;
TObject tweet;
//borders of grid attempt 1 
/*int minX = 480000;
int maxX = 580000;
int minY = 130000;
int maxY = 230000;
*/
//attempt 2
unsigned  int minX = 505000;
unsigned  int maxX = 555000;
unsigned  int minY = 155000;
unsigned  int maxY = 205000;

map<char,vector<Topics_Cells> > v_topic_cells;
map<char,vector<Topics_Cells> > past_topics;
double theta; //Jaccard threshold
double phi;//tweets/sq.km
double max_frequency;//threshold to define popular topic
bool tolerance_condition;//existance of tolerance condition for finding popular topics  
string output_file;

int main(int argc,char *argv[])
{
	//cout<<obj_cnt<<endl;	
	//in.open("/home/manolis/Dropbox/diploma/diploma_code/DataSets/london_tweets_osgb36.csv");
	//in.open("/home/manolis/Dropbox/diploma/diploma_code/DataSets/tweetsx2.csv");
	//in.open("/home/manolis/Dropbox/diploma/diploma_code/DataSets/tweetsx5.csv");
	//in.open("/home/manolis/Dropbox/diploma/diploma_code/DataSets/tweetsx10.csv");
      	 // in.open("/home/manolis/Dropbox/diploma/diploma_code/DataSets/tweets20.csv");
	//in.open("/home/manolis/Dropbox/diploma/diploma_code/DataSets/tweets50_a.csv");
	//in.open("/home/manolis/Dropbox/diploma/diploma_code/DataSets/tweets100.csv");
	if(argc<9){
		cout<<"Not enough input parameters"<<endl;
		cout<<"The correct parameters for the input is the following list:"<<endl;
		cout<<"1.Grid dimensions NxN(integer)"<<endl;
		cout<<"2.The range of the sliding window(seconds)"<<endl;
		cout<<"3.The slide of the time frame(pane)(seconds and multiple value of range)"<<endl;
		cout<<"4.Threshold for Jaccard distance measure in Clustering procedure(double)"<<endl;
		cout<<"5.Threshold to filter the popular topics(interger)"<<endl;
		cout<<"6.Choose if you want to check popular topics from previous execution cycle to find popular topics for current cycle(y/n)"<<endl; 
		cout<<"7.Initiation timestamp for windows(integer)"<<endl;
		exit(-1);
	}  	
	string file = argv[1];
	in.open(file);
	gx=atoi(argv[2]);
	range=atoi(argv[3]);
	slide=atoi(argv[4]);
	theta=atof(argv[5]);
	phi=atof(argv[6]);
	
	
	if(argv[7]==string("y")){
		 output_file = to_string(gx)+"_"+to_string(range)+"_"+to_string(slide)+"_"+to_string(theta)+"_"+to_string(phi)+"_y";
		tolerance_condition=true;
	}
	else{
		 output_file = to_string(gx)+"_"+to_string(range)+"_"+to_string(slide)+"_"+to_string(theta)+"_"+to_string(phi)+"_n";
		tolerance_condition=false;
	}
	unsigned long  t_curr=atoi(argv[8]);//initiation time for window 1388361600
	int i = file.find("tweets");
	
	string infix = file.substr(i+6,3);
	if(infix.find(".")!=string::npos)
		infix = file.substr(i+6,2);
	//cout<<infix<<endl;
	//exit(0);
	//out.open("measures_"+infix+"_"+output_file);
	results.open("results_"+infix+"_"+output_file);
	//monitor.open("monitor_"+infix+"_"+output_file);
	//intensity.open("intensity_"+infix+"_"+output_file);
	
	box B=boxInfo(minX,minY,maxX,maxY);
	/**grid initilization**/	
	RegularGrid *grid = new RegularGrid(gx,gx,B);   	
	if (!grid->Allocate(obj_cnt)) 
	{
		cout << "Memory Allocation failed for grid partitioning!" << endl; 
		return 1;
	}
	/** regular **/
	//unsigned long t_curr=1388361600;
	unsigned long  start_window=t_curr;
	unsigned long  t_front = t_curr+slide;
	double total_find_top=0;
	double total_zones=0;
	int k1=0;
	int k2=0;
	int k3=0;
	
	/***read tweet ****/ 
	in >> timestamp >> id >> x >> y >> tags;
	tweet.ts=timestamp;
	tweet.id=id;
	tweet.x=x;
	tweet.y=y;
	tweet.tags=gethashtags(tags);
	int number_of_panes=0;
	
	//max_frequency=phi*(((double)(maxX-minX)*(maxY-minY))/(gx*gx*1000*1000));//normalize area expressed in sq.kilometers
	//cout<<max_frequency<<" "<<(maxX-minX)*(maxY-minY)<<" "<<gx<<" "<<((double)(maxX-minX)*(maxY-minY))/((double)(gx*gx*1000*1000))<<endl;
	
	//out<<"Window Range;Topic Clustering;Threshold filtering;Dense areas detection;Post Processing;Number of Topics"<<endl;
	//monitor<<"Window Range;Topic;included;excluded"<<endl;
	intensity<<"Window Range;Topic;Cell;intensity"<<endl;
	int N = range/slide;
	//cout<<N<<endl;
	int batch[N];
	for(int i=0;i<N;i++)batch[i]=0;
	//int s1=0;
	while(true){
		unsigned int tweets=0;
		double total_time_clust=0;
		high_resolution_clock::time_point t1 = high_resolution_clock::now();
		while(tweet.ts<=t_front){
			
			grid->UpdateObject(tweet,t_front,theta);//put tweet on grid 
			//unsigned int stop_clustering=get_time();
			tweets++;
			
			/****read another tweet in the same pane ****/
			in >> timestamp >> id >> x >> y >> tags;
			if(in.eof())break;
			tweet.ts=timestamp;
			tweet.id=id;
			tweet.x=x;
			tweet.y=y;
			tweet.tags=gethashtags(tags);
		}
		//cout<<t_front<<" "<<tweets<<endl;
		//cout<<1<<endl;
		batch[number_of_panes]=tweets;
		number_of_panes++;
		int window_tweets=0;
		for(int i=0;i<N;i++)window_tweets+=batch[i];
		cout<<t_front<<";"<<window_tweets<<";";
		//cout<<t_front<<" "<<(phi/100.0)*window_tweets<<endl;
		high_resolution_clock::time_point t2 = high_resolution_clock::now();
		 auto duration = duration_cast<microseconds>( t2 - t1 ).count();
		total_time_clust+=duration/1000.0;
		/***the pane ended***/	
		//out<<"Window Range;Tweets Clustering;Topic Clustering;Dense Areas;Post Processing;Number of Topics";
		out<<t_front<<";";
		//clustering time
		out<<total_time_clust<<";";
		//unsigned int start_topic=get_time();
		t1 = high_resolution_clock::now();
		grid->find_top_topics(theta,window_tweets,phi,t_front-range,t_front,v_topic_cells,past_topics,t_front-slide,tolerance_condition);
		 t2 = high_resolution_clock::now();
		//unsigned int end_topic=get_time();
		  duration = duration_cast<microseconds>( t2 - t1 ).count();
		total_find_top=duration/1000.0;
		//print find topics time
		out<<total_find_top<<";";
		//unsigned int start_zones=get_time();
		 t1 = high_resolution_clock::now();		
		vector<topic_zones> Zones_per_topic=grid->findDenseCells(v_topic_cells);
		 t2 = high_resolution_clock::now();
		  duration = duration_cast<microseconds>( t2 - t1 ).count();
		
		//unsigned int end_zones=get_time();
		//total_zones=end_zones-start_zones;
		//print Topic Zones time 
		//out<<total_zones<<";";
		out<<duration/1000.0<<";";		
		//unsigned int mon_s=get_time();
		 t1 = high_resolution_clock::now();
			
		TopicMonitoring(v_topic_cells,past_topics,theta,t_front,monitor);
		 t2 = high_resolution_clock::now();
		  duration = duration_cast<microseconds>( t2 - t1 ).count();
		out<<duration/1000.0<<";";			
		//unsigned int mon_e=get_time();
	 	grid->printTopicZones(Zones_per_topic,t_front,results);
		//unsigned int total_monitor=mon_e-mon_s;
		//print Topics monitoring time(post processing)
		//out<<total_monitor<<";";
		//out<<tweets<<";";
		out<<number_of_topics(v_topic_cells)<<";";
		//out<<past_topics.size()<<endl;
		Zones_per_topic.clear();
		past_topics.clear();
		int topic_tweets=0;
		grid->printTopicCells(v_topic_cells,t_front,max_frequency,intensity,topic_tweets);
		out<<topic_tweets<<endl;
		past_topics.swap(v_topic_cells);	
		v_topic_cells.clear();
		if(number_of_panes==(range/slide)){
			
			number_of_panes=0;//initialize for next window
			start_window=t_front-range+slide;
			
		}
		
		t_front+=slide;
		if(in.eof())break;
	}
	
	cout<<"Geosocial Monitoring finished"<<endl;
	in.close();
	out.close();
	results.close();
	monitor.close();
	return 0;
}			
			
			
	

	
	
