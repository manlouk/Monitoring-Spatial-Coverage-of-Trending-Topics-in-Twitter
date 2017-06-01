#include "post_processing.h"

using namespace std;
pair< set<unsigned int>,set<unsigned int> > printCellDiffs(set< pair<unsigned int,unsigned int> > cells_C,set< pair<unsigned int,unsigned int> > cells_P)
{
	pair<set<unsigned int>,set<unsigned int> > g;
	set<unsigned int> left;
	set<unsigned int> arrived;
	set<pair<unsigned int,unsigned int> >::iterator c,p;
	for(p=cells_P.begin();p!=cells_P.end();++p){
		bool exists=false;
		for(c=cells_C.begin();c!=cells_C.end();++c){
			if((*c).first==(*p).first){
				exists=true;
				break;
			}
		}
		if(!exists){
			left.insert((*p).first);
			cells_P.erase(*p);
			//p--;
		}
		else{
			cells_P.erase(*p);
			cells_C.erase(*c);
			//p--;
			//c--;		
		}
	}
	for(c=cells_C.begin();c!=cells_C.end();++c){
		arrived.insert((*c).first);
	}
	g.first=arrived;
	g.second=left;
	return g;
}

void TopicMonitoring(map<char,vector<Topics_Cells> >Curr,map<char,vector<Topics_Cells> > past,double theta,unsigned int t_front,ofstream &monitor)
{	
	//cout<<Curr.size()<<endl;
	
	for(char c='a';c<='z';c++){
		if(Curr.find(c)!=Curr.end()){
		for(int i=0;i<Curr[c].size();i++){
			set<unsigned int> arr;
			set<unsigned int> left;
			pair<set<unsigned int>,set<unsigned int> > p;		
			Set topic = Curr[c][i].name;
			int past_index=-1;	
			if(past.find(c)!=past.end()){
				for(int j=0;j<past[c].size();j++){
					if(similarity(topic,past[c][j].name)>=theta){
						past_index=j;
						break;
					}	
				}
			}
			
			if(past_index>-1){
				//cout<<1<<endl;	
				p = printCellDiffs(Curr[c][i].cells,past[c][past_index].cells);

			}
	
			else{
				for(set<pair<unsigned int,unsigned int> >::iterator t=(Curr[c][i].cells).begin();t!=(Curr[c][i].cells).end();++t){
					arr.insert((*t).first);	
			}
				 p.first=arr;
				 p.second=left;
			
			}
		//print differences
		//monitor<<t_front<<";Topic:{ ";
		if((p.first.size()!=0)||(p.second.size()!=0)){
		monitor<<t_front<<";";
		monitor<<"{";
		for(Set::iterator n=topic.begin();n!=topic.end();++n){
			monitor<<*n<<" ";
		}
		monitor<<"};";
		monitor<<"{";
		for(set<unsigned int>::iterator t1=p.first.begin();t1!=p.first.end();++t1){
			monitor<<*t1<<" ";
		}
		monitor<<"};";
		monitor<<"{";
		for(set<unsigned int>::iterator t2=p.second.begin();t2!=p.second.end();++t2){
			monitor<<*t2<<" ";
		}
		monitor<<"}"<<endl;
		
		//p.clear();	
		arr.clear();
		left.clear();
	}
	}
	}
	}
}
unsigned int number_of_topics(map<char,vector<Topics_Cells> > A){
	int sum=0;	
	for(char c='a';c<='z';c++){
		sum+=A[c].size();
	}
	return sum;
}
