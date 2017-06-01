#include<algorithm>
#include<vector>
#include<string>

using namespace std;
typedef vector<string> Set;

bool is_member(Set A,string s){

	int i=0;
	bool flag=false;
	while(i<A.size()){
		if(s.compare(A[i])==0){
			flag=true;
			break;
		}
		i++;
	}
	return flag;	
			
}

Set set_union(Set A,Set B){
	int a=A.size();
	int b=B.size();
	Set C; 
	if(a<b){
		for(int i=0;i<a;i++){
			if(is_member(B,A[i]))continue;
			B.push_back(A[i]);
		}
	C=B;
	
	}
	else{
		for(int i=0;i<b;i++){
			if(is_member(A,B[i]))continue;
			A.push_back(B[i]);
		}
	 C=A;
	
	}
	return C;
}

Set set_intersection(Set A,Set B){
	int a=A.size();
	int b=B.size();
	Set C; 
	string s;
	if(a<b){
		for(int i=0;i<a;i++){
			 s=A[i];
			if(is_member(B,s))
				C.push_back(s);
		}
	}
	else{
		for(int i=0;i<b;i++){
			string s=B[i];
			if(is_member(A,s))
				C.push_back(s);
		}
	}
	return C;
}
bool is_subset(Set A,Set B )
{	
	bool sub=false;
	if(A.size()<B.size()){
		sub=true;
		for(int i=0;i<A.size();i++){
			if(!is_member(B,A[i])){sub=false;break;}
		}
	}
return sub;
}
double similarity(Set A,Set B)
{
	Set D = set_intersection(A,B);
	Set E = set_union(A,B);
	return D.size()/((double)(E.size()));
}

				
