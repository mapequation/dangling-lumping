#include <cstdio>
#include <cmath>
#include <cstdlib>
#include <cstring>
#include <ctime>
#include <iostream>
#include <sstream>
#include <fstream>
#include <iomanip>
#include <random>
#include <functional>
#include <map>
#include <unordered_map>
#include <unordered_set>
using namespace std;
const double epsilon = 1e-10;

unsigned stou(char *s);

vector<string> tokenize(const string& str,string& delimiters){

	vector<string> tokens;

  // skip delimiters at beginning.
	string::size_type lastPos = str.find_first_not_of(delimiters, 0);

  // find first "non-delimiter".
	string::size_type pos = str.find_first_of(delimiters, lastPos);

	while(string::npos != pos || string::npos != lastPos){

    // found a token, add it to the vector.
		tokens.push_back(str.substr(lastPos, pos - lastPos));

    // skip delimiters.  Note the "not_of"
		lastPos = str.find_first_not_of(delimiters, pos);

    // find next "non-delimiter"
		pos = str.find_first_of(delimiters, lastPos);
	}

	return tokens;

}


class StateNode{
public:
	StateNode();
	StateNode(int stateid, int physid, double outweight);
	int stateId;
	int physId;
	double outWeight;
	bool active = true;
	vector<pair<int,double> > links;
	vector<string> contexts;
};

StateNode::StateNode(){
};

StateNode::StateNode(int stateid, int physid, double outweight){
	stateId = stateid;
	physId = physid;
	outWeight = outweight;
}

class PhysNode{
public:
	PhysNode();
	vector<int> stateNodeNonDanglingIndices;
	vector<int> stateNodeDanglingIndices;
};

PhysNode::PhysNode(){
};


class StateNetwork{
public:
	StateNetwork(string infilename,string outfilename);
	string inFileName;
	string outFileName;
	int NphysNodes = 0;
	int NstateNodes = 0;
	int Nlinks = 0;
	int Ndanglings = 0;
	int Ncontexts = 0;
	int NphysDanglings = 0;
	double totWeight = 0.0;
	double entropyRate = 0.0;
	unordered_map<int,PhysNode> physNodes;
	vector<StateNode> stateNodes;
};

StateNetwork::StateNetwork(string infilename,string outfilename){
	inFileName = infilename;
	outFileName = outfilename;
}

template <class T>
inline std::string to_string (const T& t){
	std::stringstream ss;
	ss << t;
	return ss.str();
}

void calcEntropyRate(StateNetwork &statenetwork){

	for(vector<StateNode>::iterator it = statenetwork.stateNodes.begin(); it != statenetwork.stateNodes.end(); it++){
		if(it->active){
			double h = 0.0;
			for(vector<pair<int,double> >::iterator it_link = it->links.begin(); it_link != it->links.end(); it_link++){
				double p = it_link->second/it->outWeight;
				h -= p*log(p);
			}
			statenetwork.entropyRate += it->outWeight/statenetwork.totWeight*h/log(2.0);
		}
	}

}

void lumpDanglings(StateNetwork &statenetwork,std::mt19937 &mtrand){

	unordered_set<int> physDanglings;
	int Nlumpings = 0;

	cout << "Lumping dangling state nodes:" << endl;
	int updatedStateId = 0;
	for(vector<StateNode>::iterator it = statenetwork.stateNodes.begin(); it != statenetwork.stateNodes.end(); it++){
		if(it->outWeight > epsilon){
			it->stateId = updatedStateId;
			updatedStateId++;
		}
		else{
			// Lump with non-dangling state node in the same physical node
			int NnonDanglings = statenetwork.physNodes[it->physId].stateNodeNonDanglingIndices.size();
			if(NnonDanglings > 0){
				std::uniform_int_distribution<int> randInt(0,NnonDanglings-1);
				// Find random state node
				int lumpedStateId = statenetwork.physNodes[it->physId].stateNodeNonDanglingIndices[randInt(mtrand)];
				// Add context to lumped state node
				statenetwork.stateNodes[lumpedStateId].contexts.insert(statenetwork.stateNodes[lumpedStateId].contexts.begin(),it->contexts.begin(),it->contexts.end());
				
				// Update state id to point to lumped state node and make it inactive
				it->stateId = lumpedStateId;
				it->active = false;
				// Number of state nodes reduces by 1
				statenetwork.NstateNodes--;
				Nlumpings++;

			}
			else{
				// When all state nodes are dangling, lump them to the first dangling state node id of the physical node
				physDanglings.insert(it->physId);
				// Id of first dangling state node
				int lumpedStateId = statenetwork.physNodes[it->physId].stateNodeDanglingIndices[0];
				if(lumpedStateId == it->stateId){
					// The first dangling state node remains
					it->stateId = updatedStateId;
					updatedStateId++;
				}	
				else{
					// Add context to lumped state node
					statenetwork.stateNodes[lumpedStateId].contexts.insert(statenetwork.stateNodes[lumpedStateId].contexts.begin(),it->contexts.begin(),it->contexts.end());
					// Update state id to point to lumped state node and make it inactive
					it->stateId = lumpedStateId;
					it->active = false;
					// Number of state nodes reduces by 1
					statenetwork.NstateNodes--;
					Nlumpings++;
				}	
			}
		}
	}
	statenetwork.NphysDanglings = physDanglings.size();
	cout << "-->Lumped " << Nlumpings << " dangling state nodes." << endl;
	cout << "-->Found " << statenetwork.NphysDanglings << " physical nodes without non-dangling state nodes. Lumped into a single dangling state node in each physical node." << endl;
}

void loadStateNetwork(StateNetwork &statenetwork){

	string line;
	string buf;
	istringstream ss;

  // ************************* Read state network ************************* //
	ifstream net(statenetwork.inFileName.c_str());
	if(!net){
		cout << "failed to open \"" << statenetwork.inFileName << "\" exiting..." << endl;
		exit(-1);
	}
	else{
		cout << "Reading " << statenetwork.inFileName << ":" << endl;
	}

	// Skip header lines starting with #
	while(getline(net,line)){
		if(line[0] != '#')
			break;
	}

	// ************************* Read states ************************* //
	ss.clear();
	ss.str(line);
	ss >> buf;
	if(buf != "*States"){
		cout << "Expected *States but read " << buf << ", exiting..." << endl;
		exit(-1);
	}
	ss >> buf;
	statenetwork.NstateNodes = atoi(buf.c_str());
	cout << "-->Reading " << statenetwork.NstateNodes  << " state nodes..." << flush;
	statenetwork.stateNodes = vector<StateNode>(statenetwork.NstateNodes);
	for(int i=0;i<statenetwork.NstateNodes;i++){
		getline(net,line);
		if(line[0] != '#'){
			ss.clear();
			ss.str(line);
			ss >> buf;
			int stateId = atoi(buf.c_str());
			ss >> buf;
			int physId = atoi(buf.c_str());
	  	ss >> buf;
	  	double outWeight = atof(buf.c_str());
	  	statenetwork.totWeight += outWeight;
			if(outWeight > epsilon)
				statenetwork.physNodes[physId].stateNodeNonDanglingIndices.push_back(stateId);
			else{
				statenetwork.physNodes[physId].stateNodeDanglingIndices.push_back(stateId);
				statenetwork.Ndanglings++;
			}
			statenetwork.stateNodes[stateId] = StateNode(stateId,physId,outWeight);
		}
		else{
			// One extra step for each # comment.
			i--;
		}
	}
	statenetwork.NphysNodes = statenetwork.physNodes.size();
	cout << "found " << statenetwork.Ndanglings << " dangling state nodes in " << statenetwork.NphysNodes << " physical nodes, done!" << endl;

	// ************************* Read arcs ************************* //
	getline(net,line);
	ss.clear();
	ss.str(line);
	ss >> buf;
	if(buf != "*Arcs"){
		cout << "Expected *Arcs but read " << buf << ", exiting..." << endl;
		exit(-1);
	}
	ss >> buf;
	statenetwork.Nlinks = atoi(buf.c_str());
	cout << "-->Reading " << statenetwork.Nlinks  << " state links..." << flush;


	for(int i=0;i<statenetwork.Nlinks;i++){
		getline(net,line);
		if(line[0] != '#'){
			ss.clear();
			ss.str(line);
			ss >> buf;
			int source = atoi(buf.c_str());
			ss >> buf;
			int target = atoi(buf.c_str());
			ss >> buf;
			double linkWeight = atof(buf.c_str());
			statenetwork.stateNodes[source].links.push_back(make_pair(target,linkWeight));
 		}
 		else{
 			// One extra step for each # comment.
 			i--;
 		}
	}
 	cout << "done!" << endl;

	// ************************* Read contexts ************************* //
	getline(net,line);
	ss.clear();
	ss.str(line);
	ss >> buf;
	if(buf != "*Contexts" && buf != "*MemoryNodes"){
		cout << "Expected *Contexts or *MemoryNodes but read " << buf << ", exiting..." << endl;
		exit(-1);
	}
	cout << "-->Reading state node contexts..." << flush;
	while(getline(net,line)){
		if(line[0] != '#'){
			ss.clear();
			ss.str(line);
			ss >> buf;
			int stateNodeId = atoi(buf.c_str());
			string context = line.substr(buf.length()+1);
			statenetwork.stateNodes[stateNodeId].contexts.push_back(context);
			statenetwork.Ncontexts++;
 		}
	}
 	cout << "found " << statenetwork.Ncontexts << ", done!" << endl;

	net.close();

}

void printStateNetwork(StateNetwork &statenetwork){
	cout << "Writing to " << statenetwork.outFileName << ":" << endl;
  ofstream outfile(statenetwork.outFileName);
  cout << "-->Writing header comments..." << flush;
  outfile << "# Number of physical nodes: " << statenetwork.NphysNodes << "\n";
  outfile << "# Number of state nodes: " << statenetwork.NstateNodes << "\n";
  outfile << "# Number of dangling physical (and state) nodes: " << statenetwork.NphysDanglings << "\n";
  outfile << "# Number of links: " << statenetwork.Nlinks << "\n";
  outfile << "# Total weight: " << statenetwork.totWeight << "\n";
  outfile << "# Entropy rate: " << statenetwork.entropyRate << "\n";
	cout << "done!" << endl;

	cout << "-->Writing state nodes..." << flush;
	outfile << "*States " << statenetwork.NstateNodes << "\n";
	int index = 0;
	for(vector<StateNode>::iterator it = statenetwork.stateNodes.begin(); it != statenetwork.stateNodes.end(); it++){
		if(it->active){
			// The state node has not been lumped to another node (but other nodes may have been lumped to it)
			outfile << it->stateId << " " << it->physId << " " << it->outWeight << "\n";
		}
		index++;
	}
	cout << "done!" << endl;

	cout << "-->Writing links..." << flush;
	outfile << "*Arcs " << statenetwork.Nlinks << "\n";
	for(vector<StateNode>::iterator it = statenetwork.stateNodes.begin(); it != statenetwork.stateNodes.end(); it++){
		if(it->active){
			// The state node has not been lumped to another node (but other nodes may have been lumped to it)
			for(vector<pair<int,double> >::iterator it_link = it->links.begin(); it_link != it->links.end(); it_link++){
				outfile << it->stateId << " " << statenetwork.stateNodes[it_link->first].stateId << " " << it_link->second << "\n";
			}
		}
	}
	cout << "done!" << endl;

	cout << "-->Writing contexts..." << flush;
	outfile << "*Contexts \n";
	index = 0;
	for(vector<StateNode>::iterator it = statenetwork.stateNodes.begin(); it != statenetwork.stateNodes.end(); it++){
		if(it->active){
		// The state node has not been lumped to another node (but other nodes may have been lumped to it)
			for(vector<string>::iterator it_context = it->contexts.begin(); it_context != it->contexts.end(); it_context++){
				outfile << it->stateId << " " << (*it_context) << "\n";
			}
		}
		index++;
	}
	cout << "done!" << endl;

	outfile.close();
}



