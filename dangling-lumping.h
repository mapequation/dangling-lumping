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
const double epsilon = 1e-15;

unsigned stou(char *s);

template <class T>
inline string to_string (const T& t){
	stringstream ss;
	ss << t;
	return ss.str();
}

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
private:
	void calcEntropyRate();
	bool readLines(ifstream &ifs,string &line,vector<string> &lines);
	string inFileName;
	string outFileName;
	mt19937 &mtRand;
	ifstream ifs;
  string line = "First line";
  bool keepReading = true;
  int Nbatches = 0;
	int NphysNodes = 0;
	int NstateNodes = 0;
	int Nlinks = 0;
	int Ndanglings = 0;
	int Ncontexts = 0;
	int NphysDanglings = 0;
	double totWeight = 0.0;
	double entropyRate = 0.0;
	unordered_map<int,PhysNode> physNodes;
	unordered_map<int,StateNode> stateNodes;

public:
	StateNetwork(string infilename,string outfilename,mt19937 &mtrand);
	
	void lumpDanglings();
	bool loadStateNetworkBatch();
	void printStateNetworkBatch();
	void compileBatches();

};

StateNetwork::StateNetwork(string infilename,string outfilename,mt19937 &mtrand) : mtRand(mtrand){
	inFileName = infilename;
	outFileName = outfilename;
	mtRand = mtrand;
  
  // Open state network
	ifs.open(inFileName.c_str());
	if(!ifs){
		cout << "failed to open \"" << inFileName << "\" exiting..." << endl;
		exit(-1);
	}

}

void StateNetwork::calcEntropyRate(){
	
	entropyRate = 0.0;

	for(unordered_map<int,StateNode>::iterator it = stateNodes.begin(); it != stateNodes.end(); it++){
	StateNode &stateNode = it->second;
		if(stateNode.active){
			double h = 0.0;
			for(vector<pair<int,double> >::iterator it_link = stateNode.links.begin(); it_link != stateNode.links.end(); it_link++){
				double p = it_link->second/stateNode.outWeight;
				h -= p*log(p);
			}
			entropyRate += stateNode.outWeight/totWeight*h/log(2.0);
		}
	}

}

// void StateNetwork::lumpDanglings(){

// 	unordered_set<int> physDanglings;
// 	int Nlumpings = 0;
	

// }

void StateNetwork::lumpDanglings(){

	unordered_set<int> physDanglings;
	int Nlumpings = 0;

	cout << "Lumping dangling state nodes:" << endl;

	// First loop sets updated stateIds of non-dangling state nodes and state nodes in dangling physical nodes, which are lumped into one state node
	int updatedStateId = 0;
	for(unordered_map<int,StateNode>::iterator it = stateNodes.begin(); it != stateNodes.end(); it++){
		StateNode &stateNode = it->second;

		if(stateNode.outWeight > epsilon){
			// Set updated stateIds for non-dangling state nodes
			stateNode.stateId = updatedStateId;
			updatedStateId++;
		}
		else{
			// Lump all dangling state nodes into one state node in dangling physical nodes, and update the stateIds
			int NnonDanglings = physNodes[stateNode.physId].stateNodeNonDanglingIndices.size();
			if(NnonDanglings == 0){

				// When all state nodes are dangling, lump them to the first dangling state node id of the physical node
				physDanglings.insert(stateNode.physId);
				// Id of first dangling state node
				int lumpedStateIndex = physNodes[stateNode.physId].stateNodeDanglingIndices[0];
				if(lumpedStateIndex == stateNode.stateId){
					// The first dangling state node in dangling physical node remains
					stateNode.stateId = updatedStateId;
					updatedStateId++;
				}	
				else{
					// Add context to lumped state node
					stateNodes[lumpedStateIndex].contexts.insert(stateNodes[lumpedStateIndex].contexts.begin(),stateNode.contexts.begin(),stateNode.contexts.end());
					// Update state id to point to lumped state node with upodated stateId and make it inactive
					stateNode.stateId = stateNodes[lumpedStateIndex].stateId;
					stateNode.active = false;
					// Number of state nodes reduces by 1
					NstateNodes--;
					Nlumpings++;
				}	
			}
		}
	}

	// Second loop sets updated stateIds of dangling state nodes in physical nodes with non-dangling state nodes
	for(unordered_map<int,StateNode>::iterator it = stateNodes.begin(); it != stateNodes.end(); it++){
		StateNode &stateNode = it->second;

		if(stateNode.outWeight < epsilon){
			int NnonDanglings = physNodes[stateNode.physId].stateNodeNonDanglingIndices.size();
			if(NnonDanglings > 0){

				uniform_int_distribution<int> randInt(0,NnonDanglings-1);
				// Find random state node
				int lumpedStateIndex = physNodes[stateNode.physId].stateNodeNonDanglingIndices[randInt(mtRand)];
				// Add context to lumped state node
				stateNodes[lumpedStateIndex].contexts.insert(stateNodes[lumpedStateIndex].contexts.begin(),stateNode.contexts.begin(),stateNode.contexts.end());
				
				// Update state id to point to lumped state node and make it inactive
				stateNode.stateId = stateNodes[lumpedStateIndex].stateId;

				stateNode.active = false;
				// Number of state nodes reduces by 1
				NstateNodes--;
				Nlumpings++;

			}
		}
	}

	NphysDanglings = physDanglings.size();
	cout << "-->Lumped " << Nlumpings << " dangling state nodes." << endl;
	cout << "-->Found " << NphysDanglings << " dangling physical nodes. Lumped dangling state nodes into a single dangling state node." << endl;

}

bool StateNetwork::readLines(ifstream &ifs,string &line,vector<string> &lines){
	
	while(getline(ifs,line)){
		if(line[0] == '*'){
			return true;
		}
		if(line[0] != '-' && line[0] != '#'){
			lines.push_back(line);
		}
	}

	return false; // Reached end of file
}

bool StateNetwork::loadStateNetworkBatch(){

	vector<string> stateLines;
	vector<string> linkLines;
	vector<string> contextLines;
	bool readStates = false;
	bool readLinks = false;
	bool readContexts = false;
	string buf;
	istringstream ss;

	// ************************* Read statefile batch ************************* //
	
	// Read until next data label. Return false if no more data labels
	if(keepReading){
		cout << "Reading statefile batch " << Nbatches+1 << ":" << endl;
		if(line[0] != '*'){
			while(getline(ifs,line)){
				if(line[0] == '*')
					break;
			}
		}
	}
	else{
		cout << "No more statefile batches to read." << endl;
		return false;
	}

	while(!readStates || !readLinks || !readContexts){

		ss.clear();
		ss.str(line);
		ss >> buf;
		if(!readStates && buf == "*States"){
			cout << "-->Reading states..." << flush;
			readStates = true;
			keepReading = readLines(ifs,line,stateLines);
			NstateNodes = stateLines.size();
			cout << "found " << NstateNodes << " states." << endl;
		}
		else if(!readLinks && buf == "*Links"){
			cout << "-->Reading links..." << flush;
			readLinks = true;
			keepReading = readLines(ifs,line,linkLines);
			Nlinks = linkLines.size();
			cout << "found " << Nlinks << " links." << endl;
		}
		else if(!readContexts && buf == "*Contexts"){
			cout << "-->Reading contexts..." << flush;
			readContexts = true;
			keepReading = readLines(ifs,line,contextLines);
			Ncontexts = contextLines.size();
			cout << "found " << Ncontexts << " contexts." << endl;
		}
		else{
			cout << "Expected *States, *Links, or *Contexts, but found " << buf << " exiting..." << endl;
			exit(-1);
		}
	}

	// ************************* Process statefile batch ************************* //
	Nbatches++;
	cout << "Processing statefile batch " << Nbatches << ":" << endl;

	//Process states
	cout << "-->Processing " << NstateNodes  << " state nodes..." << flush;
	for(int i=0;i<NstateNodes;i++){

		ss.clear();
		ss.str(stateLines[i]);
		ss >> buf;
		int stateId = atoi(buf.c_str());
		ss >> buf;
		int physId = atoi(buf.c_str());
	  ss >> buf;
	  double outWeight = atof(buf.c_str());
	  totWeight += outWeight;
		if(outWeight > epsilon)
			physNodes[physId].stateNodeNonDanglingIndices.push_back(stateId);
		else{
			physNodes[physId].stateNodeDanglingIndices.push_back(stateId);
			Ndanglings++;
		}
		stateNodes[stateId] = StateNode(stateId,physId,outWeight);
	}
	NphysNodes = physNodes.size();
	cout << "found " << Ndanglings << " dangling state nodes in " << NphysNodes << " physical nodes, done!" << endl;

	// Process links 
	cout << "-->Processing " << Nlinks  << " links..." << flush;
	for(int i=0;i<Nlinks;i++){
			ss.clear();
			ss.str(linkLines[i]);
			ss >> buf;
			int source = atoi(buf.c_str());
			ss >> buf;
			int target = atoi(buf.c_str());
			ss >> buf;
			double linkWeight = atof(buf.c_str());
			stateNodes[source].links.push_back(make_pair(target,linkWeight));
	}
 	cout << "done!" << endl;

	// Process contexts
	cout << "-->Processing " << Ncontexts  << " contexts..." << flush;
	for(int i=0;i<Ncontexts;i++){
			ss.clear();
			ss.str(contextLines[i]);
			ss >> buf;
			int stateNodeId = atoi(buf.c_str());
			string context = contextLines[i].substr(buf.length()+1);
			stateNodes[stateNodeId].contexts.push_back(context);
	}
	cout << "done!" << endl;

 	return true;

}

void StateNetwork::printStateNetworkBatch(){

  calcEntropyRate();

	cout << "Writing to " << outFileName << ":" << endl;
  ofstream ofs(outFileName);
  cout << "-->Writing header comments..." << flush;
  ofs << "# Number of physical nodes: " << NphysNodes << "\n";
  ofs << "# Number of state nodes: " << NstateNodes << "\n";
  ofs << "# Number of dangling physical (and state) nodes: " << NphysDanglings << "\n";
  ofs << "# Number of links: " << Nlinks << "\n";
  ofs << "# Number of contexts: " << Ncontexts << "\n";
  ofs << "# Total weight: " << totWeight << "\n";
  ofs << "# Entropy rate: " << entropyRate << "\n";
	cout << "done!" << endl;

	cout << "-->Writing " << NstateNodes << " state nodes..." << flush;
	ofs << "*States " << NstateNodes << "\n";
	ofs << "#stateId ==> (physicalId, outWeight)\n";
	int index = 0;
	for(unordered_map<int,StateNode>::iterator it = stateNodes.begin(); it != stateNodes.end(); it++){
		StateNode &stateNode = it->second;
		if(stateNode.active){
			// The state node has not been lumped to another node (but other nodes may have been lumped to it)
			ofs << stateNode.stateId << " " << stateNode.physId << " " << stateNode.outWeight << "\n";
		}
		index++;
	}
	cout << "done!" << endl;

	cout << "-->Writing " << Nlinks << " links..." << flush;
	ofs << "*Links " << Nlinks << "\n";
	ofs << "#(source target) ==> weight\n";
	for(unordered_map<int,StateNode>::iterator it = stateNodes.begin(); it != stateNodes.end(); it++){
		StateNode &stateNode = it->second;
		if(stateNode.active){
			// The state node has not been lumped to another node (but other nodes may have been lumped to it)
			for(vector<pair<int,double> >::iterator it_link = stateNode.links.begin(); it_link != stateNode.links.end(); it_link++){
				ofs << stateNode.stateId << " " << stateNodes[it_link->first].stateId << " " << it_link->second << "\n";
			}
		}
	}
	cout << "done!" << endl;

	cout << "-->Writing " << Ncontexts << " contexts..." << flush;
	ofs << "*Contexts \n";
	ofs << "#stateId <== (physicalId priorId [history...])\n";
	index = 0;
	for(unordered_map<int,StateNode>::iterator it = stateNodes.begin(); it != stateNodes.end(); it++){
		StateNode &stateNode = it->second;
		if(stateNode.active){
		// The state node has not been lumped to another node (but other nodes may have been lumped to it)
			for(vector<string>::iterator it_context = stateNode.contexts.begin(); it_context != stateNode.contexts.end(); it_context++){
				ofs << stateNode.stateId << " " << (*it_context) << "\n";
			}
		}
		index++;
	}
	cout << "done!" << endl;

}

void StateNetwork::compileBatches(){

}



