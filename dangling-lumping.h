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

// ofstream with higher precision to avoid truncation errors
struct my_ofstream : std::ofstream {
  explicit my_ofstream(std::streamsize prec = 15)
  {
    this->precision(prec);
  }
};

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

enum WriteMode { STATENODES, LINKS, CONTEXTS };

template <class T>
inline std::string to_string (const T& t){
	std::stringstream ss;
	ss << t;
	return ss.str();
}

class StateNode{
public:
	StateNode();
	StateNode(int stateid, int physid, double outweight);
	int stateId;
	int updatedStateId;
	int physId;
	int prevPhysId;
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
	prevPhysId = -1;
}

class PhysNode{
public:
	PhysNode();
	vector<int> stateNodeIndices;
	vector<int> stateNodeDanglingIndices;
	unordered_map<int,vector<int> > contextStateNodeIndices;

};

PhysNode::PhysNode(){
};


class StateNetwork{
private:
	double calcEntropyRate();
	bool readLines(string &line,vector<string> &lines);
	void writeLines(ifstream &ifs_tmp, ofstream &ofs, WriteMode &writeMode, string &line,int &batchNr);

	// For all batches
	string inFileName;
	string outFileName;
	string tmpOutFileName;
	mt19937 &mtRand;
	ifstream ifs;
  string line = "First line";
  double totWeight = 0.0;
  int updatedStateId = 0;
  double entropyRate = 0.0;
  unordered_map<int,int> completeStateNodeIdMapping;
  int totNphysNodes = 0;
	int totNstateNodes = 0;
	int totNlinks = 0;
	int totNdanglings = 0;
	int totNcontexts = 0;
	int totNphysDanglings = 0;

  // For each batch
  double weight = 0.0;
	int NphysNodes = 0;
	int NstateNodes = 0;
	int Nlinks = 0;
	int Ndanglings = 0;
	int Ncontexts = 0;
	int NphysDanglings = 0;
	unordered_map<int,int> stateNodeIdMapping;
	unordered_map<int,PhysNode> physNodes;
	unordered_map<int,StateNode> stateNodes;

public:
	StateNetwork(string infilename,string outfilename,mt19937 &mtrand);
	
	void lumpDanglings();
	bool loadStateNetworkBatch();
	void printStateNetworkBatch();
	void printStateNetwork();
	void concludeBatch();
	void compileBatches();

	bool keepReading = true;
  int Nbatches = 0;

};

StateNetwork::StateNetwork(string infilename,string outfilename,mt19937 &mtrand) : mtRand(mtrand){
	inFileName = infilename;
	outFileName = outfilename;
	tmpOutFileName = string(outFileName).append("_tmp");
	mtRand = mtrand;
  
  // Open state network
	ifs.open(inFileName.c_str());
	if(!ifs){
		cout << "failed to open \"" << inFileName << "\" exiting..." << endl;
		exit(-1);
	}

}

double StateNetwork::calcEntropyRate(){
	
	double h = 0.0;

	for(unordered_map<int,StateNode>::iterator it = stateNodes.begin(); it != stateNodes.end(); it++){
		StateNode &stateNode = it->second;
		if(stateNode.active){
			double H = 0.0;
			for(vector<pair<int,double> >::iterator it_link = stateNode.links.begin(); it_link != stateNode.links.end(); it_link++){
				double p = it_link->second/stateNode.outWeight;
				H -= p*log(p);
			}
			h += stateNode.outWeight*H/log(2.0);
		}
	}

	return h;

}

void StateNetwork::lumpDanglings(){

	unordered_set<int> physDanglings;
	int Nlumpings = 0;

	cout << "Lumping dangling state nodes:" << endl;

	// First loop records updated stateIds for lumped state nodes that other state nodes can be lumping to
	for(unordered_map<int,StateNode>::iterator it = stateNodes.begin(); it != stateNodes.end(); it++){
		StateNode &stateNode = it->second;
		PhysNode &physNode = physNodes[stateNode.physId];

		if(stateNode.outWeight > epsilon){
			// Record updated stateIds for non-dangling state nodes
			stateNodeIdMapping[stateNode.stateId] = updatedStateId;
			stateNode.updatedStateId = updatedStateId;
			updatedStateId++;
		}
		else{
			// Lump all dangling state nodes into one state node in dangling physical nodes, and update the stateIds
			int NnonDanglings = physNode.stateNodeIndices.size();
			if(NnonDanglings == 0){

				// When all state nodes are dangling, lump them to the first dangling state node id of the physical node
				physDanglings.insert(stateNode.physId);
				// Id of first dangling state node
				int lumpedStateIndex = physNode.stateNodeDanglingIndices[0];
				if(lumpedStateIndex == stateNode.stateId){
					// The first dangling state node in dangling physical node remains
					stateNodeIdMapping[stateNode.stateId] = updatedStateId;
					stateNode.updatedStateId = updatedStateId;
					updatedStateId++;
				}	
			}
		}
	}

	// Second loop records updated stateIds of lumping state nodes
	int NwithContext = 0;
	int NwithoutContext = 0;
	for(unordered_map<int,StateNode>::iterator it = stateNodes.begin(); it != stateNodes.end(); it++){
		StateNode &stateNode = it->second;
		PhysNode &physNode = physNodes[stateNode.physId];

		if(stateNode.outWeight < epsilon){

			int NnonDanglings = physNode.stateNodeIndices.size();
			
			if(NnonDanglings == 0){
	
				// When all state nodes are dangling, lump them to the first dangling state node id of the physical node
				// Id of first dangling state node
				int lumpedStateIndex = physNode.stateNodeDanglingIndices[0];
				if(lumpedStateIndex != stateNode.stateId){
					// All but the first dangling state node in dangling physical node are lumping to the first dangling state node
					// Add context to lumped state node
					stateNodes[lumpedStateIndex].contexts.insert(stateNodes[lumpedStateIndex].contexts.begin(),stateNode.contexts.begin(),stateNode.contexts.end());
					// Record updated state id to point to lumped state node with updated stateId and make it inactive
					stateNodeIdMapping[stateNode.stateId] = stateNodeIdMapping[stateNodes[lumpedStateIndex].stateId];
					stateNode.updatedStateId = stateNodes[lumpedStateIndex].stateId;
					stateNode.active = false;
					// Number of state nodes reduces by 1
					NstateNodes--;
					Nlumpings++;
				}	
			}
			else{
	
				// When dangling state node can be moved to non-dangling state node
				int lumpedStateIndex = -1;
				if(stateNode.prevPhysId >= 0){
					// Third order. First try context lumping.
					unordered_map<int,vector<int> >::iterator contextStates_it = physNode.contextStateNodeIndices.find(stateNode.prevPhysId);
					if(contextStates_it != physNode.contextStateNodeIndices.end()){
						int NcontextStates = contextStates_it->second.size();
						uniform_int_distribution<int> randInt(0,NcontextStates-1);
						// Find random state node with shared context
						lumpedStateIndex = contextStates_it->second[randInt(mtRand)];
						NwithContext++;
					}
				}
				if(lumpedStateIndex < 0){
					// If no shared context withing physical node
					uniform_int_distribution<int> randInt(0,NnonDanglings-1);
					// Find random state node
					lumpedStateIndex = physNode.stateNodeIndices[randInt(mtRand)];
					NwithoutContext++;
				}
				// Add context to lumped state node
				stateNodes[lumpedStateIndex].contexts.insert(stateNodes[lumpedStateIndex].contexts.begin(),stateNode.contexts.begin(),stateNode.contexts.end());
				
				// Update state id to point to lumped state node and make it inactive
				stateNodeIdMapping[stateNode.stateId] = stateNodeIdMapping[stateNodes[lumpedStateIndex].stateId];
				stateNode.updatedStateId = stateNodes[lumpedStateIndex].stateId;

				stateNode.active = false;
				// Number of state nodes reduces by 1
				NstateNodes--;
				Nlumpings++;
	
			}
		}
	}

	NphysDanglings = physDanglings.size();
	cout << "-->Lumped " << Nlumpings << " dangling state nodes (" << NwithContext << " with second-order context and " << NwithoutContext << " with first-order context)." << endl;
	cout << "-->Found " << NphysDanglings << " dangling physical nodes. Lumped dangling state nodes into a single dangling state node." << endl;
}

bool StateNetwork::readLines(string &line,vector<string> &lines){
	
	while(getline(ifs,line)){
		if(line[0] == '*'){
			return true;
		}
		else if(line[0] != '=' && line[0] != '#'){
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

	// ************************* Read statenetwork batch ************************* //
	
	// Read until next data label. Return false if no more data labels
	if(keepReading){
		cout << "Reading statenetwork, batch " << Nbatches+1 << ":" << endl;
		if(line[0] != '*'){
			while(getline(ifs,line)){
				if(line[0] == '*')
					break;
			}
		}
	}
	else{
		cout << "-->No more statenetwork batches to read." << endl;
		return false;
	}

	while(!readStates || !readLinks || !readContexts){

		ss.clear();
		ss.str(line);
		ss >> buf;
		if(!readStates && buf == "*States"){
			cout << "-->Reading states..." << flush;
			readStates = true;
			keepReading = readLines(line,stateLines);
			NstateNodes = stateLines.size();
			cout << "found " << NstateNodes << " states." << endl;
		}
		else if(!readLinks && buf == "*Links"){
			cout << "-->Reading links..." << flush;
			readLinks = true;
			keepReading = readLines(line,linkLines);
			Nlinks = linkLines.size();
			cout << "found " << Nlinks << " links." << endl;
		}
		else if(!readContexts && buf == "*Contexts"){
			cout << "-->Reading contexts..." << flush;
			readContexts = true;
			keepReading = readLines(line,contextLines);
			Ncontexts = contextLines.size();
			cout << "found " << Ncontexts << " contexts." << endl;
		}
		else{
			cout << "Expected *States, *Links, or *Contexts, but found " << buf << " exiting..." << endl;
			exit(-1);
		}
	}

	// ************************* Process statenetwork batch ************************* //
	Nbatches++;
	cout << "Processing statenetwork, batch " << Nbatches << ":" << endl;

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
	  weight += outWeight;
		if(outWeight > epsilon)
			physNodes[physId].stateNodeIndices.push_back(stateId);
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
	string delim = " ";
	for(int i=0;i<Ncontexts;i++){
			vector<string> contextLine = tokenize(contextLines[i],delim);
			int stateId = atoi(contextLine[0].c_str());
			if(contextLine.size() > 3){
				// Third order. Save context for context lumping.
				int physId = atoi(contextLine[1].c_str());
				int prevPhysId = atoi(contextLine[2].c_str());	
				// Add non-dangling state node to lumping context
				StateNode &stateNode = stateNodes[stateId];
				stateNode.prevPhysId = prevPhysId; 
				if(stateNode.outWeight > epsilon)
					physNodes[physId].contextStateNodeIndices[prevPhysId].push_back(stateId);
			}
			
			string context = contextLines[i].substr(contextLine[0].length()+1);
			stateNodes[stateId].contexts.push_back(context);
	}
	cout << "done!" << endl;

	// // Validate out-weights
	// for(unordered_map<int,StateNode>::iterator it = stateNodes.begin(); it != stateNodes.end(); it++){
	// StateNode &stateNode = it->second;
	// 	double w = 0.0;
	// 	for(vector<pair<int,double> >::iterator it_link = stateNode.links.begin(); it_link != stateNode.links.end(); it_link++){
	// 		w += it_link->second;
	// 	}
	// 	if((w < (stateNode.outWeight-epsilon)) || (w > (stateNode.outWeight+epsilon))){
	// 		cout << "::::::::::: Warning: out-weight does not match link weights for state node " << stateNode.stateId << ": " << stateNode.outWeight << " vs " << w << ", updating. :::::::::::" << endl;
	// 		stateNode.outWeight = w;
	// 	}
	// }

 	return true;

}

void StateNetwork::printStateNetworkBatch(){

  my_ofstream ofs;
	if(Nbatches == 1){ // Start with empty file for first batch
		ofs.open(tmpOutFileName.c_str());
	}
	else{ // Append to existing file
		ofs.open(tmpOutFileName.c_str(),ofstream::app);
	}
	cout << "Writing temporary results to " << tmpOutFileName << ":" << endl;

	cout << "-->Writing " << NstateNodes << " state nodes..." << flush;
	// To order state nodes by id
	map<int,int> orderedStateNodeIds;
	for(unordered_map<int,StateNode>::iterator it = stateNodes.begin(); it != stateNodes.end(); it++){
		StateNode &stateNode = it->second;
		if(stateNode.active)
			orderedStateNodeIds[stateNode.updatedStateId] = stateNode.stateId;
	}

	ofs << "*States\n";
	ofs << "#stateId ==> (physicalId, outWeight)\n";
	for(map<int,int>::iterator it = orderedStateNodeIds.begin(); it != orderedStateNodeIds.end(); it++){
		StateNode &stateNode = stateNodes[it->second];
		ofs << stateNode.stateId << " " << stateNode.physId << " " << stateNode.outWeight << "\n";
	}
	cout << "done!" << endl;

	cout << "-->Writing " << Nlinks << " links..." << flush;
	ofs << "*Links\n";
	ofs << "#(source target) ==> weight\n";
	for(map<int,int>::iterator it = orderedStateNodeIds.begin(); it != orderedStateNodeIds.end(); it++){
		StateNode &stateNode = stateNodes[it->second];
		for(vector<pair<int,double> >::iterator it_link = stateNode.links.begin(); it_link != stateNode.links.end(); it_link++){
				ofs << stateNode.stateId << " " << it_link->first << " " << it_link->second << "\n";
		}
	}
	cout << "done!" << endl;

	cout << "-->Writing " << Ncontexts << " contexts..." << flush;
	ofs << "*Contexts \n";
	ofs << "#stateId <== (physicalId priorId [history...])\n";
	for(map<int,int>::iterator it = orderedStateNodeIds.begin(); it != orderedStateNodeIds.end(); it++){
		StateNode &stateNode = stateNodes[it->second];
		for(vector<string>::iterator it_context = stateNode.contexts.begin(); it_context != stateNode.contexts.end(); it_context++){
			ofs << stateNode.stateId << " " << (*it_context) << "\n";
		}
	}
	cout << "done!" << endl;

}

void StateNetwork::printStateNetwork(){

	entropyRate += calcEntropyRate();

  my_ofstream ofs;
  ofs.open(outFileName.c_str());
 
	cout << "No more batches, writing results to " << outFileName << ":" << endl;
	cout << "-->Writing header comments..." << flush;
  	ofs << "# Physical nodes: " << NphysNodes << "\n";
  	ofs << "# Dangling physical nodes: " << NphysDanglings << "\n";
  	ofs << "# State nodes: " << NstateNodes << "\n";
  	ofs << "# Links: " << Nlinks << "\n";
  	ofs << "# Contexts: " << Ncontexts << "\n";
  	ofs << "# Weight: " << weight << "\n";
  	ofs << "# Entropy rate: " << entropyRate/weight << "\n";	
	cout << "done!" << endl;

	cout << "-->Writing " << NstateNodes << " state nodes..." << flush;
	// To order state nodes by id
	map<int,int> orderedStateNodeIds;
	for(unordered_map<int,StateNode>::iterator it = stateNodes.begin(); it != stateNodes.end(); it++){
		StateNode &stateNode = it->second;
		if(stateNode.active)
			orderedStateNodeIds[stateNode.updatedStateId] = stateNode.stateId;
	}

	ofs << "*States\n";
	ofs << "#stateId ==> (physicalId, outWeight)\n";
	for(map<int,int>::iterator it = orderedStateNodeIds.begin(); it != orderedStateNodeIds.end(); it++){
		StateNode &stateNode = stateNodes[it->second];
		ofs << stateNodeIdMapping[stateNode.stateId] << " " << stateNode.physId << " " << stateNode.outWeight << "\n";
	}
	cout << "done!" << endl;

	cout << "-->Writing " << Nlinks << " links..." << flush;
	ofs << "*Links\n";
	ofs << "#(source target) ==> weight\n";
	for(map<int,int>::iterator it = orderedStateNodeIds.begin(); it != orderedStateNodeIds.end(); it++){
		StateNode &stateNode = stateNodes[it->second];
		int source = stateNodeIdMapping[stateNode.stateId];
		for(vector<pair<int,double> >::iterator it_link = stateNode.links.begin(); it_link != stateNode.links.end(); it_link++){
			ofs << source << " " << stateNodeIdMapping[it_link->first] << " " << it_link->second << "\n";
		}
	}
	cout << "done!" << endl;

	cout << "-->Writing " << Ncontexts << " contexts..." << flush;
	ofs << "*Contexts \n";
	ofs << "#stateId <== (physicalId priorId [history...])\n";
	for(map<int,int>::iterator it = orderedStateNodeIds.begin(); it != orderedStateNodeIds.end(); it++){
		StateNode &stateNode = stateNodes[it->second];
		for(vector<string>::iterator it_context = stateNode.contexts.begin(); it_context != stateNode.contexts.end(); it_context++){
			ofs << stateNodeIdMapping[stateNode.stateId] << " " << (*it_context) << "\n";
		}
	}
	cout << "done!" << endl;

}

void StateNetwork::concludeBatch(){

	cout << "Concluding batch:" << endl;

	entropyRate += calcEntropyRate();
	totWeight += weight;
	totNphysNodes += NphysNodes;
	totNstateNodes += NstateNodes;
	totNlinks += Nlinks;
	totNdanglings += Ndanglings;
	totNcontexts += Ncontexts;
	totNphysDanglings += NphysDanglings;
	weight = 0.0;
	NphysNodes = 0;
	NstateNodes = 0;
	Nlinks = 0;
	Ndanglings = 0;
	Ncontexts = 0;
	NphysDanglings = 0;

	cout << "-->Current estimate of the entropy rate: " << entropyRate/totWeight << endl;

	completeStateNodeIdMapping.insert(stateNodeIdMapping.begin(),stateNodeIdMapping.end());
	stateNodeIdMapping.clear();
	physNodes.clear();
	stateNodes.clear();

}

void StateNetwork::compileBatches(){

  ifstream ifs_tmp(tmpOutFileName.c_str());
  my_ofstream ofs;
  ofs.open(outFileName);
  string buf;
	istringstream ss;
	bool writeStates = false;
	bool writeLinks = false;
	bool writeContexts = false;
	int batchNr = 1;

	cout << "Writing final results to " << outFileName << ":" << endl;
  
  	cout << "-->Writing header comments..." << flush;
  	ofs << "# Physical nodes: " << totNphysNodes << "\n";
	ofs << "# Number of dangling physical nodes: " << totNphysDanglings << "\n";  
  	ofs << "# State nodes: " << totNstateNodes << "\n";
  	ofs << "# Links: " << totNlinks << "\n";
  	ofs << "# Contexts: " << totNcontexts << "\n";
  	ofs << "# Weight: " << totWeight << "\n";
  	ofs << "# Entropy rate: " << entropyRate/totWeight << "\n";
	cout << "done!" << endl;

	cout << "-->Relabeling and writing " << totNstateNodes << " state nodes, " << totNlinks << " links, and " << totNcontexts << " contexts:" << endl;
	// Copy lines directly until data format
	while(getline(ifs_tmp,line)){
		if(line[0] == '*'){
			break;	
		}
		ofs << line << "\n";
	}
	while(!ifs_tmp.eof()){

		if(!writeStates && !writeLinks && !writeContexts){
			cout << "-->Batch " << batchNr << "/" << Nbatches << endl;
		}
		ofs << line << "\n";
		ss.clear();
		ss.str(line);
		ss >> buf;
		if(buf == "*States"){
			cout << "-->Writing state nodes..." << flush;
			writeStates = true;
			WriteMode writeMode = STATENODES;
			writeLines(ifs_tmp,ofs,writeMode,line,batchNr);
		}
		else if(buf == "*Links"){
			cout << "-->Writing links..." << flush;
			writeLinks = true;
			WriteMode writeMode = LINKS;
			writeLines(ifs_tmp,ofs,writeMode,line,batchNr);
		}
		else if(buf == "*Contexts"){
			cout << "-->Writing contexts..." << flush;
			writeContexts = true;
			WriteMode writeMode = CONTEXTS;
			writeLines(ifs_tmp,ofs,writeMode,line,batchNr);
		}
		else{
			cout << "Failed on line: " << line << endl;
		}
		cout << "done!" << endl;
		if(writeStates && writeLinks && writeContexts){
			writeStates = false;
			writeLinks = false;
			writeContexts = false;
			batchNr++;
		}
	}

	remove( tmpOutFileName.c_str() );

}

void StateNetwork::writeLines(ifstream &ifs_tmp, ofstream &ofs, WriteMode &writeMode, string &line,int &batchNr){

	string buf;
	istringstream ss;

	while(getline(ifs_tmp,line)){
		if(line[0] != '*'){
			if(line[0] != '=' && line[0] != '#'){
				ss.clear();
				ss.str(line);
				ss >> buf;
				if(writeMode == STATENODES){
					int stateId = atoi(buf.c_str());
					ss >> buf;
					int physId = atoi(buf.c_str());
	 				ss >> buf;
	 				double outWeight = atof(buf.c_str());
					ofs << completeStateNodeIdMapping[stateId] << " " << physId << " " << outWeight << "\n";
				}
				else if(writeMode == LINKS){
					int source = atoi(buf.c_str());
					ss >> buf;
					int target = atoi(buf.c_str());
					ss >> buf;
					double linkWeight = atof(buf.c_str());
					ofs << completeStateNodeIdMapping[source] << " " << completeStateNodeIdMapping[target] << " " << linkWeight << "\n";					
				}
				else if(writeMode == CONTEXTS){
					int stateId = atoi(buf.c_str());
					string context = line.substr(buf.length()+1);
					ofs << completeStateNodeIdMapping[stateId] << " " << context << "\n";
				}
			}
			else{
				if(line[0] == '='){
					ofs << "=== " << batchNr << "/" << Nbatches << " ===\n";
				}
				else{
					ofs << line << "\n";
				}
			}
		}
		else{
			return;
		}
	}
}
