#pragma once


#include<iostream>
#include<fstream>
#include<vector>
#include<unordered_set>
#include<unordered_map>
#include<algorithm>
#include<functional>
#include<math.h>
#include<string>
#include<algorithm>
#include<sys/time.h>
#include<sys/stat.h>
#include<sys/types.h>

#include<sys/resource.h>
#include<time.h>
#include<chrono>

using namespace std;

#include <set>



//Taking as ordered_set. set bcoz duplicates are not allowed. ordered bcoz I need in Increasing order.
//set<char> global_vrtxlbl_set;
//unordered_map<char, unsigned> global_vrtxlbl_map;  //to map global vertex label to unsigned numeric. 

class Graph
{
	public:
		unsigned gid, vertexCount, edgeCount; // graph-id, no. of vertices, no. of edges
		vector<unsigned> vertices; // Vertex-set
		vector<unsigned> degrees; // Degree-set
		vector<char> vrtxLbl;      //will save vertex label of all vertices
		unordered_map<unsigned,unsigned> vid_to_ind; // vid to index in adjacency list of graph
		unordered_map<unsigned,unsigned> vid_to_vc;
		vector< unsigned> edges; // Edge-Set
	        vector<unsigned> edgeType;  //------>	
		//unordered_map<unsigned, unsigned> VertexLabelMap; // count of a label of vertex label
		vector<unsigned> VertexLabelVectr; //will save count of a vertex label in the graph.
		//unordered_map<unsigned, unsigned> EdgeLabelMap;  // count of a type of vertex label
		vector<unsigned> edgeTypeVectr;  //-----> // count the particular edge type in this graph.

		Graph(){
			gid = 0;
			vertexCount = 0;
			edgeCount = 0;
		}

		//void readGraph(istream &inp, unsigned vlblCount, unordered_set<unsigned>& v_label, unordered_set<unsigned>& e_label); // reads the graph from input file 
		void readGraph(istream &inp, unsigned vlblCount, unordered_map<char, unsigned> global_vrtxlbl_map, unordered_map<string, unsigned> global_edgeType_map, unordered_set<unsigned>& e_label);
		void pushEdge(unsigned u, unsigned v, unsigned ec, unordered_set<unsigned>& e_label); // adds an edge to the graph 
		void displayGraph(); // prints details of the graph

};
