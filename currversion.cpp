#include<iostream>
#include<bits/stdc++.h>
using namespace std;
class Graph
{
	public:
		int  gid, vertexCount, edgeCount ; // graph-id, no. of vertices, no. of edges
		//vector<unsigned> vertices; // Vertex-set
		vector<int> degrees{0,0,0,0,0};
	    int edgemat[5][5] ={{0,0,0,0,0} ,{0,0,0,0,0},{0,0,0,0,0},{0,0,0,0,0},{0,0,0,0,0}};
	    
		 // Degree-set
		vector< char > vrtxLbl; // a b c type
		
		      //will save vertex label of all vertices
	    map< int, int > vid_to_ind;  // id to index
	    
		 // vid to index in adjacency list of graph
		map< int ,int > vid_to_vc;  // id to label
		
		vector<pair<int,int>> edges; // Edge-Set
	    vector< int > edgeType;  //------>	
		
		
		vector< int > VertexLabelVectr; // will save count of a vertex label in the graph. c02-3=count
		
		vector< int > edgeTypeVectr;  //-----> // count the particular edge type in this graph. if edge is bw c and 0 , will be count as a type

		Graph()
		{
			gid = 0;
			vertexCount = 0;
			edgeCount = 0;
		}

		//void readGraph(istream &inp, unsigned vlblCount, unordered_set<unsigned>& v_label, unordered_set<unsigned>& e_label); // reads the graph from input file 
		void readGraph( int vlblCount, map< char, int > global_vrtxlbl_map, map<string, int> global_edgeType_map );
		// labelmap a - 1 , b - 2   etypemap c-o - 2
		void pushEdge ( int  u , int v); // adds an edge to the graph 
		void displayGraph(); // prints details of the graph
		void triangles(vector<vector<int>> & edgemat , vector<int> & edges );
	};
	void Graph::triangles(vector<vector<int>> & edgemat , vector<int> & edges)
	{
		
	}
	
	void Graph  :: pushEdge(int u, int  v)
{
	//if(u>v)
		edges.push_back(make_pair(u, v));
	//else 
		edges.push_back(make_pair(v, u));
	degrees[u]++;
	degrees[v]++;
	edgemat[u][v] = 1;
	edgemat[v][u] = 1;
}
	

	void Graph:: readGraph( int vlblCount ,   map<char, int> global_vrtxlbl_map,   map<string,  int> global_edgetype_map)
	{
	char tag; // uselss variable
	// first line should be of the format "g vertexCount(unsigned int) edgeCount(unsigned int) gid(unsigned int)"
	cin >> tag; // the tag 'g'
	cin >> vertexCount; // the no. of vertices in the graph
	cin >> edgeCount; // the no. of edges in the graph
	cin >> gid; // the graph-id of the graph
	//vertices.resize  ( vertexCount );
	//degrees.assign( vertexCount , 0);
	vrtxLbl.resize( vertexCount);
    edgeTypeVectr.resize( global_edgetype_map.size(), 0 );
	
	int ind = 0 ;
	int vid , src_vtx , dest_vtx , ec ;
	char vc ;
	//unsigned vc,sc,dc;  //Vertex type -- vc 
	for( int vtx_ind=0; vtx_ind < vertexCount; vtx_ind++ )
	{
		// each line for each vertex should be in the format like: "v vid(unsigned int)"
		cin >> tag >> vid >> vc; // the tag 'v' along with the vertex-id
		// vertices[vtx_ind] = (vc);//vid;
		vid_to_ind[vid] = vtx_ind; // mapping vertex-id to its index
		vid_to_vc[vid] = vc;
		vrtxLbl[vtx_ind] = vc;
	}
	for(int e_ind = 0; e_ind < edgeCount; e_ind++)
	{
		//each line for each edge should be in the format like: "e vid_src(unsigned int) vid_dest(unsigned int)"
		cin >> tag >> src_vtx  >> dest_vtx ; // the tag 'e' along with the source and destination vertex-ids
		// Undirected graph : adding edge source to destination and destination to source
		/*if( char(vid_to_vc[src_vtx]-65 +'A') < char(vid_to_vc[dest_vtx]- 65 + 'A' ))
                {
                        string eType {char(vid_to_vc[src_vtx]-65 +'A') , '-', char(vid_to_vc[dest_vtx]- 65 + 'A') };
                        int unsign_eType = global_edgetype_map[eType];
                        //cout<<eType<<"  "<<unsign_eType<<"    "<<edgeTypeVectr[unsign_eType]<<endl;
                        edgeTypeVectr[unsign_eType]++;
                        //cout<<edgeTypeVectr[unsign_eType]<<endl;
                }
                else
                {
                        string eType {char(vid_to_vc[dest_vtx]-65 +'A') , '-', char(vid_to_vc[src_vtx]- 65 + 'A') };
                        //edgeType[e_ind] = global_edgetype_map[eType];
                        int unsign_eType = global_edgetype_map[eType];
                        //cout<<eType<<"  "<<unsign_eType<<"    "<<edgeTypeVectr[unsign_eType]<<endl;
                        edgeTypeVectr[unsign_eType]++;
                        //cout<<edgeTypeVectr[unsign_eType]<<endl;
                }*/

	//	pushEdge(vid_to_vc[src_vtx], vid_to_vc[dest_vtx], ec);
		pushEdge( src_vtx, dest_vtx);
		//pushEdge(sc,dc);
	}

	VertexLabelVectr.resize(vlblCount, 0);
        for(int i=0; i<vrtxLbl.size(); i++)
        {
            int vlbl = global_vrtxlbl_map[vrtxLbl[i]];
            VertexLabelVectr[vlbl]++;
        }

}
int main()
{
	cout<<8;
	Graph g1 ;
	map< char, int > global_vrtxlbl_map ;
    map<string, int> global_edgeType_map; 
	g1.readGraph(5, global_vrtxlbl_map ,global_edgeType_map );
	cout<<g1.degrees[0];
	for(auto i : g1.edges)
	{
		cout<<i.first<<i.second<<" ";
	}
	int count=0;
	for(auto i: g1.edges)
	{
		if(i.first<i.second)
		{
			for(int j =i.second;j<5;j++)
			{
				if(g1.edgemat[i.first][j] && g1.edgemat[i.second][j])
				{
				   count++;
				}
			}
			
		}
	}
	cout<<count;
	
	
	return 0 ;
}
