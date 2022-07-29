#include "veo.h"
#include <cmath>
#include <chrono>




extern unordered_map<char, unsigned> global_vrtxlbl_map;  //to map global vertex label to unsigned numeric. 
extern unordered_map<string, unsigned> global_edgetype_map;
unsigned long r;



bool freqComp(pair<pair<unsigned, unsigned>,unsigned> v1, pair<pair<unsigned, unsigned>,unsigned> v2)
{ 
    return (v1.second > v2.second); 
}


bool rankComp(pair<unsigned, unsigned > g_r1, pair<unsigned, unsigned > g_r2)
{
    return (g_r1.first > g_r2.first);
}





// Edge comparator
bool edgeComp(pair<unsigned, unsigned> &a, pair<unsigned, unsigned> &b)
{
	if(a.first == b.first)
		return a.second < b.second;
	return a.first < b.first;
}

// Intersection of vertex lists of 2 graphs
/*
double intersection_vertices(vector<unsigned> &s1, vector<unsigned> &s2)
{
	unsigned s1_iter = 0;
	unsigned s2_iter = 0;
	double common = 0;

	while(s1_iter < s1.size() && s2_iter < s2.size())
	{
		if(s1[s1_iter] < s2[s2_iter])
			s1_iter++;
		else if(s1[s1_iter] > s2[s2_iter])
			s2_iter++;
		else
		{
			
			s1_iter++;
			s2_iter++;
			common++;		
		}
	}
	return common;
}



// Intersection of edge lists of 2 graphs
double intersection_edges(vector<pair<unsigned, unsigned>> &s1, vector<pair<unsigned, unsigned>> &s2)
{
	unsigned s1_iter = 0;
	unsigned s2_iter = 0;
	double common = 0;

	for(auto s1_iter = s1.begin(), s2_iter = s2.begin(); s1_iter != s1.end() && s2_iter != s2.end(); )
	{
		if(edgeComp(*s1_iter,*s2_iter))
			s1_iter++;
		else if(edgeComp(*s2_iter,*s1_iter))
			s2_iter++;
		else
		{
			
			s1_iter++;
			s2_iter++;
			common++;		
		}
	}
	return common;
}


*/

// Intersection of vertex lists of 2 graphs
double intersection_vertices(Graph &g1, Graph &g2)
{
        double common = 0;

        for(auto pr : global_vrtxlbl_map){
                unsigned vlbl = pr.second;
                common += min(g1.VertexLabelVectr[vlbl], g2.VertexLabelVectr[vlbl]);
        }
       //cout<<g1.gid <<"and "<<g2.gid<<"  In intersection_vertices common vertices = "<<common<<endl;

        return common;
}


double intersection_edges(Graph &g1, Graph &g2)
{

        double common = 0;

        /*for(auto label : e_label){

                common += min(g1.EdgeLabelMap[label], g2.EdgeLabelMap[label]);
        } */
        for(auto pr :global_edgetype_map){
                unsigned eType =pr.second;
                //cout<<pr.first<<" "<<pr.second<<endl;
                //cout<<g1.edgeTypeVectr[eType]<<"  "<<g2.edgeTypeVectr[eType]<<endl;
                common += min(g1.edgeTypeVectr[eType], g2.edgeTypeVectr[eType]);
        }
        //cout<<"In intersection_edges common edges = "<<common<<endl;

       return common;
}





// VEO Similarity computation
double VEO:: computeSimilarity(Graph &g1, Graph &g2, double &commonV)
{
/*	if(commonV == 0)
		commonV = intersection_vertices(g1.vertices, g2.vertices);
	commonV+=intersection_edges(g1.edges, g2.edges);
*/

	if(commonV == 0)
                commonV = intersection_vertices(g1,g2);
        commonV+=intersection_edges(g1,g2);


	double simScore = (double)(200.0*(commonV)/(double)(g1.vertexCount+g2.vertexCount+g1.edgeCount+g2.edgeCount));
	return simScore;
}

// Sorts vertex and edge set of graph dataset
void VEO:: sortVertexEdge(vector<Graph> &graph_dataset)
{
	for(int i = 0; i < graph_dataset.size(); i++)
	{
		sort(graph_dataset[i].vertices.begin(), graph_dataset[i].vertices.end());
		sort(graph_dataset[i].edges.begin(), graph_dataset[i].edges.end());
	}
}

void VEO:: printlist(int k)
{
	for(int i = 0; i < rankList[k].size(); i++)
		cout << i  << ": "<< rankList[k][i].first << endl;
}

//globaly ranks vertices and edges together
void VEO:: ranking(vector<Graph> &graph_dataset) 
{
	// traversing the graph-dataset
	for(int g_ind = 0; g_ind < graph_dataset.size(); g_ind++)
        {
                // traversing the vertex-set

                for(int vtx_ind = 0; vtx_ind < graph_dataset[g_ind].vrtxLbl.size(); vtx_ind++)
                {
                        pair<char, char> vtx_pair = make_pair(graph_dataset[g_ind].vertices[vtx_ind], '.');
                        if(fearture_freq.find(vtx_pair) != fearture_freq.end())
                                fearture_freq[vtx_pair]++;
                        else
                                fearture_freq[vtx_pair] = 0;    //To_CHECK
                }
                // traversing the edge-set
                for(int edge_ind = 0; edge_ind < graph_dataset[g_ind].edgeType.size(); edge_ind++)
                {
                        pair<char, char> edge_pair;
                        /*
                        if(graph_dataset[g_ind].edges[edge_ind].first > graph_dataset[g_ind].edges[edge_ind].second)
                                edge_pair = make_pair(graph_dataset[g_ind].edges[edge_ind].second, graph_dataset[g_ind].edges[edge_ind].first);
                        else
                                edge_pair = make_pair(graph_dataset[g_ind].edges[edge_ind].first, graph_dataset[g_ind].edges[edge_ind].second);
                        */
                        edge_pair = (graph_dataset[g_ind]).edgeType[edge_ind];
                        if(fearture_freq.find(edge_pair) != fearture_freq.end())
                                fearture_freq[edge_pair]++;
                        else
                                fearture_freq[edge_pair] = 0;
                }
        }

        vector<pair<pair<char, char>, unsigned long> > freqList;
        copy(fearture_freq.begin(), fearture_freq.end(), back_inserter(freqList));
        sort(freqList.begin(), freqList.end(), freqComp);

        // ranking each vertex and each edge
        r = 1;
        for(auto entry = freqList.begin(); entry!=freqList.end(); entry++)
        {
                fearture_freq[entry->first] = r;  // fearture_freq will have now rank instead of feequency of feature. as frequency of feature is now stored in freqList vector. so, no issue.
                r++;
        }

        // setting the size of invertecindex list as max value of rank - r
        InvertedIndex.resize(r);

}

void VEO:: buildPrefix(vector<Graph> &graph_dataset, int mode, bool isBucket, int no_of_buckets)
{
	
        rankList.resize(graph_dataset.size());
        //cout<<"rankList size in buildPrefix  = "<<rankList.size()<<endl;

	if(isBucket)
                bucket.resize(graph_dataset.size());

        for(int g_ind = 0; g_ind < graph_dataset.size(); g_ind++)
        {
                unsigned prefixLength = 0;

                double graph_size = graph_dataset[g_ind].vertexCount + graph_dataset[g_ind].edgeCount; // size of graph g_ind

                if(mode == 3)    // static Mode
                {
                        // For (g_ind)th graph prefix length will be this only for computation with any other graph.

                        double invUbound = (double)1.0/ubound; // inverse of ubound
                        prefixLength = 1 +(unsigned)(ceil((graph_size*(double)(1.0-invUbound)))); // Prefix Length
                }
                else            // Dynamic Mode
                {
                        // since it is dynamic, For (g_ind)th graph, prefix length will be calculated depending on other
                        //  graphs size. so for now we will compute rank list for full size.

                        prefixLength = graph_size;
                }

                // Constructing the rank-list 
                unordered_map<char, bool> vLabelOccur;           //To know whether I have already saved {rank, count} of that vertex feature or not.
                map< pair<char, char>, bool> eTypeOccur;     //To know whether I have already saved {rank, count} of that Edge feature or not. 
                // Putting vertex and edge ranks together
                vector<pair<unsigned long, unsigned long> > graph_ranks; //Will have Rank of the graph feature with it's no. of occurrence (Count) in that particular Graph.
                for(int vtx_ind = 0; vtx_ind < graph_dataset[g_ind].vrtxLbl.size(); vtx_ind++)
                {
                        char vlabl = graph_dataset[g_ind].vrtxLbl[vtx_ind];
                        if(vLabelOccur.find(vlabl) == vLabelOccur.end()) //If this vertex feature coming first time then save the corresponding {rank, count}. else ignore.
                        {
                                unsigned un_vlabl = global_vrtxlbl_map[vlabl];
                                pair<char, char> vtx_pair = make_pair(vlabl, '.');
				//Representing each graph in terms of Ranks. Ie.. relacing Vertices of the Graph using Rank. 
				//We are also saving count of that rank in this graph. 
                                graph_ranks.push_back({fearture_freq[vtx_pair], (graph_dataset[g_ind]).VertexLabelVectr[un_vlabl] });   //As, fearture_freq contains rank of feature now.
                                vLabelOccur[vlabl] =true;
                        }

                }

               for(int edge_ind = 0; edge_ind < graph_dataset[g_ind].edgeType.size(); edge_ind++)
                {
                        pair<char, char> edge_pair;
                        edge_pair = graph_dataset[g_ind].edgeType[edge_ind];
                        if(eTypeOccur.find(edge_pair) == eTypeOccur.end())    //If this Edge feature coming first time then save the corresponding {rank, count}. else ignore.
                        {
                                string etype {edge_pair.first, '-', edge_pair.second };
                                unsigned un_etype = global_edgetype_map[etype];
				//Representing each graph in terms of Ranks. Ie.. relacing Vertices of the Graph using Rank.
                                //We are also saving count of that rank in this graph.
                                graph_ranks.push_back({fearture_freq[edge_pair], (graph_dataset[g_ind]).edgeTypeVectr[un_etype]} );
                                eTypeOccur[edge_pair] =true;
                        }
                }

                // sort the ranks of the graph in descending order
                sort(graph_ranks.begin(), graph_ranks.end(), rankComp);


                // Traverse upto Whole part instead prefix_length. Bcoz we will apply suffix also.
                for(int pref = 0; pref < graph_size; pref++)
                        rankList[g_ind].push_back(graph_ranks[pref]);

                // rankList is till prefix length

                if(isBucket)
                {
                        bucket[g_ind].resize(no_of_buckets);
                        for(int buck_ind = 0; buck_ind < no_of_buckets; buck_ind++)
                                bucket[g_ind][buck_ind].resize(graph_size, 0);

                        vector<unsigned> sumBucket(no_of_buckets, 0);

                        // traversing graph's rank-list in ascending order
                        for(int grank_ind = graph_size-1; grank_ind >= 0; grank_ind--)
                        {
                                sumBucket[((graph_ranks[grank_ind]).first)%no_of_buckets]++;
                                for(int buck_ind = 0; buck_ind < no_of_buckets; buck_ind++)
                                {
                                        bucket[g_ind][buck_ind][grank_ind] = sumBucket[buck_ind];
                                }
                        }
                }
        }

}




////////////////////////////////////////////////////////////////////////////////////////
//Sparse table calculation


void VEO:: calculate_sparse_table(vector<Graph> &graph_dataset, int g_ind, long double minPrevSize)  // INVERTED INDEXING FOR RANKLISTS
{
	//cout<<"Reached to calculate_sparse_table "<<endl;
        // INVERTED INDEXING FOR RANKLISTS

                sparse_table.clear();   // re-setting the sparse table
                double graph_size = graph_dataset[g_ind].vertexCount + graph_dataset[g_ind].edgeCount;
		double invUbound = (double)1.0/ubound; // inverse of ubound
                int     prefixLength = 1 +(unsigned)(ceil((graph_size*(double)(1.0-invUbound)))); // Prefix Length
                //int prefixLength = rankList[g_ind].size();
                
		//cout<<" prefixLength   =  "<<prefixLength<<endl;
		//Finding exact index to go for this graph g_ind till what part of rankList comes in prefixLength. As it is Labeled Graph. so, a rank has more than one count.		
		int jthIndex =0, hop =0;
		//cout<<" rankList[g_ind].size() " <<rankList[g_ind].size()<<endl;
		for(int i =0; i< rankList[g_ind].size() ; i++)
		{
			hop = hop +(rankList[g_ind][i]).second;
			//cout<<" (rankList[g_ind][i]).second "<<(rankList[g_ind][i]).second<<endl;
			jthIndex = i;
			if(hop >= prefixLength) 
				break;
		}

		
		///////////////////////////////////////////////////
                // for this particular graph, make a sparse table with help of inverted list

		//cout<<"jthIndex = "<<jthIndex<<endl;
                for(int i = 0; i < jthIndex; i++){ // traversing a graph

                        unsigned long  rank = (rankList[g_ind][i]).first;
                        unsigned long rankCount_g_ind = (rankList[g_ind][i]).second;
			//cout<<" InvertedIndex[rank].size() = "<<InvertedIndex[rank].size()<<endl;
                        for(int j = InvertedIndex[rank].size()-1; j >=0; j--){ // traversing a rank's inverted list

                                int gr = InvertedIndex[rank][j].first.first;
                                long double PrevSize = graph_dataset[gr].vertexCount + graph_dataset[gr].edgeCount;
	                        if(PrevSize < minPrevSize)              // loose filter condition to break and stop unnecessary computation
        	                        break;

                	        if(gr == g_ind)         // skipping itself in inverted list
                        	        continue;


				unsigned long  rankCount_gr = InvertedIndex[rank][j].first.second;

                                //ss.insert ( InvertedIndex[rank][j] );
                                if( sparse_table.count( gr ) == 0)  // if occuring for the first time
                                        sparse_table[ gr ].first = min(rankCount_g_ind, rankCount_gr);
                                else
                                {
                                        sparse_table[ gr ].first += min(rankCount_g_ind, rankCount_gr);          // incrementing the count of commons in 2 graphs list upto a length
                                        sparse_table[ gr ].second = InvertedIndex[rank][j].second;  // storing that length of gr (shorter graph) upto which commons have been counted
                                }
                        }
                }
		
		                ///////////////////////////////////////////////////
                                //Now lets create InvertedIndex for (g_ind)th graph..

                for(int pref = 0; pref < jthIndex; pref++)
                {
                        unsigned long rank =(rankList[g_ind][pref]).first;
                        unsigned long rankCount_g_ind = (rankList[g_ind][pref]).second;
                        InvertedIndex[ rank ].push_back({make_pair(g_ind, rankCount_g_ind),  pref+1}); // pushing graph_id and position of the particular attribute/rank in that graph's list

                }


} //Sparse table function.









// index each input graphs in dataset
void VEO:: index(vector<Graph> &graph_dataset, int mode, bool isBucket, int no_of_buckets)
{	
	ranking(graph_dataset); // ranking vertices and edges together
	buildPrefix(graph_dataset, mode, isBucket, no_of_buckets);
}







// Applies prefix filter on Graph g1 and g2
bool VEO:: PrefixFilter(Graph &g1, Graph &g2, int index1, int index2, int mode, bool isBucket, int no_of_buckets, long unsigned &indexCount, double threshold)
{
        bool out = true;
        // atleast 1 common you got in the prefix part, then .
    if(sparse_table.count(index2) != 0)// or sparse_table[index2].count(index1) != 0)  // this sparse table is of g1/index1 graph
    {
                out=false;
    }

        if(out)
                return out;
        indexCount++;

        return out;
}




bool VEO:: PositioningFilter(Graph &g1, Graph &g2, int index1, int index2, int mode, bool isBucket, int no_of_buckets, long unsigned &indexCount, double threshold)
{


	return false;


}





/*


// Applies index filter on Graph g1 and g2
bool VEO:: indexFilter(Graph &g1, Graph &g2, int index1, int index2, int mode, bool isBucket, int no_of_buckets, long unsigned &indexCount, long unsigned &partitionCount, double threshold)
{
	unsigned size1 = g1.vertexCount + g1.edgeCount; // Size of Graph g1
	unsigned size2 = g2.vertexCount + g2.edgeCount; // Size of Graph g2

	unsigned prefix1;
	unsigned prefix2;

	if(mode == 3) // static Mode
	{
		double invUbound = (double)1.0/ubound; // inverse of ubound
		//prefix1 = rankList[index1].size(); // prefix-length of graph g1
		prefix1 = 1 +(unsigned)(ceil((size1*(double)(1.0-invUbound)))); // Prefix Length
		//prefix2 = rankList[index2].size(); // prefix-length of graph g2
		prefix2 = 1 +(unsigned)(ceil((size2*(double)(1.0-invUbound)))); // Prefix Length
		
	}
	if(mode == 4)	// Dynamic Mode
	{
		unsigned common = (unsigned)floor(((double)(threshold/200.0))*(double)(size1+size2));
		prefix1 = size1 - common + 1; // prefix-length of graph g1
		prefix2 = size2 - common + 1; // prefix-length of graph g2
	}

	unsigned start1 = 0;
	unsigned start2 = 0;
	bool out = true;
	long double commonTotal=0;

	// atleast 1 common you got in the prefix part, then .
    if(sparse_table[index1].count(index2) != 0 or sparse_table[index2].count(index1) != 0)
        //{out = false; cout<<"f\n";}
    {
		out=false;
        int partial_score, remaining1, remaining2;

        if(sparse_table[index1].count(index2) != 0)
            partial_score = sparse_table[index1][index2].first;
        else
            partial_score = sparse_table[index2][index1].first;

//cout<<partial_score<<"\n";
        remaining1 = size1 - prefix1;  // big graph
        remaining2 = size2 - sparse_table[index1][index2].second;   // sparse_table[index1][index2].second is the position upto which partial_score is stored in sparse table
        
        long double sizeSum = (long double)(size1 + size2);
        long double Common = (long double)1.0*(partial_score + min(remaining1, remaining2));
		commonTotal = partial_score;
        long double veoEstimate = (long double)(200.0*Common)/(sizeSum);

	//	if(g1.gid == 5 and g2.gid == 33)
       // cout<<index1<<" "<<index2<<" "<<" "<<veoEstimate<<"\n\n";
        if((long double)veoEstimate < (long double)threshold)
            out = true;
    }


	/*while(start1 < prefix1 && start2 < prefix2)
	{
		if(rankList[index1][start1] == rankList[index2][start2] && isBucket)
		{
			out = false;
			start1++;
			start2++;
			commonTotal++;		
		}
		else if(rankList[index1][start1] == rankList[index2][start2] && !isBucket)
		{
			out = false;
			break;
		}
		else if(rankList[index1][start1] > rankList[index2][start2])
			start1++;
		else
			start2++;
	}*/
	
	
/*	if(out)
		return out;
	indexCount++;
	//cout << "prefix1: " << prefix1 << " , prefix2: " << prefix2 << endl;
	//cout << "start1: " << start1 << " , start2: " << start2 << endl;
	start1 = prefix1;
	start2 = sparse_table[index1][index2].second;

	if(start1 != 0)
	 start1--;
	if(start2 != 0)
	 start2--;

	if(start1 < 0)
	{start1=0;cout<<"neg ";}
	if(start2 < 0)
	start2=0;
	
	//cout<<"done1 "<<start1<<" "<<start2<<"\n";
	if(isBucket)
	{
		for(int i = 0; i < no_of_buckets; i++)
		{
			//cout << bucket[index1][i][start1] << " , " << bucket[index2][i][start2] << endl;
			commonTotal += (long double)min(bucket[index1][i][start1], bucket[index2][i][start2]);
		}

		long double sizeSum = (long double)(g1.vertexCount + g1.edgeCount + g2.vertexCount + g2.edgeCount);
		long double veoEstimate =(long double)200.0*((long double)commonTotal/(long double)(sizeSum));
		if(veoEstimate < (long double)threshold)
			out = true;
		else
		{
			out = false;
			partitionCount++;
		}
	}
	return out;
}

*/


// MisMatching Filter
bool VEO:: mismatchingFilter(Graph &g1, Graph &g2, double &common, double threshold)
{
	int index1 = 0;
	int index2 = 0;
	unsigned degreeSum12 = 0;
	unsigned degreeSum21 = 0;
	while(index1 < g1.vertexCount && index2 < g2.vertexCount)
	{
		if(g1.vertices[index1] == g2.vertices[index2])
		{
			common++;
			index1++;
			index2++;
		}
		else if(g1.vertices[index1] > g2.vertices[index2])
		{
			degreeSum21 += g2.degrees[g2.vid_to_ind[g2.vertices[index2]]];
			index2++;
		}
		else
		{
			degreeSum12 += g1.degrees[g1.vid_to_ind[g1.vertices[index1]]];
			index1++;
		}
	}

	double ES = (double)min((g1.edgeCount - degreeSum12), (g2.edgeCount - degreeSum21) );
	ES = max(0.0,ES);
	double simEstimate = (double)(200.0*(double)((common + ES)/(double)(g1.edgeCount + g1.vertexCount + g2.edgeCount + g2.vertexCount)));
	return (simEstimate <= threshold);
}
