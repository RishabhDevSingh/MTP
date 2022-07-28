
#include <iostream>
#include<bits/stdc++.h>
using namespace std;


int main ()
{
    int n =5;
  int vertexcount;
 // cin >> vertexcount;
  int edgecount;
  //cin >> edgecount;
  vector < vector < int >>adj;
  /*adj[0].push_back (1);
  adj[0].push_back (2);
  adj[1].push_back (3);
  adj[1].push_back (4);
  adj[1].push_back (0);
  adj[2].push_back (3);
  adj[2].push_back (0);
  adj[2].push_back (4);
  adj[3].push_back (1);
  adj[3].push_back (2);
  adj[4].push_back (1);
  adj[4].push_back (2);
  adj[1].push_back (2);
  adj[2].push_back (1);
  int vc;
  int ec;*/
  vector<vector<int>> edges(5 , vector<int>(5,0));
 // cin>>vc>>ec;
  //pair< int ,int >p;
  //vector<int> edge;
  for( int k=0; k<5; ++k )
  {
    int i ,j;
    cin>>i>>j;
    edges[i][j]=1;
    //p  =make_pair(i,j);
    //edge.push_back(p);
    /*edge.push_back(i);
    edge.push_back(j);
    edges.push_back(edge);
    edge.clear();*/
  }
  cout<<edges[1][2]<<endl;
  for(auto e: edges)
  {
      for(int i=0;i<5;i++)
      {
          cout<<e[i];
      }
      cout<<endl;
  }
  
  int count =0;
  //cout<<edges[2][3]<<edges[0][1]<<edges[4][1];
  //map<int, int,int> mp;
  /*for (auto i: edges)
  {
      for(auto j:adj[i[0]])
      {
          for(auto k: adj[i[1]])
          {
              if( j == k )
              mp[i[0]] = k;
              mp[i[1]] = k;
          }
      }
      
  }
  for(  )
  {
      if( m == 1 )
      {
          count++;
      }
  }*/
  //cout<<count;
  
  for(int i =0;i<5;i++)
  for(int j=i+1;j<5;j++)
  {
      if(edges[i][j])
      {
          for(int k=j;k<n;k++)
          {
              if(edges[i][k] && edges[j][k])
              {
                  count++;
              }
          }
      }
  }
  cout<< count;


  return 0;
}


