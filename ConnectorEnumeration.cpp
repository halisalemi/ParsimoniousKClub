#include "ConnectorEnumeration.h"
#include <map>



using namespace std;


vector<long> sortnodes(long i, long j, long k)
{
	vector<long> sorted;
	sorted.push_back(i);
	sorted.push_back(j);
	sorted.push_back(k);

	sort(sorted.begin(), sorted.end());

	return sorted;
}

map<vector<long>, long, classcomp> EnumerateLength3Connector(KGraph &g)
{
	//create map, mapsize, and path
	map<vector<long>, long, classcomp>map;
	long mapsize = 0;
	vector<long> path;

	for (long a = 0; a < g.n; a++)
	{
		//compute SP from node a 
		vector<long> dist_from_a = g.ShortestPathsUnweighted(a);
		for (long b = a + 1; b < g.n; b++)
		{
			// If dist(a,b)=1 then nothing is needed to connect them. 
			// If dist(a,b)>3 then they cannot be connected with short path.
			if (dist_from_a[b] == 2 || dist_from_a[b] == 3)
			{
				//compute SP from node b 
				vector<long> dist_from_b = g.ShortestPathsUnweighted(b);

				for (long a_neighbors_iterator = 0; a_neighbors_iterator < g.adj[a].size(); a_neighbors_iterator++)
				{
					long i = g.adj[a][a_neighbors_iterator];
					if (dist_from_b[i] == 2)  // if i belongs to V_{12}
					{
						for (long i_neighbors_iterator = 0; i_neighbors_iterator < g.adj[i].size(); i_neighbors_iterator++)
						{
							long j = g.adj[i][i_neighbors_iterator];
							if (dist_from_a[j] == 2 && dist_from_b[j] == 1) // if j\in N(i) belongs to V_{21}
							{
								//Now, we have found the a,i,j,b path

								//checking whether a,i,j,b path does not belong to map and if not add it to the map		
								long minij = min(i, j);
								long maxij = max(i, j);
								std::map<std::vector<long>, long>::iterator it = map.find({ minij,maxij });
								if (it == map.end())
								{
									path = { minij, maxij };
									//cerr << minij << " " << maxij << endl;
									map.insert(pair<vector<long>, long>(path, mapsize));
									mapsize++;
								}
							}
						}
					}
				}
			}
		}
	}
	return map;
}

std::map<std::vector<long>, long, classcomp> EnumerateLength4Connector(KGraph &g)
{
	map<vector<long>, long, classcomp>map;
	long mapsize = 0;
	EnumerateLength4ConnectorType1(g, map, mapsize);
	EnumerateLength4ConnectorType2(g, map, mapsize);
	EnumerateLength4ConnectorType3(g, map, mapsize);
	EnumerateLength4ConnectorType4(g, map, mapsize);

	return map;
}

void EnumerateLength4ConnectorType1(KGraph &g, map<vector<long>, long, classcomp>&map, long &mapsize)
{

	//create map, mapsize, and path
	//map<vector<long>, long, classcomp>map;
	//long mapsize = 0;
	vector<long> path;

	for (long a = 0; a < g.n; a++)
	{
		//compute SP from node a 
		vector<long> dist_from_a = g.ShortestPathsUnweighted(a);
		for (long b = a + 1; b < g.n; b++)
		{
			// If dist(a,b)=1 then nothing is needed to connect them. 
			// If dist(a,b)>4 then they cannot be connected with short path.
			if (dist_from_a[b] == 2 || dist_from_a[b] == 3 || dist_from_a[b] == 4)
			{
				//compute SP from node b 
				vector<long> dist_from_b = g.ShortestPathsUnweighted(b);
				for (long a_neighbors_iterator = 0; a_neighbors_iterator < g.adj[a].size(); a_neighbors_iterator++)
				{
					long i = g.adj[a][a_neighbors_iterator];
					if (dist_from_b[i] == 3) // if i belongs to V_{13}
					{
						for (long i_neighbors_iterator = 0; i_neighbors_iterator < g.adj[i].size(); i_neighbors_iterator++)
						{
							long j = g.adj[i][i_neighbors_iterator];
							if (dist_from_a[j] == 2 && dist_from_b[j] == 2)  // if j belongs to V_{22}
							{
								for (long j_neighbors_iterator = 0; j_neighbors_iterator < g.adj[j].size(); j_neighbors_iterator++)
								{
									long k = g.adj[j][j_neighbors_iterator];
									if (dist_from_a[k] == 2 && dist_from_b[k] == 1)  // if k belongs to V_{21}
									{
										//Now, we have found the a,i,j,,k,b path

										//checking whether a,i,j,b path does not belong to map and if not add it to the map	

										vector<long> sortedsubset = sortnodes(i, j, k);
										std::map<std::vector<long>, long>::iterator it = map.find(sortedsubset);
										if (it == map.end())
										{
											map.insert(pair<vector<long>, long>(sortedsubset, mapsize));
											//cerr << "When a = " << a << " b = " << b << ", we have connectors: " << sortedsubset[0] << " " << sortedsubset[1] << " " << sortedsubset[2] << "\n";
											mapsize++;
										}
									}
								}
							}
						}
					}
				}
			}
		}
	}
	//return map;
}

void EnumerateLength4ConnectorType2(KGraph &g, map<vector<long>, long, classcomp>&map, long &mapsize)
{
	//create map, mapsize, and path
	//map<vector<long>, long, classcomp>map;
	//long mapsize = 0;
	vector<long> path;

	for (long a = 0; a < g.n; a++)
	{
		//compute SP from node a 
		vector<long> dist_from_a = g.ShortestPathsUnweighted(a);
		for (long b = a + 1; b < g.n; b++)
		{
			// If dist(a,b)=1 then nothing is needed to connect them. 
			// If dist(a,b)>4 then they cannot be connected with short path.
			if (dist_from_a[b] == 2 || dist_from_a[b] == 3 || dist_from_a[b] == 4)
			{
				//compute SP from node b 
				vector<long> dist_from_b = g.ShortestPathsUnweighted(b);
				for (long a_neighbors_iterator = 0; a_neighbors_iterator < g.adj[a].size(); a_neighbors_iterator++)
				{
					long i = g.adj[a][a_neighbors_iterator];
					if (dist_from_b[i] == 2)  //if i belongs to V_{12}
					{
						for (long i_neighbors_iterator = 0; i_neighbors_iterator < g.adj[i].size(); i_neighbors_iterator++)
						{
							long j = g.adj[i][i_neighbors_iterator];
							if (dist_from_a[j] == 2 && dist_from_b[j] == 2) //if j belongs to V_{22}
							{
								for (long j_neighbors_iterator = 0; j_neighbors_iterator < g.adj[j].size(); j_neighbors_iterator++)
								{
									long k = g.adj[j][j_neighbors_iterator];
									if (dist_from_a[k] == 3 && dist_from_b[k] == 1) //if k belongs to V_{31}
									{
										//Now, we have found the a,i,j,,k,b path

										//checking whether a,i,j,b path does not belong to map and if not add it to the map	

										vector<long> sortedsubset = sortnodes(i, j, k);
										std::map<std::vector<long>, long>::iterator it = map.find(sortedsubset);
										if (it == map.end())
										{
											map.insert(pair<vector<long>, long>(sortedsubset, mapsize));
											//cerr << "When a = " << a << " b = " << b << ", we have connectors: " << sortedsubset[0] << " " << sortedsubset[1] << " " << sortedsubset[2] << "\n";
											mapsize++;
										}
									}
								}
							}
						}
					}
				}
			}
		}
	}
	//return map; 
}

void EnumerateLength4ConnectorType3(KGraph &g, map<vector<long>, long, classcomp>&map, long &mapsize)
{
	//create map, mapsize, and path
	//map<vector<long>, long, classcomp>map;
	//long mapsize = 0;
	vector<long> path;

	for (long a = 0; a < g.n; a++)
	{
		//compute SP from node a 
		vector<long> dist_from_a = g.ShortestPathsUnweighted(a);
		for (long b = a + 1; b < g.n; b++)
		{
			// If dist(a,b)=1 then nothing is needed to connect them. 
			// If dist(a,b)>4 then they cannot be connected with short path.
			if (dist_from_a[b] == 2 || dist_from_a[b] == 3 || dist_from_a[b] == 4)
			{
				//compute SP from node b 
				vector<long> dist_from_b = g.ShortestPathsUnweighted(b);
				for (long a_neighbors_iterator = 0; a_neighbors_iterator < g.adj[a].size(); a_neighbors_iterator++)
				{
					long i = g.adj[a][a_neighbors_iterator];
					if (dist_from_b[i] == 3) //if i belongs to V_{13}
					{
						for (long i_neighbors_iterator = 0; i_neighbors_iterator < g.adj[i].size(); i_neighbors_iterator++)
						{
							long j = g.adj[i][i_neighbors_iterator];
							if (dist_from_a[j] == 2 && dist_from_b[j] == 2) //if j belongs to V_{22}
							{
								for (long j_neighbors_iterator = 0; j_neighbors_iterator < g.adj[j].size(); j_neighbors_iterator++)
								{
									long k = g.adj[j][j_neighbors_iterator];
									if (dist_from_a[k] == 3 && dist_from_b[k] == 1) //if k belongs to V_{31}
									{
										//Now, we have found the a,i,j,,k,b path

										//checking whether a,i,j,b path does not belong to map and if not add it to the map	

										vector<long> sortedsubset = sortnodes(i, j, k);
										std::map<std::vector<long>, long>::iterator it = map.find(sortedsubset);
										if (it == map.end())
										{
											map.insert(pair<vector<long>, long>(sortedsubset, mapsize));
											//cerr << "When a = " << a << " b = " << b << ", we have connectors: " << sortedsubset[0] << " " << sortedsubset[1] << " " << sortedsubset[2] << "\n";
											mapsize++;
										}
									}
								}
							}
						}
					}
				}
			}
		}
	}
	//return map;
}


void EnumerateLength4ConnectorType4(KGraph &g, map<vector<long>, long, classcomp>&map, long &mapsize)
{
	//create map, mapsize, and path
	//map<vector<long>, long, classcomp>map;
	//long mapsize = 0;
	vector<long> path;

	for (long a = 0; a < g.n; a++)
	{
		//compute SP from node a 
		vector<long> dist_from_a = g.ShortestPathsUnweighted(a);
		for (long b = a + 1; b < g.n; b++)
		{
			// If dist(a,b)=1 then nothing is needed to connect them. 
			// If dist(a,b)>4 then they cannot be connected with short path.
			if (dist_from_a[b] == 2 || dist_from_a[b] == 3 || dist_from_a[b] == 4)
			{
				//compute SP from node b 
				vector<long> dist_from_b = g.ShortestPathsUnweighted(b);
				for (long a_neighbors_iterator = 0; a_neighbors_iterator < g.adj[a].size(); a_neighbors_iterator++)
				{
					long i = g.adj[a][a_neighbors_iterator];

					if (dist_from_b[i] == 2) //if i belongs to V_{12}
					{
						// boolify neighborhood of i
						vector<bool> N_i(g.n, false);
						for (long NeighborNumber = 0; NeighborNumber < g.degree[i]; NeighborNumber++)
						{
							N_i[g.adj[i][NeighborNumber]] = true;
						}

						// look for next node j in our a-i-j-k-b path
						for (long i_neighbors_iterator = 0; i_neighbors_iterator < g.adj[i].size(); i_neighbors_iterator++)
						{
							long j = g.adj[i][i_neighbors_iterator];
							if (dist_from_a[j] == 2 && dist_from_b[j] == 2) //if j belongs to V_{22}
							{
								for (long j_neighbors_iterator = 0; j_neighbors_iterator < g.adj[j].size(); j_neighbors_iterator++)
								{
									long k = g.adj[j][j_neighbors_iterator];
									if (N_i[k] == false) // k cannot be neighbor of i 
									{
										if (dist_from_a[k] == 2 && dist_from_b[k] == 1) //if k belongs to V_{21}
										{
											//Now, we have found the a,i,j,,k,b path

											//checking whether a,i,j,b path does not belong to map and if not add it to the map	

											vector<long> sortedsubset = sortnodes(i, j, k);
											std::map<std::vector<long>, long>::iterator it = map.find(sortedsubset);
											if (it == map.end())
											{
												map.insert(pair<vector<long>, long>(sortedsubset, mapsize));
												//cerr << "When a = " << a << " b = " << b << ", we have connectors: " << sortedsubset[0] << " " << sortedsubset[1] << " " << sortedsubset[2] << "\n";
												mapsize++;
											}
										}
									}
								}
							}
						}
					}
				}
			}
		}
	}
	//return map;
}





