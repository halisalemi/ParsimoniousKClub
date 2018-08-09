#include <iostream>
#include <map>
#include <vector>
#include <algorithm>
#include "KGraph.h"


// Use to impose a total ordering on paths
struct classcomp
{
	// if lhs comes before rhs return true
	bool operator()(const std::vector<long>& lhs, const std::vector<long>& rhs) const
	{
		//smaller path comes first. As an example {7,9} comes befor {4,6,7} 
		if (lhs.size() < rhs.size()) return 1;
		if (lhs.size() > rhs.size()) return 0;

		// if paths have same length lexicographically smaller ones comes first
		for (int i = 0; i < lhs.size(); i++)
		{
			if (lhs[i] < rhs[i]) return 1;
			if (lhs[i] > rhs[i]) return 0;
		}

		return 0;
	}
};

std::map<std::vector<long>, long, classcomp> EnumerateLength3Connector(KGraph &g);

std::map<std::vector<long>, long, classcomp> EnumerateLength4Connector(KGraph &g);

void EnumerateLength4ConnectorType1(KGraph &g, map<vector<long>, long, classcomp>&map, long &mapsize);

void EnumerateLength4ConnectorType2(KGraph &g, map<vector<long>, long, classcomp>&map, long &mapsize);

void EnumerateLength4ConnectorType3(KGraph &g, map<vector<long>, long, classcomp>&map, long &mapsize);

void EnumerateLength4ConnectorType4(KGraph &g, map<vector<long>, long, classcomp>&map, long &mapsize);



vector<long> sortnodes(long i, long j, long k);

