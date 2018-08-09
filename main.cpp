#include "GRBInterface.h"
#include "KGraph.h"
#include "ConnectorEnumeration.h"
#include <sstream>
#include <string>
#include <ctime>
#include <vector>
#include <cmath>
#include <cstdlib>
#include <cstdio>

using namespace std;

int main(int argc, char *argv[])
{
	// change the precision for outputting the running time. I want it to be rounded to 2 decimal places.
	cout.setf(ios::fixed);
	cout.precision(2);

	time_t start_time = clock();

	if (argc<2)
		cerr << "ERROR: Not enough arguments.";
	else if (strcmp(argv[1], "Preprocess") == 0)  // only does the heuristic and preprocessing. Used to make table in paper.
	{
		// read the input graph g and the value for k
		KGraph g(argv[3], argv[3], argv[2]);
		long k = atol(argv[4]);
		cout << g.name << " " << g.n << " " << g.m << " " << k << " ";

		// start the timer and do heuristic & preprocessing
		time_t start = clock();
		vector<long> DROP_Solution = HeuristicAndPreprocess(g, k);

		// output # non-isolated vertices in preprocessed graph; # remaining edges; and preprocess+heuristic time:
		cout << DROP_Solution.size() << " " << g.ConnectedVertices() << " " << g.m << " " << (double)(clock() - start) / CLOCKS_PER_SEC;
		cout << "\n";
	}

	else if (strcmp(argv[1], "S2") == 0)
	{
		// read the input graph g and the value for k
		KGraph g(argv[3], argv[3], argv[2]);
		long k = 2;
		cout << g.name << " " << k << " " << g.n << " ";

		// start the timer and do heuristic & preprocessing
		time_t start = clock();
		vector<long> DROP_Solution = HeuristicAndPreprocess(g, k);


		bool subOpt;
		vector<long> Optimal_Solution = solve2Club(g, k, DROP_Solution, subOpt);

		// output solve info
		if (!subOpt) cout << Optimal_Solution.size() << " " << (double)(clock() - start) / CLOCKS_PER_SEC << endl;
		else cout << "? >3600" << endl;
		cout << "\n";
		if (subOpt) cerr << "Model is infeasible or other problem." << endl;
		else
		{
			cerr << "Maximum " << k << "club nodes are: ";
			for (long i = 0; i < Optimal_Solution.size(); i++)
				cerr << Optimal_Solution[i] << " ";
			cerr << endl;
		}
		cerr << "Is it actually a k-club? " << g.IsKClub(Optimal_Solution, k);
	}



	// Our main method
	else if (strcmp(argv[1], "CutLike") == 0)
	{
		// read the input graph g and the value for k
		KGraph g(argv[3], argv[3], argv[2]);
		long k = atol(argv[4]);
		cout << g.name << " " << k << " " << g.n << " ";

		// start the timer and do heuristic & preprocessing
		time_t start = clock();
		vector<long> DROP_Solution = HeuristicAndPreprocess(g, k);

		// solve max k-club using the cut-like formulation
		bool subOpt;
		vector<long> Optimal_Solution = solveMaxKClub_CutLike(g, k, DROP_Solution, subOpt);

		/*KGraph g1 = g.CreateInducedGraph(Optimal_Solution);
		for (long t = 0; t < g1.n; t++)
		{
		cerr << "it is : " << g1.degree[t] << endl;
		}*/

		// output solve info
		if (!subOpt) cout << Optimal_Solution.size() << " " << (double)(clock() - start) / CLOCKS_PER_SEC << endl;
		else cout << "? >3600" << endl;
		cout << "\n";
		if (subOpt) cerr << "Model is infeasible or other problem." << endl;
		else
		{
			cerr << "# nodes in " << k << "club = " << Optimal_Solution.size() << endl;
			cerr << "Maximum " << k << "club nodes are: ";
			for (long i = 0; i < Optimal_Solution.size(); i++)
				cerr << Optimal_Solution[i] << " ";
			cerr << endl;
		}
		cerr << "Is it actually a k-club? " << g.IsKClub(Optimal_Solution, k);
	}

	//3clubpathlike
	else if (strcmp(argv[1], "3ClubPathLike") == 0)
	{
		// read the input graph g and the value for k
		KGraph g(argv[3], argv[3], argv[2]);
		long k = 3;
		cout << g.name << " " << g.n << " ";

		// start the timer and do heuristic & preprocessing
		time_t start = clock();
		vector<long> DROP_Solution = HeuristicAndPreprocess(g, k);

		// solve max k-club using the path-like formulation
		bool subOpt;
		vector<long> Optimal_Solution = Pathlike_3Club(g, DROP_Solution, subOpt);

		// output solve info
		if (!subOpt) cout << Optimal_Solution.size() << " " << (double)(clock() - start) / CLOCKS_PER_SEC << endl;
		else cout << "? >3600" << endl;
		cout << "\n";
		if (subOpt) cerr << "Model is infeasible or other problem." << endl;
		else
		{
			cerr << "Maximum " << k << "club nodes are: ";
			for (long i = 0; i < Optimal_Solution.size(); i++)
				cerr << Optimal_Solution[i] << " ";
			cerr << endl;
		}
		cerr << "Is it actually a k-club? " << g.IsKClub(Optimal_Solution, k);
	}

	// 4clubpathlike
	else if (strcmp(argv[1], "4ClubPathLike") == 0)
	{
		// read the input graph g and the value for k
		KGraph g(argv[3], argv[3], argv[2]);
		long k = 4;
		cout << g.name << " " << g.n << " ";

		// start the timer and do heuristic & preprocessing
		time_t start = clock();
		vector<long> DROP_Solution = HeuristicAndPreprocess(g, k);

		// solve max k-club using the path-like formulation
		bool subOpt;
		vector<long> Optimal_Solution = Pathlike_4Club(g, DROP_Solution, subOpt);

		// output solve info
		if (!subOpt) cout << Optimal_Solution.size() << " " << (double)(clock() - start) / CLOCKS_PER_SEC << endl;
		else cout << "? >3600" << endl;
		cout << "\n";
		if (subOpt) cerr << "Model is infeasible or other problem." << endl;
		else
		{
			cerr << "Maximum " << k << "club nodes are: ";
			for (long i = 0; i < Optimal_Solution.size(); i++)
				cerr << Optimal_Solution[i] << " ";
			cerr << endl;
		}
		cerr << "Is it actually a k-club? " << g.IsKClub(Optimal_Solution, k);
	}


	// Solving Maximum k-club using Veremyev et al. formulation.
	else if (strcmp(argv[1], "Veremyev") == 0)
	{
		// read the input graph g and the value for k
		KGraph g(argv[3], argv[3], argv[2]);
		long k = atol(argv[4]);
		cout << g.name << " " << k << " " << g.n << " ";

		// start the timer and do heuristic & preprocessing
		time_t start = clock();
		vector<long> DROP_Solution = HeuristicAndPreprocess(g, k);

		// solve max k-club using the Veremyev et al. formulation
		bool subOpt;
		vector<long> Optimal_Solution = MaxKclubRevisedVeremyevFormulation(g, k, DROP_Solution, subOpt);

		// output solve info
		if (!subOpt) cout << Optimal_Solution.size() << " " << (double)(clock() - start) / CLOCKS_PER_SEC << endl;
		else cout << "? >3600" << endl;
		cout << "\n";
		if (subOpt) cerr << "Model is infeasible or other problem." << endl;
		else
		{
			cerr << "Maximum " << k << "club nodes are: ";
			for (long i = 0; i < Optimal_Solution.size(); i++)
				cerr << Optimal_Solution[i] << " ";
			cerr << endl;
		}
		cerr << "Is it actually a k-club? " << g.IsKClub(Optimal_Solution, k);
	}

	// Solve max k-club, using CHC cuts in callback
	else if (strcmp(argv[1], "CHC") == 0)
	{
		// read the input graph g and the value for k
		KGraph g(argv[3], argv[3], argv[2]);
		long k = atol(argv[4]);
		cout << g.name << " " << k << " " << g.n << " ";

		// start the timer and do heuristic & preprocessing
		time_t start = clock();
		vector<long> DROP_Solution = HeuristicAndPreprocess(g, k);


		bool subOpt;
		vector<long> Optimal_Solution = solveMaxClubCHC(g, k, DROP_Solution, subOpt);

		// output solve info
		if (!subOpt) cout << Optimal_Solution.size() << " " << (double)(clock() - start) / CLOCKS_PER_SEC << endl;
		else cout << "? >3600" << endl;
		cout << "\n";
		if (subOpt) cerr << "Model is infeasible or other problem." << endl;
		else
		{
			cerr << "Maximum " << k << "club nodes are: ";
			for (long i = 0; i < Optimal_Solution.size(); i++)
				cerr << Optimal_Solution[i] << " ";
			cerr << endl;
		}
		cerr << "Is it actually a k-club? " << g.IsKClub(Optimal_Solution, k);
	}

	else
	{
		cout << "ERROR: Your command is not valid.";
	}
	return EXIT_SUCCESS;
}
