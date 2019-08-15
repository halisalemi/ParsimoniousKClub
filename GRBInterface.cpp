#include "GRBInterface.h"
#include "ConnectorEnumeration.h"
#include <sstream>
#include <ctime>
#include <iostream>
#include <algorithm>
#include <queue>
#include <limits.h>
#include <string.h>

using namespace std;

string itos(int i) { stringstream s; s << i; return s.str(); }

vector<bool> boolify(vector<long> &S, long n)
{
	vector<bool> Sbool(n, false);
	for (long i = 0; i < S.size(); i++) Sbool[S[i]] = true;
	return Sbool;
}

vector<long> HeuristicAndPreprocess(KGraph &g, long k)
{
	// find k-th power of graph g and call it gk
	KGraph gk = g.CreatePowerGraph(k);

	vector<long> rd;  // right-degree of degeneracy ordering.
	vector<long> degeneracyordering = gk.FindDegeneracyOrdering(rd);
	vector<long> kclique = gk.FindHeuristicClique(degeneracyordering, rd);
	//cerr << "k-clique (heuristic) size is = " << kclique.size() << endl;

	// perform DROP heuristic, with kclique as starting point.
	vector<bool> HeuristicSolution = boolify(kclique, g.n); // make the heuristic kclique into a boolean form
	vector<long> DROP_Solution = DROP(g, HeuristicSolution, k);
	long lb = DROP_Solution.size();
	//cerr << "Drop heuristic gives LB = " << lb << endl;

	// perform preprocessing
	vector<long> kcorevertices = gk.FindVerticesOfKCore(degeneracyordering, rd, lb - 1);
	//cerr << "Preprocessed instances has this many vertices = " << kcorevertices.size() << endl;
	g.FindInducedGraph(kcorevertices);

	return DROP_Solution;
}

vector<long> DROP(KGraph &g, long k)
{
	vector<bool> W(g.n, true);	// our k-club (eventually)
	return DROP(g, W, k);
}
vector<long> DROP(KGraph &g, vector<bool> &W, long k)
{
	vector<long> near(g.n, 0);	// the number of vertices "nearby" to a vertex (within k-hops)

	long Wsize = count(W.begin(), W.end(), true);
	// while W is not an k-club, delete some "bad" vertex. The number of vertices in W is size.
	for (long size = Wsize; size >= 0; size--)
	{
		// compute how many vertices are "nearby" to each vertex.
		for (long i = 0; i < g.n; i++)
		{
			if (!W[i]) continue;
			near[i] = KNeighborhoodSize(g, W, k, i);
		}

		// pick vertex w (in W) with smallest "nearby" vertices
		long w;	// vertex with smallest number of nearby vertices
		long smallestNearby = size; // number of vertices nearby to w.
		for (long i = 0; i < g.n; i++)
		{
			if (!W[i]) continue;
			if (near[i] < smallestNearby)
			{
				w = i;
				smallestNearby = near[i];
			}
		}

		// check if k-club. If not, then remove vertex w.
		if (smallestNearby == size) break;
		W[w] = false;
	}

	// convert W to vector<long> form, and return
	vector<long> Wvector;
	for (long i = 0; i < g.n; i++)
	{
		if (W[i])
		{
			Wvector.push_back(i);
		}
	}
	return Wvector;
}

long Kclub_callback::numCallbacks = 0;
double Kclub_callback::TotalCallbackTime = 0;
long Kclub_callback::numLazyCuts = 0;

long CHC::numCallbacks = 0;
double CHC::TotalCallbackTime = 0;
long CHC::numLazyCuts = 0;


// calback function for our main method - using MINIMALIZE 
// to obtain minimal subset of length-k a,b-separator
void Kclub_callback::callback()
{
	try {
		if (where == GRB_CB_MIPSOL)
		{
			numCallbacks++;
			time_t start = clock();

			// get the ''solution (S)'' from Gurobi
			double *x = new double[g1.n];
			x = getSolution(vars, g1.n);

			// make it boolean and call it D
			vector<bool> D(g1.n, false);

			for (long i = 0; i < g1.n; i++)
			{
				if (x[i] > 0.5) D[i] = true;
			}

			if (count(D.begin(), D.end(), false) != g1.n)
			{
				vector<long> SelectedVertices;
				vector<long> C_Prime;

				for (long i = 0; i < g1.n; i++)
				{
					if (D[i]) SelectedVertices.push_back(i);
					else C_Prime.push_back(i);
				}

				// create G[D], which we call g2
				vector<long> Rmap;
				KGraph g2 = g1.CreateInducedGraph(D, Rmap);
				for (long i = 0; i < g2.n; i++)
				{
					vector<long> dist_from_i_to = g2.ShortestPathsUnweighted(i);
					for (long j = i + 1; j < g2.n; j++)
					{
						if (dist_from_i_to[j] > k1)
						{
							long a = SelectedVertices[i];
							long b = SelectedVertices[j];

							// C_Prime is a length-s a,b-separator. now make it minimal
							vector<long> minimal_length_k_ab_separator = MINIMALIZE(g1, a, b, C_Prime, k1);
							GRBLinExpr expr = vars[a] + vars[b];
							for (long q = 0; q < minimal_length_k_ab_separator.size(); q++)
							{
								long v = minimal_length_k_ab_separator[q];
								expr -= vars[v];
							}
							////additions for lifted inequality x(D)-x(S)<=1
							//vector<long> D;
							//D.push_back(a);
							//D.push_back(b);
							//vector<long> dist_from_D = g1.MultiSourceShortestPaths(D);
							//for (long q = 0; q < g1.n; q++)
							//{
								//if (dist_from_D[q] <= k1) continue;
								//D.push_back(q);
								//dist_from_D = g1.MultiSourceShortestPaths(D);
								//expr += vars[q];
							//}

							////// new additions for lifting
							//vector<bool> complementOfC(g1.n, true);
							//for (long w = 0; w < minimal_length_k_ab_separator.size(); w++)
							//	complementOfC[minimal_length_k_ab_separator[w]] = false;
							//vector<long> dist_from_a = g1.ShortestPathsUnweighted(a, complementOfC);
							//vector<long> dist_from_b = g1.ShortestPathsUnweighted(b, complementOfC);
							//vector<bool> CanLift(g1.n, true);

							//// check if belongs to same component in G
							//vector<long> dist_in_G = g1.ShortestPathsUnweighted(a);
							//for (long d = 0; d < g1.n; d++) 
							//	if (dist_in_G[d] == g1.n) 
							//		CanLift[d] = false;

							//// check if satisifes property 1:
							//for (long d = 0; d < g1.n; d++)
							//	if (dist_from_a[d] <= k1 || dist_from_b[d] <= k1) 
							//		CanLift[d] = false;

							//// now check property 2:
							//for (long c = 0; c < minimal_length_k_ab_separator.size(); c++)
							//{
							//	complementOfC[minimal_length_k_ab_separator[c]] = true;

							//	vector<long> dist_from_c = g1.ShortestPathsUnweighted(c, complementOfC);

							//	// check if star-like k-club containing a,b,d with c as a sort of hub.
							//	for (long d = 0; d < g1.n; d++)
							//		if (dist_from_c[a] + dist_from_c[d] <= k1 
							//			&& dist_from_c[b] + dist_from_c[d] <= k1)
							//			CanLift[d] = false;

							//	complementOfC[minimal_length_k_ab_separator[c]] = false;
							//}

							//bool lifted = false;
							//for (long d = 0; d < g1.n && !lifted; d++)
							//{
							//	if (CanLift[d])
							//	{
							//		expr += vars[d];
							//		lifted = true;
							//		cerr << ".";
							//	}
							//}
							////// end new things for lifting
							addLazy(expr <= 1);
							numLazyCuts++;
						}
					}
				}
			}
			TotalCallbackTime += (double)(clock() - start) / CLOCKS_PER_SEC;
			delete[] x;
		}
	}
	catch (GRBException e) {
		cout << "Error number: " << e.getErrorCode() << endl;
		cout << e.getMessage() << endl;
	}
	catch (...) {
		cout << "Error during callback" << endl;
	}
}


vector<long> MINIMALIZE(KGraph &g, long a, long b, vector<long> ab_Separator, long k)
{
	vector<long> MinimalCut; // what we return at end of function
	vector<bool> Cut(g.n, false); // a boolean representation of the cut. 

								  // first do some linear-time preprocessing to remove vertices that are not on length-k a,b-path
								  //		for example, this removes vertices that belong to a different component in G.
	vector<long> dist_from_a = g.ShortestPathsUnweighted(a);
	vector<long> dist_from_b = g.ShortestPathsUnweighted(b);
	for (long i = 0; i < ab_Separator.size(); i++)
	{
		long v = ab_Separator[i];
		if (dist_from_a[v] + dist_from_b[v] <= k)
		{
			Cut[v] = true;  // vertex v was in ab_Separator AND belongs to a length-k a,b-path in graph G
		}
	}

	// initialize VerticesInGc = V \ Cut
	vector<bool> VerticesInGc(g.n, true);
	for (long i = 0; i < g.n; i++) if (Cut[i]) VerticesInGc[i] = false;

	// now run the normal minimalize algorithm.
	for (long c = 0; c<g.n; c++)
	{
		if (!Cut[c]) continue; // only test for cut vertices
		if (dist_from_a[c] == 1 && dist_from_b[c] == 1) continue; // if c \in N(a)\cap N(b), then c belongs to every minimal cut.

																  // put vertex c in G_c
		VerticesInGc[c] = true;

		// compute distances from c in G_c
		vector<long> dist_from_c = g.ShortestPathsUnweighted(c, VerticesInGc);

		// check if vertex c would close the distance from a to b to be at most k.
		if (dist_from_c[a] + dist_from_c[b] <= k) // vertex c remains in cut
		{
			VerticesInGc[c] = false;
		}
		else // vertex c is pulled from the cut.
		{
			Cut[c] = false;
		}
	}
	for (long i = 0; i < g.n; i++) if (Cut[i]) MinimalCut.push_back(i);
	return MinimalCut;
}

// using canonical hypercube cuts (CHC)
void CHC::callback()
{
	try {
		if (where == GRB_CB_MIPSOL)
		{
			//cerr << "callback";
			numCallbacks++;
			time_t start = clock();

			// get the ''solution - Vector S'' from Gurobi
			double *x = new double[g1.n];
			x = getSolution(vars, g1.n);

			// boolify S and call it D
			vector<bool> D(g1.n, false);
			for (long i = 0; i < g1.n; i++)	if (x[i] > 0.5) D[i] = true;

			if (count(D.begin(), D.end(), false) != g1.n)
			{
				if (g1.IsKClub(D, k1) == false)
				{
					GRBLinExpr expr;
					for (long i = 0; i < g1.n; i++)
					{
						if (D[i]) expr += (1 - vars[i]);
						else expr += vars[i];
					}
					addLazy(expr >= 1);
				}
			}

			TotalCallbackTime += (double)(clock() - start) / CLOCKS_PER_SEC;
			delete[] x;
		}
	}
	catch (GRBException e) {
		cout << "Error number: " << e.getErrorCode() << endl;
		cout << e.getMessage() << endl;
	}
	catch (...) {
		cout << "Error during callback" << endl;
	}
}



long KNeighborhoodSize(KGraph &g, vector<bool> &W, long k, long v)
{
	if (!W[v])
	{
		cerr << "\n ERROR: K-neighborhoodSize. You are calculating distances across W nodes starting from some vertex v which does not belong to W.";
		return 0;
	}
	vector<long> dist = g.ShortestPathsUnweighted(v, W);
	long nsize = 0;		// the number of nodes j with dist(v,j) <= s (including node v).
	for (long i = 0; i < g.n; i++) if (dist[i] <= k) nsize++;
	return nsize;
}
long KNeighborhoodSize(KGraph &g, long k, long v)
{
	vector<bool> W(g.n, true);
	return KNeighborhoodSize(g, W, k, v);
}

vector<long> solve2Club(KGraph &g, long k, vector<long> &HeuristicSolution, bool &subOpt)
{

	vector<long> S;
	subOpt = true;
	vector< vector< long> > components;
	vector<long> degreeZero;
	g.FindConnectedComponents(components, degreeZero);
	long lb = HeuristicSolution.size();

	try
	{
		GRBEnv env = GRBEnv();
		env.set(GRB_IntParam_OutputFlag, 0);
		env.set(GRB_IntParam_Method, 3); //concurrent  
		env.set(GRB_DoubleParam_TimeLimit, 3600);
		env.set(GRB_DoubleParam_MIPGap, 0.0);
		GRBModel model = GRBModel(env);
		GRBVar *X = model.addVars(g.n, GRB_BINARY);
		GRBVar *Y = model.addVars(components.size(), GRB_BINARY);
		model.update();


		// Adding objective function
		for (long w = 0; w < g.n; w++)
			X[w].set(GRB_DoubleAttr_Obj, 1);
		model.set(GRB_IntAttr_ModelSense, GRB_MAXIMIZE);


		// Fixing singletons to zero.
		for (long i = 0; i < degreeZero.size(); i++)
		{
			long v = degreeZero[i];
			model.addConstr(X[v] == 0);
		}

		// Fixing y[i]=0 when |V(G_i)| < lb.
		for (long i = 0; i < components.size(); i++)
		{
			if (components[i].size() < lb)
			{
				model.addConstr(Y[i] == 0);
			}
		}

		// Adding X[v] + X[q] - X[N(a) \cap N(b)] <= 1 
		for (long i = 0; i < components.size(); i++)
		{
			// add all constraints within component i
			for (long j = 0; j < components[i].size(); j++)
			{
				long v = components[i][j]; // v is a vertex in component i.
				vector<bool> neighbors(g.n, false);
				long w;
				for (long q = 0; q < g.degree[v]; q++)
				{
					w = g.adj[v][q];
					neighbors[w] = true;
				}

				for (long q = v + 1; q < g.n; q++)
				{
					if (neighbors[q]) continue;
					vector<long> p = g.CommonNeighborsList(v, q);

					GRBLinExpr expr = 0;
					for (long k = 0; k < p.size(); k++)
					{
						expr += X[p[k]];
					}
					expr = X[v] + X[q] - expr;
					model.addConstr(expr <= 1);
				}
			}
		}

			 
		// Adding \sigma_i y_i <= 1 constraints.
		GRBLinExpr expr = 0;
		for (long i = 0; i < components.size(); i++)
		{
			expr += Y[i];
		}
		model.addConstr(expr <= 1);



		// Adding X[v] <= Y[i] constraints.
		for (long i = 0; i < components.size(); i++)
		{
			for (long j = 0; j < components[i].size(); j++)
			{
				long v = components[i][j];  // v is a vertex in component i
				model.addConstr(X[v] <= Y[i]);
			}
		}
		model.update();

		// Providing initial solution
		for (long i = 0; i<g.n; i++)
			X[i].set(GRB_DoubleAttr_Start, 0);

		for (long i = 0; i < HeuristicSolution.size(); i++)
		{
			long v = HeuristicSolution[i];
			X[v].set(GRB_DoubleAttr_Start, 1);
		}

		model.optimize();

		double bestLB = model.get(GRB_DoubleAttr_ObjVal);
		double bestUB = model.get(GRB_DoubleAttr_ObjBound);
		cout << bestLB << " " << bestUB << " ";

		int status = model.get(GRB_IntAttr_Status);
		if (status == GRB_OPTIMAL)
		{
			subOpt = false;
			for (long l = 0; l<g.n; l++)
				if (X[l].get(GRB_DoubleAttr_X)>0.5)
					S.push_back(l);
		}
	}
	catch (GRBException e) {
		cout << "Error code = " << e.getErrorCode() << endl;
		cout << e.getMessage() << endl;
	}
	catch (...) {
		cout << "Error during optimization" << endl;
	}
	return S;
}

vector<long> solveMaxKClub_CutLike(KGraph &g, long k, vector<long> &HeuristicSolution, bool &subOpt)
{
	vector<long> S;
	subOpt = true;
	vector< vector< long> > components;
	vector<long> degreeZero;
	g.FindConnectedComponents(components, degreeZero);
	long lb = HeuristicSolution.size();

	try
	{
		GRBEnv env = GRBEnv();
		env.set(GRB_IntParam_Method, 3); //concurrent 
		env.set(GRB_IntParam_OutputFlag, 0);
		env.set(GRB_DoubleParam_TimeLimit, 3600);
		env.set(GRB_IntParam_LazyConstraints, 1);
		env.set(GRB_DoubleParam_MIPGap, 0.0);
		GRBModel model = GRBModel(env);
		GRBVar *X = model.addVars(g.n, GRB_BINARY);
		GRBVar *Y = model.addVars(components.size(), GRB_BINARY);
		model.update();


		// Adding objective function
		for (long w = 0; w < g.n; w++)
			X[w].set(GRB_DoubleAttr_Obj, 1);
		model.set(GRB_IntAttr_ModelSense, GRB_MAXIMIZE);
		model.update();

		// Fixing singletons to zero.
		for (long i = 0; i < degreeZero.size(); i++)
		{
			long v = degreeZero[i];
			model.addConstr(X[v] == 0);
		}

		// Fixing y[i]=0 when |V(G_i)| < lb.
		for (long i = 0; i < components.size(); i++)
		{
			if (components[i].size() < lb)
			{
				model.addConstr(Y[i] == 0);
			}
		}

		// Adding x_u + x_v <= 1 constraints
		// only need to add such constraints when u and v belong to same component.
		for (long i = 0; i < components.size(); i++)
		{
			// add all constraints within component i
			for (long j = 0; j < components[i].size(); j++)
			{
				long v = components[i][j]; // v is a vertex in component i.
				vector<long> dist = g.ShortestPathsUnweighted(v);
				for (long q = j + 1; q < components[i].size(); q++)
				{
					long u = components[i][q];  // u is another vertex in component i
												// vertices u and v belong to same component i
												// if dist(v,u) > s, then add constraint x_u + x_v <= 1.
					if (dist[u] > k)
					{
						model.addConstr(X[u] + X[v] <= 1);
					}
				}
			}
		}
 		model.update();


		// Adding \sigma_i y_i <= 1 constraints.
		GRBLinExpr expr = 0;
		for (long i = 0; i < components.size(); i++)
		{
			expr += Y[i];
		}
		model.addConstr(expr <= 1);
		model.update();


		// Adding X[v] <= Y[i] constraints.
		for (long i = 0; i < components.size(); i++)
		{
			for (long j = 0; j < components[i].size(); j++)
			{
				long v = components[i][j];  // v is a vertex in component i
				model.addConstr(X[v] <= Y[i]);
			}
		}
		model.update();


		//Adding lazy constraints.
		Kclub_callback cb = Kclub_callback(X, g, k);
		model.setCallback(&cb);

		// Providing initial solution
		for (long i = 0; i<g.n; i++)
			X[i].set(GRB_DoubleAttr_Start, 0);

		for (long i = 0; i < HeuristicSolution.size(); i++)
		{
			long v = HeuristicSolution[i];
			X[v].set(GRB_DoubleAttr_Start, 1);
		}

		model.optimize();

		long NumOfBandBNodes = (long)model.get(GRB_DoubleAttr_NodeCount);
		cerr << "# B&B nodes = " << NumOfBandBNodes << endl;
		cerr << "# callbacks = " << Kclub_callback::numCallbacks << endl;
		cerr << "# lazy cuts = " << Kclub_callback::numLazyCuts << endl;

		double bestLB = model.get(GRB_DoubleAttr_ObjVal);
		double bestUB = model.get(GRB_DoubleAttr_ObjBound);
		cout << bestLB << " " << bestUB << " ";

		int status = model.get(GRB_IntAttr_Status);
		if (status == GRB_OPTIMAL)
		{
			subOpt = false;
			for (long i = 0; i < g.n; i++)
				if (X[i].get(GRB_DoubleAttr_X) > 0.5)
					S.push_back(i);
		}

		delete[] X;
	}

	catch (GRBException e) {
		cout << "Error code = " << e.getErrorCode() << endl;
		cout << e.getMessage() << endl;
	}
	catch (...) {
		cout << "Error during optimization" << endl;
	}
	return S;
}


vector<long> ICUT_Subproblem(KGraph &g, long k, long v_i, long LowerBound, long &SubProblemLB, long &SubProblemUB)
{
	vector<long> S;
	vector< vector< long> > components;
	vector<long> degreeZero;
	g.FindConnectedComponents(components, degreeZero);

	try
	{
		GRBEnv env = GRBEnv();
		env.set(GRB_IntParam_Method, 3); //concurrent 
		env.set(GRB_IntParam_OutputFlag, 0);
		env.set(GRB_DoubleParam_TimeLimit, 3600);
		env.set(GRB_IntParam_LazyConstraints, 1);
		env.set(GRB_DoubleParam_MIPGap, 0.0);
		env.set(GRB_DoubleParam_Cutoff, LowerBound);
		GRBModel model = GRBModel(env);
		GRBVar *X = model.addVars(g.n, GRB_BINARY);
		GRBVar *Y = model.addVars(components.size(), GRB_BINARY);
		model.update();


		// Adding objective function
		for (long w = 0; w < g.n; w++)
			X[w].set(GRB_DoubleAttr_Obj, 1);
		model.set(GRB_IntAttr_ModelSense, GRB_MAXIMIZE);
		model.update();

		// Fixing singletons to zero.
		for (long i = 0; i < degreeZero.size(); i++)
		{
			long v = degreeZero[i];
			model.addConstr(X[v] == 0);
		}

		// Fixing y[i]=0 when |V(G_i)| < lb.
		for (long i = 0; i < components.size(); i++)
		{
			if (components[i].size() < LowerBound)
			{
				model.addConstr(Y[i] == 0);
			}
		}

		// Adding x_u + x_v <= 1 constraints
		// only need to add such constraints when u and v belong to same component.
		for (long i = 0; i < components.size(); i++)
		{
			// add all constraints within component i
			for (long j = 0; j < components[i].size(); j++)
			{
				long v = components[i][j]; // v is a vertex in component i.
				vector<long> dist = g.ShortestPathsUnweighted(v);
				for (long q = j + 1; q < components[i].size(); q++)
				{
					long u = components[i][q];  // u is another vertex in component i
												// vertices u and v belong to same component i
												// if dist(v,u) > s, then add constraint x_u + x_v <= 1.
					if (dist[u] > k)
					{
						model.addConstr(X[u] + X[v] <= 1);
					}
				}
			}
		}
		model.update();

		// Adding \sigma_i y_i <= 1 constraints.
		GRBLinExpr expr = 0;
		for (long i = 0; i < components.size(); i++)
		{
			expr += Y[i];
		}
		model.addConstr(expr <= 1);
		model.update();


		// Adding X[v] <= Y[i] constraints.
		for (long i = 0; i < components.size(); i++)
		{
			for (long j = 0; j < components[i].size(); j++)
			{
				long v = components[i][j];  // v is a vertex in component i
				model.addConstr(X[v] <= Y[i]);
			}
		}

		// Fixing X[v_i] to 1.
		model.addConstr(X[v_i] == 1);

		model.update();


		//Adding lazy constraints.
		Kclub_callback cb = Kclub_callback(X, g, k);
		model.setCallback(&cb);


		model.optimize();


		SubProblemLB = model.get(GRB_DoubleAttr_ObjVal);
		SubProblemUB = model.get(GRB_DoubleAttr_ObjBound);
		
		int status = model.get(GRB_IntAttr_Status);
		if (status == GRB_CUTOFF)
		{
			SubProblemLB = LowerBound;
			SubProblemUB = LowerBound;
			return S;
		}
				
		else if (status == GRB_OPTIMAL)
		{
			SubProblemLB = model.get(GRB_DoubleAttr_ObjVal);
			SubProblemUB = model.get(GRB_DoubleAttr_ObjBound);
			for (long i = 0; i < g.n; i++)
				if (X[i].get(GRB_DoubleAttr_X) > 0.5)
					S.push_back(i);
		}

		delete[] X;
	}

	catch (GRBException e) {
		cout << "Error code = " << e.getErrorCode() << endl;
		cout << e.getMessage() << endl;
	}
	catch (...) {
		cout << "Error during optimization" << endl;
	}
	return S;
}


void ICUT(KGraph &g, long k, vector<long> &BestKClub)
{

	time_t start = clock();
	KGraph gk = g.CreatePowerGraph(k);

	long BestLowerBound = BestKClub.size();
	long BestUpperBound = BestKClub.size();

	vector<long> degeneracyordering = gk.FindDegeneracyOrdering();
	vector<bool> T(g.n, false);

	for (long i = g.n - 1; i >= 0; i--)
	{
		vector<long> SubproblemVertices;

		long v = degeneracyordering[i];
		T[v] = true;
		vector<long> dist_from_v = g.ShortestPathsUnweighted(v, T);

		for (long j = 0; j < g.n; j++)
			if (dist_from_v[j] <= k)
				SubproblemVertices.push_back(j);

		if (SubproblemVertices.size() <= BestLowerBound) continue;

		vector<long> Rmap;
		KGraph g_subproblem = g.CreateInducedGraph(SubproblemVertices,Rmap);

		long SubProblemLB;
		long SubProblemUB;
		vector<long> SubproblemSolution = ICUT_Subproblem(g_subproblem, k, Rmap[v], BestLowerBound, SubProblemLB, SubProblemUB);

		BestLowerBound = max(BestLowerBound, SubProblemLB);
		BestUpperBound = max(BestUpperBound, SubProblemUB);
		 
		if (SubproblemSolution.size() >= BestLowerBound)
		{
			BestKClub = SubproblemSolution;
			for (long q = 0; q < BestLowerBound; q++)
			{
				long v = BestKClub[q];
				long w = SubproblemVertices[v];
				BestKClub[q] = w;   
			}
		}
	}
	cout << "bestLB = " << BestLowerBound << ", bestUB = " << BestUpperBound << " ";
}


vector<long> MaxKclubRevisedVeremyevFormulation(KGraph &g1, long s, vector<long> &hsoln, bool &subOpt)
{
	vector<long> D;
	subOpt = true;

	//create induced graph g
	vector<long> nonzero;
	for (long i = 0; i < g1.n; i++)
	{
		if (g1.degree[i] > 0) nonzero.push_back(i);
	}
	vector<long> Rmap;
	KGraph g = g1.CreateInducedGraph(nonzero, Rmap);

	try {

		long start_formulation_creation = clock();

		GRBEnv env = GRBEnv();
		env.set(GRB_IntParam_OutputFlag, 0);
		env.set(GRB_IntParam_Method, 3);  //concurrent
		env.set(GRB_DoubleParam_TimeLimit, 3600);
		env.set(GRB_DoubleParam_MIPGap, 0.0);
		GRBModel model = GRBModel(env);
		GRBVar *X = model.addVars(g.n, GRB_BINARY);
		GRBVar ***Y = new GRBVar**[g.n];   // Y[i][j][k] denotes whether there is a path from i to j of length k+1 (in the solution)
		for (long i = 0; i<g.n; i++)
		{
			GRBVar **Y_temp = new GRBVar*[g.n];
			for (long j = 0; j<g.n; j++)
			{
				Y_temp[j] = model.addVars(s, GRB_BINARY);
			}
			Y[i] = Y_temp;
		}
		model.update();
		
		// Adding objective function. 
		for (long i = 0; i<g.n; i++)
			X[i].set(GRB_DoubleAttr_Obj, 1);
		model.set(GRB_IntAttr_ModelSense, GRB_MAXIMIZE);
		model.update();

		// Adding constraints for t=1.
		for (long i = 0; i<g.n; i++)
		{
			vector<bool> neighbors(g.n, false);
			for (long k = 0; k<g.degree[i]; k++)
			{
				long j = g.adj[i][k];
				neighbors[j] = true;
			}
			for (long j = 0; j<g.n; j++)
			{
				if (neighbors[j])
				{
					model.addConstr(X[i] + X[j] <= Y[i][j][0] + 1);	// constraint (1)
					model.addConstr(Y[i][j][0] <= X[i]);			// constraint (2)
					model.addConstr(Y[i][j][0] <= X[j]);			// constraint (3)
				}
				else
				{
					Y[i][j][0].set(GRB_DoubleAttr_UB, 0);
				}
			}
		}
		
		// Adding constraints for t>1.
		for (long t = 1; t<s; t++)
		{
			for (long k = 0; k<g.n; k++)
			{
				for (long v = 0; v<g.degree[k]; v++)
				{
					long j = g.adj[k][v];
					// now we have an arc (j,k)
					for (long i = 0; i<g.n; i++)
					{
						if (i == j || i == k) continue;
						model.addConstr(Y[i][j][t - 1] + X[k] <= Y[i][k][t] + 1);	// constraint (4)
					}
				}
			}
		}
		for (long t = 1; t<s; t++)
		{
			for (long k = 0; k<g.n; k++)
			{
				for (long i = 0; i<g.n; i++)
				{
					if (i == k) continue;
					model.addConstr(Y[i][k][t] <= X[k]);	// constraint (5)

					GRBLinExpr expr = 0;
					for (long v = 0; v<g.degree[k]; v++)
					{
						long j = g.adj[k][v];
						expr += Y[i][j][t - 1];
					}
					model.addConstr(Y[i][k][t] <= expr);	// constraint (6)
				}
			}
		}
		for (long i = 0; i<g.n; i++)
		{
			for (long j = 0; j<g.n; j++)
			{
				if (i == j) continue;
				GRBLinExpr expr1 = 0;
				for (long t = 0; t<s; t++)
				{
					expr1 += Y[i][j][t];
				}
				model.addConstr(X[i] + X[j] <= 1 + expr1);  // constraint (7)
			}
		}
		// Since graph is undirected, we are fixing y^t_{ij} = y^t_{ji}.
		for (long i = 0; i<g.n; i++)
		{
			for (long j = i + 1; j<g.n; j++)
			{
				for (long t = 0; t<s; t++)
				{
					model.addConstr(Y[i][j][t] == Y[j][i][t]);
				}
			}
		}
		model.update();

		
		// Providing initial solution.
		for (long i = 0; i<g.n; i++)
			X[i].set(GRB_DoubleAttr_Start, 0);
		for (long i = 0; i < hsoln.size(); i++)
		{
			long v = Rmap[hsoln[i]];
			X[v].set(GRB_DoubleAttr_Start, 1);
		}

		long stop_formulation_creation = clock();
		cerr << "formulation creation time: " << (stop_formulation_creation - start_formulation_creation) / double(CLOCKS_PER_SEC) << endl;

		cerr << "optimizing" << endl;
		model.optimize();

		double bestLB = model.get(GRB_DoubleAttr_ObjVal);
		double bestUB = model.get(GRB_DoubleAttr_ObjBound);
		cout << bestLB << " " << bestUB << " ";

		long status = model.get(GRB_IntAttr_Status);
		cerr << "Status = " << status << endl;
		if (status == GRB_OPTIMAL)
		{
			subOpt = false;
			for (long i = 0; i < g.n; i++)
				if (X[i].get(GRB_DoubleAttr_X) > 0.5)
					D.push_back(i);
		}

	}
	catch (GRBException e) {
		cout << "Error code = " << e.getErrorCode() << endl;
		cout << e.getMessage() << endl;
	}
	catch (...) {
		cout << "Error during optimization" << endl;
	}
	return D;

}


vector<long> solveMaxClubCHC(KGraph &g, long k, vector<long> &HeuristicSolution, bool &subOpt)
{
	vector<long> S;
	subOpt = true;
	vector< vector< long> > components;
	vector<long> degreeZero;
	g.FindConnectedComponents(components, degreeZero);
	long lb = HeuristicSolution.size();

	try
	{
		GRBEnv env = GRBEnv();
		env.set(GRB_IntParam_Method, 3); //concurrent 
		env.set(GRB_IntParam_OutputFlag, 0);
		env.set(GRB_DoubleParam_TimeLimit, 3600);
		env.set(GRB_IntParam_LazyConstraints, 1);
		env.set(GRB_DoubleParam_MIPGap, 0.0);
		GRBModel model = GRBModel(env);
		GRBVar *X = model.addVars(g.n, GRB_BINARY);
		GRBVar *Y = model.addVars(components.size(), GRB_BINARY);
		model.update();
		cerr << "# connected components = " << components.size() << endl;
		cerr << "# non-isolated nodes = " << g.n - degreeZero.size() << endl;
		cerr << "# remaining edges = " << g.m << endl;

		// Adding objective function.
		for (long w = 0; w < g.n; w++)     //
			X[w].set(GRB_DoubleAttr_Obj, 1);
		model.set(GRB_IntAttr_ModelSense, GRB_MAXIMIZE);
		model.update();

		// Fixing singletons to zero.
		for (long i = 0; i < degreeZero.size(); i++)
		{
			long v = degreeZero[i];
			model.addConstr(X[v] == 0);
		}

		// Fixing y[i]=0 when |V(G_i)| < lb.
		for (long i = 0; i < components.size(); i++)
		{
			if (components[i].size() < lb)
			{
				model.addConstr(Y[i] == 0);
			}
		}

		// Adding x_u + x_v <= 1 constraints.
		// only need to add such constraints when u and v belong to same component.
		for (long i = 0; i < components.size(); i++)
		{
			// add all constraints within component i
			for (long j = 0; j < components[i].size(); j++)
			{
				long v = components[i][j]; // v is a vertex in component i.
				vector<long> dist = g.ShortestPathsUnweighted(v);
				for (long q = j + 1; q < components[i].size(); q++)
				{
					long u = components[i][q];  // u is another vertex in component i
												// vertices u and v belong to same component i
												// if dist(v,u) > s, then add constraint x_u + x_v <= 1.
					if (dist[u] > k)
					{
						model.addConstr(X[u] + X[v] <= 1);
					}
				}
			}
		}
		model.update();

		// Adding \sigma_i y_i <= 1 constraints.
		GRBLinExpr expr = 0;
		for (long i = 0; i < components.size(); i++)
		{
			expr += Y[i];
		}
		model.addConstr(expr <= 1);
		model.update();


		// Adding X[v] <= Y[i] constraints.
		for (long i = 0; i < components.size(); i++)
		{
			for (long j = 0; j < components[i].size(); j++)
			{
				long v = components[i][j];  // v is a vertex in component i
				model.addConstr(X[v] <= Y[i]);
			}
		}
		model.update();

		// Adding lazy constraints.
		CHC cb = CHC(X, g, k);
		model.setCallback(&cb);

		// Providing initial solution.
		for (long i = 0; i<g.n; i++)
			X[i].set(GRB_DoubleAttr_Start, 0);

		for (long i = 0; i < HeuristicSolution.size(); i++)
		{
			long v = HeuristicSolution[i];
			X[v].set(GRB_DoubleAttr_Start, 1);
		}
		model.optimize();

		double bestLB = model.get(GRB_DoubleAttr_ObjVal);
		double bestUB = model.get(GRB_DoubleAttr_ObjBound);
		cout << bestLB << " " << bestUB << " ";

		long NumOfBandBNodes = (long)model.get(GRB_DoubleAttr_NodeCount);
		cerr << "# B&B Nodes = " << NumOfBandBNodes << endl;
		cerr << "# callbacks = " << CHC::numCallbacks << endl;

		int status = model.get(GRB_IntAttr_Status);
		if (status == GRB_OPTIMAL)
		{
			subOpt = false;
			for (long i = 0; i < g.n; i++)
				if (X[i].get(GRB_DoubleAttr_X) > 0.5)
					S.push_back(i);
		}

		delete[] X;
	}
	catch (GRBException e) {
		cout << "Error code = " << e.getErrorCode() << endl;
		cout << e.getMessage() << endl;
	}
	catch (...) {
		cout << "Error during optimization" << endl;
	}
	return S;
}



vector<long> Pathlike_3Club(KGraph &g, vector<long> &hsoln, bool &subOpt)
{
	vector<long> D;
	long k = 3;
	subOpt = true;
	vector< vector< long> > components;
	vector<long> degreeZero;
	g.FindConnectedComponents(components, degreeZero);
	long lb = hsoln.size();

	map<vector<long>, long, classcomp>map = EnumerateLength3Connector(g);
	long mapsize = 0;

	try
	{
		long start_formulation_creation = clock();

		GRBEnv env = GRBEnv();
		env.set(GRB_IntParam_Method, 3); //concurrent 
		env.set(GRB_IntParam_OutputFlag, 0);
		env.set(GRB_DoubleParam_TimeLimit, 3600);
		env.set(GRB_DoubleParam_MIPGap, 0.0);
		GRBModel model = GRBModel(env);
		GRBVar *X = model.addVars(g.n, GRB_BINARY);
		GRBVar *Y = model.addVars(components.size(), GRB_BINARY);
		GRBVar *Z = model.addVars(map.size(), GRB_BINARY);
		model.update();
		cerr << "# connected components = " << components.size() << endl;
		cerr << "# non-isolated nodes = " << g.n - degreeZero.size() << endl;
		cerr << "# remaining edges = " << g.m << endl;

		// Adding objective function.
		for (long w = 0; w < g.n; w++)     
			X[w].set(GRB_DoubleAttr_Obj, 1);
		model.set(GRB_IntAttr_ModelSense, GRB_MAXIMIZE);
		model.update();

		// Fixing singletons to zero.
		for (long i = 0; i < degreeZero.size(); i++)
		{
			long v = degreeZero[i];
			model.addConstr(X[v] == 0);
		}

		// Fixing y[i]=0 when |V(G_i)| < lb.
		for (long i = 0; i < components.size(); i++)
		{
			if (components[i].size() < lb)
			{
				model.addConstr(Y[i] == 0);
			}
		}

		// Adding x_u + x_v <= 1 constraints.
		// only need to add such constraints when u and v belong to same component.
		for (long i = 0; i < components.size(); i++)
		{
			// add all constraints within component i
			for (long j = 0; j < components[i].size(); j++)
			{
				long v = components[i][j]; // v is a vertex in component i.
				vector<long> dist = g.ShortestPathsUnweighted(v);
				for (long q = j + 1; q < components[i].size(); q++)
				{
					long u = components[i][q];  // u is another vertex in component i
												// vertices u and v belong to same component i
												// if dist(v,u) > s, then add constraint x_u + x_v <= 1.
					if (dist[u] == 1) continue;
					else if (dist[u] > k)  // dist > 3
					{
						// do nothing. rhs will be one.
						model.addConstr(X[u] + X[v] <= 1);
					}
					else
					{
						GRBLinExpr rhs = 1;
						vector<long> cn = g.CommonNeighborsList(u, v);
						for (long p = 0; p < cn.size(); p++)
						{
							long q = cn[p];
							rhs += X[q];
						}

						vector<long> dist_from_u = g.ShortestPathsUnweighted(u);
						vector<long> dist_from_v = dist;
						for (long u_neighbors_iterator = 0; u_neighbors_iterator < g.adj[u].size(); u_neighbors_iterator++)
						{
							long ii = g.adj[u][u_neighbors_iterator];
							if (dist_from_v[ii] == 2)
							{
								for (long ii_neighbors_iterator = 0; ii_neighbors_iterator < g.adj[ii].size(); ii_neighbors_iterator++)
								{
									long jj = g.adj[ii][ii_neighbors_iterator];
									if (dist_from_u[jj] == 2 && dist_from_v[jj] == 1)
									{
										long minij = min(ii, jj);
										long maxij = max(ii, jj);
										std::map<std::vector<long>, long>::iterator it = map.find({ minij,maxij });
										rhs += Z[it->second];
									}
								}
							}
						}
						model.addConstr(X[u] + X[v] <= rhs);
					}
				}
			}
		}

		// Adding coupling constraints.
		for (std::map<std::vector<long>, long>::iterator it = map.begin(); it != map.end(); ++it)
		{
			for (int i = 0; i < it->first.size(); i++)
			{
				model.addConstr(Z[it->second] <= X[it->first[i]]);
			}
		}
		model.update();

		// Adding \sigma_i y_i <= 1 constraints.
		GRBLinExpr expr = 0;
		for (long i = 0; i < components.size(); i++)
		{
			expr += Y[i];
		}
		model.addConstr(expr <= 1);
		model.update();


		// Adding X[v] <= Y[i] constraints.
		for (long i = 0; i < components.size(); i++)
		{
			for (long j = 0; j < components[i].size(); j++)
			{
				long v = components[i][j];  // v is a vertex in component i
				model.addConstr(X[v] <= Y[i]);
			}
		}
		model.update();

		// Providing initial solution.
		for (long i = 0; i<g.n; i++)
			X[i].set(GRB_DoubleAttr_Start, 0);
		for (long i = 0; i < hsoln.size(); i++)
		{
			long v = hsoln[i];
			X[v].set(GRB_DoubleAttr_Start, 1);
		}

		long stop_formulation_creation = clock();
		cerr << "formulation creation time: " << (stop_formulation_creation - start_formulation_creation) / double(CLOCKS_PER_SEC) << endl;

		model.optimize();

		double bestLB = model.get(GRB_DoubleAttr_ObjVal);
		double bestUB = model.get(GRB_DoubleAttr_ObjBound);
		cout << bestLB << " " << bestUB << " ";

		long NumBBNodes = (long)model.get(GRB_DoubleAttr_NodeCount);
		cerr << "# B&B nodes = " << NumBBNodes << endl;

		int status = model.get(GRB_IntAttr_Status);
		if (status == GRB_OPTIMAL)
		{
			subOpt = false;
			for (long i = 0; i < g.n; i++)
				if (X[i].get(GRB_DoubleAttr_X) > 0.5)
					D.push_back(i);
		}

		delete[] X;
	}
	catch (GRBException e) {
		cout << "Error code = " << e.getErrorCode() << endl;
		cout << e.getMessage() << endl;
	}
	catch (...) {
		cout << "Error during optimization" << endl;
	}
	return D;
}



vector<long> Pathlike_4Club(KGraph &g, vector<long> &hsoln, bool &subOpt)
{
	vector<long> D;
	long k = 4;
	subOpt = true;
	vector< vector< long> > components;
	vector<long> degreeZero;
	g.FindConnectedComponents(components, degreeZero);
	long lb = hsoln.size();

	map<vector<long>, long, classcomp>map1 = EnumerateLength3Connector(g);
	cerr << "# length 3 connectors = " << map1.size() << endl;
	map<vector<long>, long, classcomp>map2 = EnumerateLength4Connector(g);
	cerr << "# length 4 connectors = " << map2.size() << endl;

	try
	{
		long start_formulation_creation = clock();

		GRBEnv env = GRBEnv();
		env.set(GRB_IntParam_Method, 3); //concurrent 
		env.set(GRB_IntParam_OutputFlag, 0);
		env.set(GRB_DoubleParam_TimeLimit, 3600);
		env.set(GRB_IntParam_LazyConstraints, 1);
		env.set(GRB_DoubleParam_MIPGap, 0.0);
		GRBModel model = GRBModel(env);
		GRBVar *X = model.addVars(g.n, GRB_BINARY);
		GRBVar *Y = model.addVars(components.size(), GRB_BINARY);
		GRBVar *Z = model.addVars(map1.size(), GRB_CONTINUOUS);
		GRBVar *W = model.addVars(map2.size(), GRB_CONTINUOUS);
		model.update();
		cerr << "# connected components = " << components.size() << endl;
		cerr << "# non-isolated nodes = " << g.n - degreeZero.size() << endl;
		cerr << "# remaining edges = " << g.m << endl;

		// cerr << "Adding objective function.
		for (long w = 0; w < g.n; w++)     //
			X[w].set(GRB_DoubleAttr_Obj, 1);
		model.set(GRB_IntAttr_ModelSense, GRB_MAXIMIZE);
		model.update();

		// Fixing singletons to zero.
		for (long i = 0; i < degreeZero.size(); i++)
		{
			long v = degreeZero[i];
			model.addConstr(X[v] == 0);
		}

		// cerr << Fixing y[i]=0 when |V(G_i)| < lb.
		for (long i = 0; i < components.size(); i++)
		{
			if (components[i].size() < lb)
			{
				model.addConstr(Y[i] == 0);
			}
		}

		// Adding x_u + x_v <= 1 constraints.
		// only need to add such constraints when u and v belong to same component.
		for (long i = 0; i < components.size(); i++)
		{
			// add all constraints within component i
			for (long j = 0; j < components[i].size(); j++)
			{
				long v = components[i][j]; // v is a vertex in component i.
				vector<long> dist = g.ShortestPathsUnweighted(v);
				for (long q = j + 1; q < components[i].size(); q++)
				{
					long u = components[i][q];  // u is another vertex in component i
												// vertices u and v belong to same component i
												// if dist(v,u) > s, then add constraint x_u + x_v <= 1.
					if (dist[u] == 1) continue;
					else if (dist[u] > k)  // dist > 4
					{
						// do nothing. rhs will be one.
						model.addConstr(X[u] + X[v] <= 1);
					}
					else
					{
						GRBLinExpr rhs = 1;
						// add 2-hop connectors
						vector<long> cn = g.CommonNeighborsList(u, v);
						for (long p = 0; p < cn.size(); p++)
						{
							long q = cn[p];
							rhs += X[q];
						}

						// add 3-hop connecotrs
						vector<long> dist_from_u = g.ShortestPathsUnweighted(u);
						if (dist_from_u[v] == 2 || dist_from_u[v] == 3)
						{
							vector<long> dist_from_v = dist;
							for (long u_neighbors_iterator = 0; u_neighbors_iterator < g.adj[u].size(); u_neighbors_iterator++)
							{
								long ii = g.adj[u][u_neighbors_iterator];
								if (dist_from_v[ii] == 2)
								{
									for (long ii_neighbors_iterator = 0; ii_neighbors_iterator < g.adj[ii].size(); ii_neighbors_iterator++)
									{
										long jj = g.adj[ii][ii_neighbors_iterator];
										if (dist_from_u[jj] == 2 && dist_from_v[jj] == 1)
										{
											long minij = min(ii, jj);
											long maxij = max(ii, jj);
											std::map<std::vector<long>, long>::iterator it = map1.find({ minij,maxij });
											if (it == map1.end()) cerr << "ERROR: Could not find this path in the map.\n";
											rhs += Z[it->second];
										}
									}
								}
							}
						}

						// add 4-hop connectors
						vector<long> dist_from_v = dist;
						//type1
						for (long u_neighbors_iterator = 0; u_neighbors_iterator < g.adj[u].size(); u_neighbors_iterator++)
						{
							long ii = g.adj[u][u_neighbors_iterator];
							if (dist_from_v[ii] == 3)
							{
								for (long ii_neighbors_iterator = 0; ii_neighbors_iterator < g.adj[ii].size(); ii_neighbors_iterator++)
								{
									long jj = g.adj[ii][ii_neighbors_iterator];
									if (dist_from_u[jj] == 2 && dist_from_v[jj] == 2)
									{
										for (long jj_neighbors_iterator = 0; jj_neighbors_iterator < g.adj[jj].size(); jj_neighbors_iterator++)
										{
											long kk = g.adj[jj][jj_neighbors_iterator];
											if (dist_from_u[kk] == 2 && dist_from_v[kk] == 1)
											{
												vector<long> sortedsubset = sortnodes(ii, jj, kk);
												std::map<std::vector<long>, long>::iterator it = map2.find(sortedsubset);
												if (it == map2.end()) cerr << "ERROR: Could not find this path in the map.\n";
												rhs += W[it->second];
											}
										}
									}
								}
							}
						}
						//type 2
						for (long u_neighbors_iterator = 0; u_neighbors_iterator < g.adj[u].size(); u_neighbors_iterator++)
						{
							long ii = g.adj[u][u_neighbors_iterator];
							if (dist_from_v[ii] == 2)
							{
								for (long ii_neighbors_iterator = 0; ii_neighbors_iterator < g.adj[ii].size(); ii_neighbors_iterator++)
								{
									long jj = g.adj[ii][ii_neighbors_iterator];
									if (dist_from_u[jj] == 2 && dist_from_v[jj] == 2)
									{
										for (long jj_neighbors_iterator = 0; jj_neighbors_iterator < g.adj[jj].size(); jj_neighbors_iterator++)
										{
											long kk = g.adj[jj][jj_neighbors_iterator];
											if (dist_from_u[kk] == 3 && dist_from_v[kk] == 1)
											{
												vector<long> sortedsubset = sortnodes(ii, jj, kk);
												std::map<std::vector<long>, long>::iterator it = map2.find(sortedsubset);
												if (it == map2.end()) cerr << "ERROR: Could not find this path in the map.\n";
												rhs += W[it->second];
											}
										}
									}
								}
							}
						}
						//type 3
						for (long u_neighbors_iterator = 0; u_neighbors_iterator < g.adj[u].size(); u_neighbors_iterator++)
						{
							long ii = g.adj[u][u_neighbors_iterator];
							if (dist_from_v[ii] == 3)
							{
								for (long ii_neighbors_iterator = 0; ii_neighbors_iterator < g.adj[ii].size(); ii_neighbors_iterator++)
								{
									long jj = g.adj[ii][ii_neighbors_iterator];
									if (dist_from_u[jj] == 2 && dist_from_v[jj] == 2)
									{
										for (long jj_neighbors_iterator = 0; jj_neighbors_iterator < g.adj[jj].size(); jj_neighbors_iterator++)
										{
											long kk = g.adj[jj][jj_neighbors_iterator];
											if (dist_from_u[kk] == 3 && dist_from_v[kk] == 1)
											{
												vector<long> sortedsubset = sortnodes(ii, jj, kk);
												std::map<std::vector<long>, long>::iterator it = map2.find(sortedsubset);
												if (it == map2.end()) cerr << "ERROR: Could not find this path in the map.\n";
												rhs += W[it->second];
											}
										}
									}
								}
							}
						}
						//type 4
						for (long u_neighbors_iterator = 0; u_neighbors_iterator < g.adj[u].size(); u_neighbors_iterator++)
						{
							long ii = g.adj[u][u_neighbors_iterator];
							if (dist_from_v[ii] == 2)
							{
								// boolify neighborhood of i
								vector<bool> N_ii(g.n, false);
								for (long NeighborNumber = 0; NeighborNumber < g.degree[ii]; NeighborNumber++)
								{
									N_ii[g.adj[ii][NeighborNumber]] = true;
								}

								for (long ii_neighbors_iterator = 0; ii_neighbors_iterator < g.adj[ii].size(); ii_neighbors_iterator++)
								{
									long jj = g.adj[ii][ii_neighbors_iterator];
									if (dist_from_u[jj] == 2 && dist_from_v[jj] == 2)
									{
										for (long jj_neighbors_iterator = 0; jj_neighbors_iterator < g.adj[jj].size(); jj_neighbors_iterator++)
										{
											long kk = g.adj[jj][jj_neighbors_iterator];
											if (N_ii[kk] == false)
											{
												if (dist_from_u[kk] == 2 && dist_from_v[kk] == 1)
												{
													vector<long> sortedsubset = sortnodes(ii, jj, kk);
													std::map<std::vector<long>, long>::iterator it = map2.find(sortedsubset);
													if (it == map2.end()) cerr << "ERROR: Could not find this path in the map.\n";
													rhs += W[it->second];
												}
											}
										}
									}
								}
							}
						}
						model.addConstr(X[u] + X[v] <= rhs);
					}
				}
			}
		}
		model.update();

		// Adding coupling constraints.
		for (std::map<std::vector<long>, long>::iterator it = map1.begin(); it != map1.end(); ++it)
		{
			for (long i = 0; i < it->first.size(); i++)
			{
				model.addConstr(Z[it->second] <= X[it->first[i]]);
			}
		}


		for (std::map<std::vector<long>, long>::iterator it = map2.begin(); it != map2.end(); ++it)
		{
			for (long i = 0; i < it->first.size(); i++)
			{
				model.addConstr(W[it->second] <= X[it->first[i]]);
			}
		}


		// Adding \sigma_i y_i <= 1 constraints.
		GRBLinExpr expr = 0;
		for (long i = 0; i < components.size(); i++)
		{
			expr += Y[i];
		}
		model.addConstr(expr <= 1);
		model.update();


		// Adding X[v] <= Y[i] constraints.
		for (long i = 0; i < components.size(); i++)
		{
			for (long j = 0; j < components[i].size(); j++)
			{
				long v = components[i][j];  // v is a vertex in component i
				model.addConstr(X[v] <= Y[i]);
			}
		}
		model.update();


		// Providing initial solution.
		for (long i = 0; i<g.n; i++)
			X[i].set(GRB_DoubleAttr_Start, 0);
		for (long i = 0; i < hsoln.size(); i++)
		{
			long v = hsoln[i];
			X[v].set(GRB_DoubleAttr_Start, 1);
		}

		long stop_formulation_creation = clock();
		cerr << "formulation creation time: " << (stop_formulation_creation - start_formulation_creation) / double(CLOCKS_PER_SEC) << endl;

		model.optimize();

		double bestLB = model.get(GRB_DoubleAttr_ObjVal);
		double bestUB = model.get(GRB_DoubleAttr_ObjBound);
		cout << bestLB << " " << bestUB << " ";

		long NumBBNodes = (long)model.get(GRB_DoubleAttr_NodeCount);
		cerr << "# B&B nodes are: " << NumBBNodes << endl;

		int status = model.get(GRB_IntAttr_Status);
		if (status == GRB_OPTIMAL)
		{
			subOpt = false;
			for (long i = 0; i < g.n; i++)
				if (X[i].get(GRB_DoubleAttr_X) > 0.5)
					D.push_back(i);
		}

		delete[] X;
	}
	catch (GRBException e) {
		cout << "Error code = " << e.getErrorCode() << endl;
		cout << e.getMessage() << endl;
	}
	catch (...) {
		cout << "Error during optimization" << endl;
	}
	return D;
}
