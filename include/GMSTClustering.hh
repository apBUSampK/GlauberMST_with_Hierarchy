#include <iostream>
#include <vector>
#include <memory>
#include "../TGlauber/TGlauberMC.hh"
#include "../TGlauber/TGlauNucleon.hh"
#include "TObject.h"
using namespace std;

// Creating shortcut for an integer pair 
typedef  pair<int, int> iPair;

struct GNode {
    pair<shared_ptr<GNode>, shared_ptr<GNode>> children;
    double height;
    int* V;
    int size;

    //constructor from children nodes
    GNode(shared_ptr<GNode>, shared_ptr<GNode>, double);

    //initial nodes constructor
    GNode(int);

    GNode(const GNode&);
    GNode& operator=(const GNode&);

    //destructor
    ~GNode();
};


//the tree handler. Isn't intended to be constructed outside the Graph methods.
class GTree {
    public:
    //constructor
    GTree(int);

    GTree(const GTree&);
    GTree& operator= (const GTree&);

    //merge 2 nodes
    void merge(int, int, double);

    //get clusterization
    vector<GNode> get_cluster(double CD);

    //get node
    inline const GNode* get_node(int a) {return nodes[a].get();}

    //get size
    int get_size() {return size;};

    //width print (for testing purposes)
    friend ostream& operator<< (ostream&, const GTree&);

    private:
    int size;
    vector<shared_ptr<GNode>> nodes;
};

// Structure to represent a graph 
struct Graph
{
	// Vert and edges
    int V, E;
    vector< pair<double, iPair> > edges;
    vector<vector<double>> adj;

    // Constructors 
    Graph(int V, int E);
    Graph();

    // Destructor
	~Graph();

    // Utility function to add an edge
    void addEdge(int u, int v, double w);
    // Function to find MST using Kruskal's 
    // MST algorithm
    // hierarchical algorithms
    GTree AdvancedKruskalMST_Dendro();
};

class GMSTCluster{

	public: 
	GMSTCluster(int Z_in, int A_in);
	~GMSTCluster();

	public: 
	inline int GetZ() {return Z;};
	inline int GetA() {return A;};
	inline int SetZ(int Z_in) {Z = Z_in;}
	inline int SetA(int A_in) {A = A_in;}

	private: 
	int Z;
	int A;
};

typedef std::vector<GMSTCluster> GMSTClusterVector;

class GMSTClustering{

	public:
	GMSTClustering();
	GMSTClustering(double CD_in, double single_silh);
	~GMSTClustering();

	inline double SetCD(double CD_in) {CritDist = CD_in;}
    void SetUp(TObjArray*);
	
    GTree GetTree();
	GMSTClusterVector GetClusters();
    GMSTClusterVector GetClusters_HSilhouette();
	private:
	//Work with graphs is up to Nepeyvoda Roman, even data types. I prefer GraphToCluster and ClusterToGraph methods to be private
	double CritDist, single_silh;
    int A{};
    Graph g;
    TObjArray* nucleons{};
    Graph ClusterToGraph();
    GMSTClusterVector CompileVector(const vector<vector<int>>&);
};