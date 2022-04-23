#include "../include/GMSTClustering.hh"
#include <queue>


bool cd_comp(const GNode& left, const GNode& right) {
    return left.height > right.height;
}


GMSTCluster::GMSTCluster(int Z_in, int A_in): 
A(A_in), Z(Z_in){}

GMSTCluster::~GMSTCluster() = default;

GMSTClustering::GMSTClustering(){
CritDist = 0;
}

GMSTClustering::GMSTClustering(double CD_in, double single_silh, double variation) :
    CritDist(CD_in), single_silh(single_silh), variation(variation) {}

GMSTClustering::~GMSTClustering(){
    delete nucleons;
}

Graph GMSTClustering::ClusterToGraph(){
    Graph g(A, A*(A-1)/2);
    //  making full graph of nucleons
    for(int iArray = 0; iArray < nucleons->GetEntries(); iArray++){
    	TGlauNucleon *nucleon=(TGlauNucleon*)(nucleons->At(iArray));
    	for(int iArray_pairs = iArray + 1; iArray_pairs < nucleons->GetEntries(); iArray_pairs++){
    		TGlauNucleon *nucleon_pair=(TGlauNucleon*)(nucleons->At(iArray_pairs));
    		g.addEdge(iArray, iArray_pairs, std::sqrt(pow(nucleon->GetX() - nucleon_pair->GetX(),2) + pow(nucleon->GetY() - nucleon_pair->GetY(),2) + pow(nucleon->GetZ() - nucleon_pair->GetZ(),2)));
    	}
	}
    return g;
}

void GMSTClustering::SetUp(TObjArray* nucleons_in){
    A = 0;
	nucleons = new TObjArray(1000); // nucleons from Side A
	for(int iArray = 0; iArray < nucleons_in->GetEntries(); iArray++){
		TGlauNucleon *nucleon=(TGlauNucleon*)(nucleons_in->At(iArray));
		if(nucleon->IsSpectator() && nucleon->IsInNucleusA()){
			A++;
			nucleons->AddLast((TObject*)(nucleons_in->At(iArray)));
		}
	}
    
    g = ClusterToGraph();
    
}

GMSTClusterVector GMSTClustering::CompileVector(const vector<vector<int>>& clusters) {
    GMSTClusterVector output_vector;
    for(int i = 0; i < clusters.size(); ++i) {
    	int Z_clust = 0;
    	int A_clust = 0;
    	for(int j = 0; j < clusters[i].size(); ++j) {
        	TGlauNucleon *nucleon=(TGlauNucleon*)(nucleons->At((clusters[i])[j] - 1));
        	if(nucleon->IsProton())
        	{
        		Z_clust += 1;
        	}
        	A_clust += 1;
        }
        output_vector.push_back(GMSTCluster(Z_clust, A_clust));
    }
	return output_vector;
}

inline GTree GMSTClustering::GetTree() {
    return g.AdvancedKruskalMST_Dendro();
}

GMSTClusterVector GMSTClustering::GetClusters(cut) {
    //check for empty nucleons input
    if(!A)
        return CompileVector(vector<vector<int>>());
    //Get the c-link dendrogram
    auto tr = g.AdvancedKruskalMST_Dendro();

    //get the clustering with corresponding CD:
    vector<GNode> clustering = tr.get_cluster(CritDist);
    //compile output vector
    vector<vector<int>> clusters;
    for (const auto & iter : clustering)
        clusters.emplace_back(vector<int>(iter.V, iter.V + iter.size));
    return CompileVector(clusters);
}

GMSTClusterVector GMSTClustering::GetClusters(silhouette) {
    //check for empty nucleons input
    if(!A)
        return CompileVector(vector<vector<int>>());
    //Get the c-link dendrogram
    auto tr = g.AdvancedKruskalMST_Dendro();
    
    //get the clustering with the biggest applicable CD:
    vector<GNode> current = tr.get_cluster(CritDist * (1.0 + variation));
    sort(current.begin(), current.end(), cd_comp);
    double max_silh = -1;
    vector<GNode> best = current;
    //calculate mean cluster silhouettes until min applicable CD
    //separate nucleon clusters are given silhouette of single_silh
    while (current.front().height > CritDist * (1.0 - (variation > 1 ? 1 : variation))) {
        double silh = 0;
        for (auto iter = current.cbegin(); iter != current.cend(); iter++) {
            if (iter->size > 1)
                for (int i = 0; i < iter->size; i++) {
                    double inner = 0, outer = -1;
                    for (int j = 0; j < iter->size; j++)
                        inner += g.adj[iter->V[i] - 1][iter->V[j] - 1];
                    inner /= (iter->size - 1);
                    for (auto jter = current.cbegin(); jter != current.cend(); jter++)
                        if (iter != jter) {
                            double buff = 0;
                            for (int j = 0; j < jter->size; j++)
                                buff += g.adj[iter->V[i] - 1][jter->V[j] - 1];
                            buff /= jter->size;
                            if (outer < 0 || buff < outer)
                                outer = buff;
                        }
                    silh += ((outer - inner) / max(outer, inner));
                }
            else
                silh += single_silh;
            }
        silh /= A;
        if (silh > max_silh) {
            max_silh = silh;
            best = current;
        }
        //divide the biggest cluster into two
        current.push_back(*current.front().children.first);
        current.push_back(*current.front().children.second);
        current.erase(current.begin());
        sort(current.begin(), current.end(), cd_comp);
    }
    //compile output vector
    vector<vector<int>> clusters;
    for (const auto & iter : best)
        clusters.emplace_back(vector<int>(iter.V, iter.V + iter.size));
	return CompileVector(clusters);
};

GMSTClusterVector GMSTClustering::GetClusters(max_alpha) {
    //check for empty nucleons input
    if(!A)
        return CompileVector(vector<vector<int>>());
    //Get the c-link dendrogram
    auto tr = g.AdvancedKruskalMST_Dendro();

    //get the clustering with the biggest applicable CD:
    vector<GNode> current = tr.get_cluster(CritDist * (1.0 + variation));
    sort(current.begin(), current.end(), cd_comp);
    int max_alpha = 0;
    vector<GNode> best = current;
    //find the cluster with the largest alpha particles count
    while (current.front().height > CritDist * (1.0 - (variation > 1 ? 1 : variation))) {
        int alpha = 0;
        for (auto & iter : current) {
            int z_count = 0;
            for (int i = 0; i < iter.size; i++)
                if (((TGlauNucleon*)nucleons->At(iter.V[i] - 1))->IsProton())
                    z_count++;
            if (z_count == 2)
                alpha++;
        }
        if (alpha > max_alpha) {
            max_alpha = alpha;
            best = current;
        }
        //divide the biggest cluster into two
        current.push_back(*current.front().children.first);
        current.push_back(*current.front().children.second);
        current.erase(current.begin());
        sort(current.begin(), current.end(), cd_comp);
    }
    //compile output vector
    vector<vector<int>> clusters;
    for (const auto & iter : best)
        clusters.emplace_back(vector<int>(iter.V, iter.V + iter.size));
    return CompileVector(clusters);
}

Graph::Graph(int V, int E)
{
    this->V = V;
    this->E = E;
    adj = vector<vector<double>>(V, vector<double>(V, 0));
}

Graph::Graph()
{
    this->V = 0;
    this->E = 0;
    adj = vector<vector<double>>(V, vector<double>(V, 0));
}

Graph::~Graph() = default;

void Graph::addEdge(int u, int v, double w)
{
    edges.push_back({ w, {u, v} });
    adj[u][v] = w;
    adj[v][u] = w;
}

GTree Graph::AdvancedKruskalMST_Dendro()
{
    // Sort edges in increasing order on basis of cost
    sort(edges.begin(), edges.end());

    // Create a tree
    GTree tr(V);

    // Iterate through all sorted edges
    vector< pair<double, iPair> >::iterator it;
    for (it = edges.begin(); it != edges.end(); it++)
    {
        int u = it->second.first;
        int v = it->second.second;

        // Check if the selected edge is creating
        // a cycle or not (Cycle is created if u
        // and v belong to same set)
        if (tr.get_node(u) != tr.get_node(v))
        {
            // create new node
            tr.merge(u, v, it->first);
        }
    }

    return tr;
}


GNode::GNode(int n) {
    size = 1;
    V = new int[1];
    V[0] = n;
    height = 0;
    children = make_pair(nullptr, nullptr);
}

GNode::GNode(shared_ptr<GNode> first, shared_ptr<GNode> second, double height) {
    size = first->size + second->size;
    V = new int[size];
    for (int i = 0; i < size; i++) {
        int buff = (i < first->size) ? first->V[i] : second->V[i - first->size];
        V[i] = buff;
       // cout << V[i] << endl;
    }
    this->height = height;
    children = make_pair(first, second);
}

GNode::GNode(const GNode& right) {
    size = right.size;
    V = new int[size];
    for (int i = 0; i < size; i++)
        V[i] = right.V[i];
    height = right.height;
    if (right.children.first == nullptr)
        children = make_pair(shared_ptr<GNode>(nullptr),shared_ptr<GNode>(nullptr));
    else
        children = make_pair(make_shared<GNode>(*right.children.first), make_shared<GNode>(*right.children.second));
}

GNode& GNode::operator=(const GNode& right) {
    size = right.size;
    delete[] V;
    V = new int[size];
    for (int i = 0; i < size; i++)
        V[i] = right.V[i];
    height = right.height;
    if (right.children.first == nullptr)
        children = make_pair(shared_ptr<GNode>(nullptr),shared_ptr<GNode>(nullptr));
    else
        children = make_pair(make_shared<GNode>(*right.children.first), make_shared<GNode>(*right.children.second));
    return *this;
}

GNode::~GNode() {
    delete[] V;
}

GTree::GTree(int size) {
    this->size = size;
    nodes = vector<shared_ptr<GNode>>(size);
    for (int i = 0; i < size; i++)
        nodes[i] = make_shared<GNode>(i + 1);
}

GTree::GTree(const GTree& right) {
    size = right.size;
    nodes = vector<shared_ptr<GNode>>(size);
    bool* cpd = new bool[size];
    for (int i = 0; i < size; i++)
        cpd[i] = false;
    for (int i = 0; i < size; i++) {
        if (!cpd[i]) {
            nodes[i] = make_shared<GNode>(*right.nodes[i]);
            for (int j = 0; j < nodes[i]->size; j++) {
                cpd[nodes[i]->V[j]] = true;
                nodes[nodes[i]->V[j]] = nodes[i];
            }
        }
    }
    delete[] cpd;
}

GTree& GTree::operator=(const GTree& right) {
    size = right.size;
    nodes = vector<shared_ptr<GNode>>(size);
    bool *cpd = new bool[size];
    for (int i = 0; i < size; i++)
        cpd[i] = false;
    for (int i = 0; i < size; i++) {
        if (!cpd[i]) {
            nodes[i] = make_shared<GNode>(*right.nodes[i]);
            for (int j = 0; j < nodes[i]->size; j++) {
                cpd[nodes[i]->V[j]] = true;
                nodes[nodes[i]->V[j]] = nodes[i];
            }
        }
    }
    delete[] cpd;
    return *this;
}

void GTree::merge(int first, int second, double height) {
    int t = nodes[first]->size + nodes[second]->size;
    int* chngPtr = new int[t];
    for (int i = 0; i < nodes[first]->size; i++)
        chngPtr[i] = nodes[first]->V[i] - 1;
    for (int i = 0; i < nodes[second]->size; i++)
        chngPtr[i + nodes[first]->size] = nodes[second]->V[i] - 1;
    nodes[first] = make_shared<GNode>(nodes[first], nodes[second], height);
    for (int i = 0; i < t; i++)
        nodes[chngPtr[i]] = nodes[first];
    delete[] chngPtr;
}

vector<GNode> GTree::get_cluster(double CD) {
    vector<GNode> cluster;
    queue<GNode*> width;
    width.push(nodes[0].get());
    while(!width.empty()) {
        GNode* node = width.front();
        width.pop();
        if (node->height <= CD)
            cluster.push_back(*node);
        else if (node->children.first != nullptr) {
            width.push(node->children.first.get());
            width.push(node->children.second.get());
        }
    }
    return cluster;
}

ostream& operator<<(ostream& out, const GTree& right) {
    queue<GNode*> width;
    width.push(right.nodes[0].get());
    while (!width.empty()) {
        GNode* node = width.front();
        width.pop();
        out << '(';
        for (int i = 0; i < node->size; i++)
            out << node->V[i] << ", ";
        out << ") height=" << node->height << endl;
        if (node->children.first != nullptr) {
            width.push(node->children.first.get());  
            width.push(node->children.second.get());
        }
    }
    return out;
}
