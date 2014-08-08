#include<vector>
#include<set>
#include<queue>
#include "fibonacci.hpp"
using namespace std;
/*
 * Class Edge used to store edges in Graph. Each edge has source destination
 * as node1 and node2 respectively. Edges are assigned cost stored with each
 * edge.
 *
 */
class Edge
{
    public:
        int node1;
        int node2;
        double cost;
        Edge(){
        }
        Edge(int node1, int node2, double cost)
        {
            this->node1 = node1;
            this->node2 = node2;
            this->cost = cost;
        }
        inline bool operator< (const Edge& rhs)
        {
            if(this->cost != rhs.cost)
                return this->cost<rhs.cost;
            else
                return this->node2 < rhs.node2;
        }
        inline bool operator> (const Edge& rhs)
        {
            if(this->cost != rhs.cost)
                return this->cost>rhs.cost;
            else
                return this->node2 > rhs.node2;
        }
        inline bool static compare(const Edge& lhs,const Edge& rhs)
        {
            if(lhs.cost != rhs.cost)
                return lhs.cost> rhs.cost;
            else
                return lhs.node1> rhs.node2;

        }
};
/**
 * class EdgeWrapper wrapper around edge class used in Kruskal's concurrent execution
 * where each priority queue stores object of type Edgewrapper to keep a track of which 
 * priority queue edge belonged.
 *
 */
class EdgeWrapper
{
    public:
        Edge edge;

        //identifier to priority queue to which edge belongs.
        //
        int pq_num;
        EdgeWrapper()
        {
        }
        EdgeWrapper(Edge edge, int pq_num)
        {
            this->edge = edge;
            this->pq_num = pq_num;
        }
        EdgeWrapper(const EdgeWrapper& ew)
        {
            this->edge = ew.edge;
            this->pq_num = ew.pq_num;
        }
        bool operator()(const EdgeWrapper& lhs, const EdgeWrapper& rhs) const
        {
            return (Edge::compare(lhs.edge, rhs.edge)); 
        }
        inline bool operator< (const EdgeWrapper& rhs)
        {
            return this->edge < rhs.edge;
        }
        inline bool operator> (const EdgeWrapper& rhs)
        {
            return this->edge > rhs.edge;
        }

};
/*
 * Class EdgeCompare declared for comparison of Edges in Priority Queue.
 * It overrides () operator which is used by STL priority queue for comparison.
 *
 */
class EdgeCompare
{
    public:
    bool operator()(const Edge& lhs, const Edge& rhs) const
    {
        if(lhs.cost != rhs.cost)
            return lhs.cost> rhs.cost;
        else
            return lhs.node1> rhs.node2;
    }
};
/**
 * Graph representation. Each graph is stored as array of vectors. Each vector is
 * adjacency list for nodes where sparse matrix optimization is used to store only
 * edges with weights.
 *
 */
class VectorGraph
{
    private:
        vector<Edge>* graph;
        int total_nodes;

        //method to initialise array of vectors.
        void initialiseGraph(int total_nodes);
        
        //method to add edges to graph for existing nodes.
        void add(Edge edge);
        
        void joinPQ(FibonacciHeap<Edge>* pq, int min, int max);
       
        //method to join two priority queues corresponding to different sets.
        void joinPQ(priority_queue<Edge, vector<Edge>,EdgeCompare>* pq, int min, int max);  

        //method to join two sets.
        void mergesets(set<int> sets[],int set_lookup[],int min, int max);

    public:

        VectorGraph(char *filename);
        VectorGraph(int total_nodes);
        void print();
        
        //Sequential Implementation of Prim.
        VectorGraph* mst_seq_prim();

        //MST implemented using Priority queue of STL and Boruvka's Algorithm.
        VectorGraph* mst_p_boruvka_pq();

        //MST using Fibonacci heap implementation and Boruvka's Algorithm.
        VectorGraph* mst_p_boruvka_fibo();

        //Sequential implementation of Boruvka Algorithm.
        VectorGraph* mst_seq_boruvka(); 
        
        //Sequential implementation of Kruskal's Algorithm.
        VectorGraph* mst_seq_kruskal();
        
        //Concurrent Implementation of Kruskal's Algorithm.
        VectorGraph* mst_p_kruskal();
        VectorGraph* mst_p_kruskal1();

        //sequential implementation of dijekstra.
        int dijkstras(int source);
        
        //parallel implementation of dijekstra.
        int p_dijkstras(int source);
};
