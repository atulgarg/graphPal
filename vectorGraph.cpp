#include "vectorGraph.h"
#include<iostream>
#include<fstream>
#include<omp.h>
#include<cmath>
#include<limits>
#include<cstdlib>
/**
 * @method add to add a edge to existing nodes in a graph.
 * Method adds edges considering the graph is undirected.
 * @param Edge edge which needs to be added.
 *
 */
void VectorGraph::add(Edge edge)
{
    graph[edge.node1].push_back(Edge(edge.node1,edge.node2,edge.cost));
    graph[edge.node2].push_back(Edge(edge.node2,edge.node1,edge.cost));
}
/**
 * @method initialiseGraph to initialise Graph with total 
 * node size of passed parameter and declare an array of 
 * vectors of size equal to total_nodes.
 * @param int total number of nodes in graph.
 *
 */
void VectorGraph::initialiseGraph(int total_nodes)
{
    this->total_nodes = total_nodes;
    graph = new vector<Edge>[total_nodes];
}
/**
 * Constructor which in turn calls initialiseGraph method 
 * to initialise total number of nodes in graph.
 *
 */
VectorGraph::VectorGraph(int total_nodes)
{
    initialiseGraph(total_nodes);
}
/**
 * Constructor to initialise a graph reading input from a 
 * file which is stored as form of 
 * First line >> Total Nodes.
 * Further lines >> node1 node2 cost.
 * @param char* name of the file to read.
 *
 */
VectorGraph::VectorGraph(char* filename)
{
    ifstream infile(filename);
    infile>>total_nodes;

    initialiseGraph(total_nodes);
    int source, destination;
    double cost;
    //read node1, node2 and cost. store undirected graphs.
    while(infile>>source>>destination>>cost)
    {
        graph[source].push_back(Edge(source,destination,cost));
        graph[destination].push_back(Edge(destination,source,cost));
    }

    infile.close();
}
/**
 * @method print() to print Graph as adjacency list. Prints
 * each of the list in one line for each
 * of the nodes. Used for testing
 *
 */
void VectorGraph::print()
{
    for(int i =0;i<total_nodes;i++)
    {
        cout<<i<<":: ";
        for(int j=0;j<graph[i].size();j++)
        {
            cout<<graph[i].at(j).node2<<", "<<graph[i].at(j).cost<<" ";
        }
        cout<<endl;
    }
}
/**
 * @method VectorGraph::mst_seq_prim to implement sequential Prim's 
 * implementation.
 *
 */
VectorGraph* VectorGraph::mst_seq_prim()
{
    VectorGraph* vg = new VectorGraph(total_nodes);
    priority_queue<Edge, vector<Edge>, EdgeCompare> pq;   
    double total_cost = 0.0;
    int start_node = 0;
    double start_time = omp_get_wtime();
    //number of nodes covered.
    int num_nodes = 1;
    bool visited[total_nodes];
    //push all the edges for source in pq.
    for(int i=0;i<graph[start_node].size();i++)
    {
        pq.push(graph[start_node].at(i));
    }

    //this and the other loop can be executed in parallel.
    for(int j=0;j<total_nodes;j++)
        visited[j] = false;

    visited[start_node] = true;
    //while all the nodes are not covered in Minimum Spanning Tree.
    while(num_nodes <total_nodes)
    {
        Edge minEdge = pq.top();
        pq.pop();
        if(!visited[minEdge.node1] || !visited[minEdge.node2])
        {
            vg->add(minEdge);
            total_cost +=minEdge.cost;
            int node;
            //add edges for second node.
            if(!visited[minEdge.node1])
            {
                node = minEdge.node1;
                visited[minEdge.node1] = true;
            }
            else
            {
                node = minEdge.node2;
                visited[minEdge.node2] = true;
            }
            //push edges from this graph to priority queue.
            for(int k = 0;k<graph[node].size();k++)
            {
                pq.push(graph[node].at(k));
            }
            num_nodes++;
        }
    }
    double end_time = omp_get_wtime();
    cout<<"Total Nodes :: "<<total_nodes<<endl;
    cout<<"Time Prim's Sequential Algorithm Implementation :: "
		<<end_time-start_time<<endl<<endl;
    return vg;

}
/**
 * @method mst_p_kruskal to implement concurrent implementation for prim
 * This method does not have much of concurrency constructs and hence
 * is not discussed in the report and was not mentioned in the poster.
 * This implementation also did not give enough of speed up compared 
 * with other implementation.
 *
 */
VectorGraph* VectorGraph::mst_p_kruskal()
{
    VectorGraph* vg = new VectorGraph(total_nodes);
    int num_proc = omp_get_max_threads() - 1;
    FibonacciHeap<Edge> pq[total_nodes];
    double total_cost = 0.0;
    int start_node = 0;
    int num_nodes = 1;
    double start_time = omp_get_wtime();
    //number of nodes covered.
    bool visited[total_nodes];

    visited[start_node] = true;
   
    //initialise priority queue for each of the nodes
    //schedule is selected as dynamic since number of edges 
    //in vertex may vary.

#pragma omp parallel for schedule(dynamic, 4) 
    for(int i=0;i<total_nodes;i++)
    {
        for(int j=0;j<graph[i].size();j++)
        {
            pq[i].insert(graph[i].at(j));
        }
    }
    cout<<omp_get_wtime()-start_time<<endl;
   
    //Only after all threads have completed initialisation 
   //of there particular respective queues.

#pragma omp barrier
    while(num_nodes<total_nodes)
    {
        if(!pq[start_node].isEmpty())
        {
            Edge minEdge = pq[start_node].removeMinimum();
            if(!visited[minEdge.node2] || !visited[minEdge.node1])
            {
                if(!visited[minEdge.node2])
                {
                    visited[minEdge.node2] = true;
                    pq[start_node].merge(pq[minEdge.node2]);
                }
                else
                {
                    visited[minEdge.node1] = true;
                    pq[start_node].merge(pq[minEdge.node1]);
                }
                num_nodes++;
                total_cost+=minEdge.cost;
                vg->add(minEdge);
            }
        }else
            break;
    }

    double end_time = omp_get_wtime();
    cout<<"Total Nodes :: "<<total_nodes<<endl;
    cout<<"Time Parallel Prim's Algorithm Implementation"
		<<end_time-start_time<<endl<<endl;
    return vg;
}

/***
 * @method mst_p1() parallel implementation of Boruvka' Algorithm
 * @return Graph* minimum spanning tree initialised.
 *
 */
VectorGraph* VectorGraph::mst_p_boruvka_pq()
{
    VectorGraph* vg = new VectorGraph(total_nodes);
    int num_components = total_nodes;
    priority_queue<Edge, vector<Edge>, EdgeCompare> pq[num_components];

    set<int> sets[num_components];
    int set_lookup[num_components];
    double start_time = omp_get_wtime();
    double total_cost = 0.0;

#pragma omp parallel for
    for(int i=0;i<num_components;i++)
    {
        //create set of single vertex each and 
	//add edge corresponding to each of the vertex's
        sets[i].insert(i);
        for(int j=0;j<graph[i].size();j++)
            pq[i].push(graph[i].at(j));
        set_lookup[i] =i;
    }
    //#pragma omp barrier
    while(num_components > 1)
    {
        //While there are more than 1 components in a 
	//connected graph keep iterating.For each of 
	//the sets find the minimum weighted edge which
	// is out side the component. This can be done 
	//in parallel for each of the component.

        Edge* edge_array[num_components];
#pragma omp parallel for
        for(int i=0;i<num_components;i++)
        {
            edge_array[i] = NULL;
    
            //find the min weight edge from pq for the set.
            //looped here since the edge may belong inside the component.
            
            while(!pq[i].empty())
            {
                Edge edge = pq[i].top();
                pq[i].pop();
            
                //if this is edge going out of the component insert the edge
                if(sets[i].find(edge.node2) == sets[i].end())
                {
                    edge_array[i] = new Edge(edge.node1, edge.node2, edge.cost);
                    break;
                }
            }
        }
        double temp_time = omp_get_wtime();
        int x = num_components;
        
	//Once all the min edges for the components are known add them
        //if there exists no cycle and self loops..
#pragma omp barrier
        for(int i=0;i<x;i++)
        {
            if(edge_array[i] != NULL)
            {
                int set_node1 = set_lookup[edge_array[i]->node1];
                int set_node2 = set_lookup[edge_array[i]->node2];
                if(set_node1!=set_node2)
                {
                    int min = std::min(set_node1, set_node2);
                    int max = std::max(set_node1, set_node2);

                    //the below can be done in parallel tasks.
                    //also can be decided on the basis if size of max set 
                    //has some threshold elements.

#pragma omp parallel sections
                    {
#pragma omp section
                        {
                            //Section for merging the two sets.
                            mergesets(sets,set_lookup, min, max);
                            vg->add(*edge_array[i]);
                            total_cost+=edge_array[i]->cost;
                        }
#pragma omp section
                        {
                            //Section for merging respective priority queues.
                            joinPQ(pq, min, max);
                        }
                    }

                    //Decreament the total number of components remaining
                    //since two are merged.
                    num_components--;
                }
            }
        }
    }
    double end_time = omp_get_wtime();
    cout<<"Total Nodes :: "<<total_nodes<<endl;
    cout<<"Parallel Boruvka's Algorithm PQ Implementation Total Time :: "
		<<end_time-start_time<<" Total Cost:: "<<total_cost<<endl;
    return vg;
}
/**
 * @method mergesets to merge two given sets.
 * @param set<int> sets array of sets in which merge in needs to be done.
 * @param int[] integer array for look up.
 * @param min set to which merging needs to be done.
 * @param max set from which merging needs to be done. 
 *
 */
void VectorGraph::mergesets(set<int> sets[],int set_lookup[], int min, int max)
{
    set<int>::iterator set_iter = sets[max].begin();
    for(;set_iter!=sets[max].end();set_iter++)
    {
        sets[min].insert(*set_iter);
        set_lookup[*set_iter] = set_lookup[min];
    }
}
/**
 * @method joinPQ to join Priority Queue for two sets.
 * @param priority_queue<Edge,vector<Edge>,EdgeCompare> array or priority queue.
 * @param min index to which priority queue needs to be merged.
 * @param max index of which priority queue needs to be merged.
 */
void VectorGraph::joinPQ(priority_queue<Edge, vector<Edge>,EdgeCompare>* pq, int min, int max)
{
    while(!pq[max].empty())
    {
        pq[min].push(pq[max].top());
        pq[max].pop();
    }
}
/**
 * @method to merge Fibonacci heaps belonging to each of the sets.
 * @param array of FibonacciHeaps
 * @param heap index to which merge needs to be done.
 * @param heap index from where merging needs to be done.
 *
 */
void VectorGraph::joinPQ(FibonacciHeap<Edge>* pq, int min, int max)
{
        pq[min].merge(pq[max]);
}
/**
 * @method mst_seq_boruvka method to evaluate Minimum spanning tree 
 * using Boruvka's algorithm and Fibonacci heap as priority queue.
 * @returns VectorGraph* pointer to Graph represented as MST.
 */
VectorGraph* VectorGraph::mst_seq_boruvka()
{
    VectorGraph* vg = new VectorGraph(total_nodes);
    int num_components = total_nodes;
    FibonacciHeap<Edge> pq[num_components];

    set<int> sets[num_components];
    int set_lookup[num_components];
    double start_time = omp_get_wtime();
    double total_cost = 0.0;
    for(int i=0;i<num_components;i++)
    {
        //create set of single vertex each and add edge corresponding 
        //to each of the vertex's
        sets[i].insert(i);
        for(int j=0;j<graph[i].size();j++)
            pq[i].insert(graph[i].at(j));
        set_lookup[i] =i;
    }

    while(num_components > 1)
    {
        //While there are more than 1 components in a 
	//connected graph keep iterating. For each of 
	//the sets find the minimum weighted edge which
	//is out side the component. This can be done in
	//parallel for each of the component.

        Edge* edge_array[num_components];
        for(int i=0;i<num_components;i++)
        {
            edge_array[i] = NULL;
          
	    //find the min weight edge from pq for the set.
            //looped here since the edge may belong inside the component.
            while(!pq[i].isEmpty())
            {
                Edge edge = pq[i].removeMinimum();
                //if this is edge going out of the component insert the edge
                if(sets[i].find(edge.node2) == sets[i].end())
                {
                    edge_array[i] = new Edge(edge.node1, edge.node2, edge.cost);
                    break;
                }
            }
        }
        //copied in another variable since should iterate till last component.
        int x = num_components;
        
        //Once all the min edges for the components are known add 
        //them if there exists no cycle and self loops..
        for(int i=0;i<x;i++)
        {
            if(edge_array[i] != NULL)
            {
                int set_node1 = set_lookup[edge_array[i]->node1];
                int set_node2 = set_lookup[edge_array[i]->node2];
                if(set_node1!=set_node2)
                {
                    int min = std::min(set_node1, set_node2);
                    int max = std::max(set_node1, set_node2);
                    
                    //merge the sets for nodes.
                    mergesets(sets,set_lookup, min, max);

                    //merge the respective priority Queue
                    joinPQ(pq, min, max);
                    vg->add(*edge_array[i]);
                    total_cost+=edge_array[i]->cost;
                    num_components--;
                }
            }
        }
    }
    double end_time = omp_get_wtime();
    cout<<"Total Nodes :: "<<total_nodes<<endl;
    cout<<"Time for Sequential Boruvka' Algorithm: "<<end_time-start_time<<endl;
    return vg;
}
/**
 * @method mst_p_boruvka_fibo to determine Minimum 
 * spanning tree using concurrent algorithm of Boruvka.
 * Concurrency is exploited in parallel for loop
 *
 */
VectorGraph* VectorGraph::mst_p_boruvka_fibo()
{
    VectorGraph* vg = new VectorGraph(total_nodes);
    int num_components = total_nodes;
    FibonacciHeap<Edge> pq[num_components];

    //array of sets for each node.
    set<int> sets[num_components];
    
    //lookup array to save time.
    int set_lookup[num_components];
    double start_time = omp_get_wtime();
    double total_cost = 0.0;

    //distribute edges to each priority queue equally to 
    //share work load.
    int div = ceil((double)num_components/omp_get_max_threads());
    
    //iterate this for loop in separate threads independently.
#pragma omp parallel for schedule(static, div) 
    for(int i=0;i<num_components;i++)
    {
        //create set of single vertex each and add edge 
        //corresponding to each of the vertex's
       
        sets[i].insert(i);
       
        //initialise priority queue for each vertex.
        for(int j=0;j<graph[i].size();j++)
            pq[i].insert(graph[i].at(j));
        set_lookup[i] =i;
    }

#pragma omp barrier
    //iterate till number of components are more than 1.
    while(num_components > 1)
    {
        //While there are more than 1 components in a connected 
        //graph keep iterating. For each of the sets find the 
        //minimum weighted edge which is out side the component. 
        //This can be done in parallel for each of the component.
    
        Edge* edge_array[num_components];
#pragma omp parallel for schedule(dynamic,4)
        for(int i=0;i<num_components;i++)
        {
            edge_array[i] = NULL;
            
            //find the min weight edge from pq for the set.
            //looped here since the edge may belong inside the component.
            while(!pq[i].isEmpty())
            {
                Edge edge = pq[i].removeMinimum();

                //if this is edge going out of the component insert the edge
                if(sets[i].find(edge.node2) == sets[i].end())
                {
                    edge_array[i] = new Edge(edge.node1, edge.node2, edge.cost);
                    break;
                }
            }
        }
        int x = num_components;
        //Once all the min edges for the components are known add 
        //them if there exists no cycle and self loops.
#pragma omp barrier
        for(int i=0;i<x;i++)
        {
            if(edge_array[i] != NULL)
            {
                int set_node1 = set_lookup[edge_array[i]->node1];
                int set_node2 = set_lookup[edge_array[i]->node2];
                if(set_node1!=set_node2)
                {
                    int min = std::min(set_node1, set_node2);
                    int max = std::max(set_node1, set_node2);

                    //the below can be done in parallel tasks. Also can be 
                    //dynamically decided on the basis if size of max set 
                    //has some threshold elements. Since the number of 
                    //parallel sections is 2 set the max active levels to 2.
                    omp_set_max_active_levels(2);
#pragma omp parallel sections
                    {
#pragma omp section
                        {
                            //Section for merging the two sets respective to 
                            //each node.
                            mergesets(sets,set_lookup, min, max);
                        }
#pragma omp section
                        {
                            //Merge Fibonacci heap corresponding to
                            //each of the component in parallel.
                            pq[min].merge(pq[max]);
                            vg->add(*edge_array[i]);
                            total_cost+=edge_array[i]->cost;
                        }
                    }
                    num_components--;
                }
            }
        }
    }
    double end_time = omp_get_wtime();
    cout<<end_time-start_time<<endl;
    cout<<"Total Nodes :: "<<total_nodes<<endl;
    cout<<"Time Boruvka's Parallel Algorithm Implementation:: "
		<<end_time-start_time<<endl<<endl;
    return vg;
}
/**
 * @method mst_seq_kruskal to determine Minimum spanning tree for
 * given Graph using Kruskal's algorithm. Method uses Fibonacci 
 * heap as priority queue for best performances.
 * @returns VectorGraph* pointer to graph initialised with 
 * minimum spanning tree.
 *
 */
VectorGraph* VectorGraph::mst_seq_kruskal()
{
    VectorGraph* mst = new VectorGraph(total_nodes);
    double cost = 0.0;
    double start_time = omp_get_wtime();
    int num_nodes = 1;
    set<int> sets[total_nodes];
    int node_lookup[total_nodes];
    FibonacciHeap<Edge> pq;  

    //Initialisation of priority queue.
    for(int i=0;i<total_nodes;i++)
    {
        for(int j=0;j<graph[i].size();j++)
        {
            pq.insert(graph[i].at(j));
        }
        sets[i].insert(i);
        node_lookup[i] = i;
    }
    
    //Iterate till all the forests merge in one.
    while(num_nodes < total_nodes && !pq.isEmpty())
    {
        //remove the minimum weighted edge.
        Edge edge = pq.removeMinimum();
        int min = std::min(edge.node1,edge.node2);
        int max = std::max(edge.node1,edge.node2);

        //if the nodes does not belong to same set i.e.
        //no cycles are formed.
        if(node_lookup[min]!=node_lookup[max])
        {
            mst->add(edge);
            num_nodes++;
            mergesets(sets,node_lookup, min, max);
        }
    }
    double end_time = omp_get_wtime();
    cout<<"Total Nodes :: "<<total_nodes<<endl; 
    cout<<"Time Kruskal's Algorithm Sequential Implementation :: "
		<<end_time - start_time<<endl<<endl; 
    return mst;
}
/**
 * @method mst_p_kruskal to implement concurrent version for 
 * Kruskal's algorithm. Method aims to concurrently determine
 * Minimum spanning tree for the given Graph represented as 
 * array of vectors. Each of the iterations in initialisation 
 * is done in Open MP parallel for loops. Multiple Priority 
 * queues equal to total number of processors are used to store
 * edges of Graph. The edges of graph are distributed over each
 * queue equally in number. Array of Minimum Edge over each
 * queue is maintained in iteration (edgew). Over each iteration
 * the task of merging node's sets and initialising next edgew 
 * array index from the corresponding iteration id done in parallel.
 * @returns VertexGraph* Graph representation for Minimum spanning tree.
 *
 */
VectorGraph* VectorGraph::mst_p_kruskal1()
{
    VectorGraph* mst = new VectorGraph(total_nodes);
    int num_nodes = 1;
    set<int> sets[total_nodes];
    int node_lookup[total_nodes];
    int num_proc = omp_get_max_threads();
    FibonacciHeap<EdgeWrapper> pq[num_proc];  
    
    double start_time = omp_get_wtime();
    int count = 0; 
    int div = total_nodes/(float)omp_get_max_threads();
    
    //Initialisation for distributing edges across priority queues 
    //is done in parallel. Schedule is choose to be static so as each
    //each thread is assigned equal number of iterations depending on
    //number of processors. Chunk size of div equal to total num 
    //nodes divided by total number of threads is used to distribute work load equally.
#pragma omp parallel for schedule(static,div)
    for(int i=0;i<total_nodes;i++)
    {
        for(int j=0;j<graph[i].size();j++)
        {
            pq[count%div].insert(EdgeWrapper(graph[i].at(j),count%num_proc));
            count++;
        }
        sets[i].insert(i);
        node_lookup[i] = i;
    }
    
    //Barrier so as each of the threads have intialised respective priority queue.
#pragma omp barrier
    EdgeWrapper* edgew[num_proc];

    //Parallel for loop construct to initialise 
    //edgew[i] corresponding to ith priority queue.
#pragma omp parallel for schedule(dynamic, 4)
    for(int i =0;i<num_proc;i++)
    {
        edgew[i] = NULL;
        if(!pq[i].isEmpty())
            edgew[i] = new EdgeWrapper(pq[i].removeMinimum());
    }

    Edge edge(edgew[0]->edge);
    int last_pq = 0;

    //Iterate till all the nodes are not discovered in MST.
    while(num_nodes < total_nodes)
    {
        //find min of the edges. 
        //This can be implemented using reduction if cost 
        //of each edge is distinct and reverse mapping for 
        //cost to edge is used, which will enable use of 
        //parallel for loop construct.
        for(int i=0;i<num_proc;i++)
        {
            if(edgew[i]!=NULL && edgew[i]->edge < edge)
            {
                edge = edgew[i]->edge;
                last_pq = edgew[i]->pq_num;
            }
        }
        
	//Total number of threads set to 2 
        //since total parallel sections are two which 
        //needs to work simultaneously.
        omp_set_max_active_levels(2);
#pragma omp parallel sections
        { 
            //Parallel sections while  one section works to
            //merge the two sets, other section initialises the last 
#pragma omp section
            {
                //Section merging the two sets of source 
                //and destination. Min and max are determined to keep
                //a ordering in merge and hence correctness.
                int min = std::min(edge.node1,edge.node2);
                int max = std::max(edge.node1,edge.node2);
                if(node_lookup[min]!=node_lookup[max])
                {
                    mst->add(edge);
                    num_nodes++;
                    mergesets(sets,node_lookup, min, max);
                }
            }
#pragma omp section
            {
                //Section to initialise edgew element with 
                //next minimum from the priority queue.
                if(!pq[last_pq].isEmpty())
                    edgew[last_pq] = new EdgeWrapper(pq[last_pq].removeMinimum());
                else
                    edgew[last_pq] = NULL;
                edge.cost = numeric_limits<double>::max();
            }
        }
    }
    double end_time = omp_get_wtime();
    cout<<"Total Nodes :: "<<total_nodes<<endl;
    cout<<"Time Kruskal's Algorithm Parallel Implementation :: "
			<<end_time-start_time<<endl;
    return mst;

}
/**
 * @method dijkstras to implement Dijkstra's Sequential Algorithm 
 * to determine single source shortest path for provided source.
 * @param int identifier for source.
 * @returns total weight of minimum distance's from source vertex
 *
 */
int VectorGraph::dijkstras(int source)
{
    bool visited[total_nodes]; 
    int min_dist[total_nodes];
    int p=0;
    int u;  

    priority_queue<Edge, vector<Edge>, EdgeCompare> pq;
    double start_time = omp_get_wtime();
    //initialising minimum distance array for all nodes
    for (int i = 0; i < total_nodes; i++){
        min_dist[i] = std::numeric_limits<int>::max();
        visited[i] = false;
    }
    //filling up the pq with adjacent nodes fo the source
    for(int i=0;i<graph[source].size();i++){
        pq.push(graph[source].at(i));
    }
    min_dist[source] = 0;
    visited[source] = true;
    for (int j = 0; j < total_nodes;j++){

        if(pq.size()!=0){
            Edge minEdge = pq.top();
            pq.pop();
            u = minEdge.node2;
            if(min_dist[u]>min_dist[minEdge.node1]+minEdge.cost)
                min_dist[u] = min_dist[minEdge.node1]+minEdge.cost;
            visited[u] = true;
            p=0;
            //adding the min among neighbours to the pq and updating
            //the minimum idstance array for all relevant nodes
            for (int k = 0; k< graph[u].size(); k++)
            {
                p = graph[u].at(k).node2;
                if (!visited[p] && min_dist[u] !=std::numeric_limits<int>::max() 
                        && min_dist[u]+graph[u].at(k).cost < min_dist[p])
                {
                    min_dist[p] = min_dist[u] + graph[u].at(k).cost;
                    for(int w=0;w<graph[p].size();w++)
                    {
                        pq.push(graph[p].at(w));
                    }
                }
            }
        }
    }

    double end_time = omp_get_wtime();
    cout<<"sequential Dijstras Implementation Time :: "<<end_time-start_time<<endl;
    cout<<min_dist[total_nodes-1]<<endl;
    return min_dist[total_nodes-1];
}
/**
 * @method p_dijkstras to implement Dijkstras's algorithm for 
 * single source shortest path algorithm.
 * @param int identifier for source vertex.
 * @returns total weight of all minimum distances.
 *
 */
int VectorGraph::p_dijkstras(int source){
	 bool visited[total_nodes];
	     int min_dist[total_nodes];
	     int p=0;
	     int u;
	     int local_first;
	     int local_last;
	     int local_min;
	     int local_min_location;
	     int global_min;
	   	 int global_min_location;
	     int local_step;
	     int thread_id;
	     int thread_num;
	     int thread_no;
	     #pragma omp parallel
	    { thread_no= omp_get_num_threads ( );
	    	}
	    	#pragma omp barrier
	     int large_num =std::numeric_limits<int>::max();
	      priority_queue<Edge, vector<Edge>, EdgeCompare> pq[thread_no];
		   double start_time = omp_get_wtime();
		   Edge min_edge = Edge(large_num,large_num,large_num);
		    Edge local_edge = Edge(large_num,large_num,large_num);
		   
		 //initialising the minimum distance array for all the nodes
	     for (int i = 0; i < total_nodes; i++){
	        min_dist[i] = large_num;
	         visited[i] = false;
		}
		
		//pushing all the neighbours of source into the priority queue
		for(int i=0;i<graph[source].size();i++){
			
				pq[0].push(graph[source].at(i));
			
		}
	     min_dist[source] = 0;
	     visited[source] = true;

		# pragma omp parallel private ( local_first, thread_id , local_edge,\
                local_last ,local_min, local_min_location, local_step) \
                shared ( visited, global_min, min_dist, global_min_location,\
                        min_edge, thread_num )
         {
             //receiving current thread id
             thread_id = omp_get_thread_num ( );
             thread_num= omp_get_num_threads ( );

             local_first =   (    thread_id      * total_nodes ) /thread_num;
             if(thread_id!= thread_num-1){
                 local_last  =   (( ( thread_id+ 1 ) * total_nodes ) / thread_num )- 1;}
             else{
                 local_last = total_nodes -1;
             }

#pragma omp barrier
             for ( local_step = 0; local_step <total_nodes; local_step++ )
             {	
                 //initialising global minimum
# pragma omp single
                 {
                     global_min = large_num ;
                     global_min_location = -1;
                     min_edge = Edge(large_num,large_num,large_num);
                 }

# pragma omp barrier		

                 //updating local min
                 //since while the global minimum is updated no other threads should interfere  
#pragma omp critical
                 {	  
                     if(pq[thread_id].size()!=0){
                         if (pq[thread_id].top().cost< global_min)
                         {  global_min = pq[thread_id].top().cost;
                             global_min_location = pq[thread_id].top().node2 ;
                             min_edge =pq[thread_id].top();
                             local_edge = pq[thread_id].top();

                         }
                     }
                 }


# pragma omp barrier

                 if(local_edge.cost==min_edge.cost
                         && local_edge.node1==min_edge.node1
                         && local_edge.node2==min_edge.node2)
                 {
                     pq[thread_id].pop();
                     if(min_dist[min_edge.node2] > min_dist[min_edge.node1]+min_edge.cost)
                         min_dist[min_edge.node2] = min_dist[min_edge.node1]+min_edge.cost;
                     visited[min_edge.node2] = true;
                 }

# pragma omp barrier
                 {
                     if(local_first<= global_min_location 
                             && global_min_location<=local_last
                             && global_min!=large_num)
                     {
                         int u = min_edge.node2;
                         for (int k = 0; k< graph[u].size(); k++)
                         {
                             p = graph[u].at(k).node2;
                             if (!visited[p] && min_dist[u] !=large_num 
                                     && min_dist[u]+graph[u].at(k).cost < min_dist[p])
                             {
                                 min_dist[p] = min_dist[u] + graph[u].at(k).cost;
                                 for(int w=0;w<graph[p].size();w++)
                                 {
                                     for(int l=0;l<thread_num;l++)
                                     {
                                         int temp = std::min(p,w);
                                         if((l*total_nodes/thread_num <=temp 
                                                     && temp<(l+1)*total_nodes/thread_num 
                                                     && l!=thread_num-1)
                                                 ||((l==thread_num-1) 
                                                     &&(l*total_nodes/thread_num <=temp)
                                                     &&(temp<total_nodes)))
                                         {
                                             pq[l].push(graph[p].at(w));
                                         }
                                     }

                                 }
                             }
                         }

                     }	
                 }
#pragma omp barrier
             }
         }
	  double end_time = omp_get_wtime();
	  cout<<"parallel Dijstras Implementation Time :: "<<end_time-start_time<<endl;
	  cout<<min_dist[total_nodes-1]<<endl;
	 return min_dist[total_nodes-1];
}
int main(int argc,char *argv[])
{
    if(argc !=2)
    {
        cerr<<"Incorrect Usage"<<endl<<"Correct Usage: ./prog input_file"<<endl;
        exit(0);
    }
    VectorGraph g(argv[1]);
    int choice = 0;
    do
    {
        cout<<endl<<endl<<"Select Algorithm to Run :: "<<endl
		<<"1. Sequential Kruskal's Implementation"<<endl
		<<"2. Concurrent Kruskal's Implementation "<<endl
		<<"3. Sequential Boruvka's Implementation"<<endl
		<<"4. Concurrent Boruvka's Implementation"<<endl
		<<"5. Sequential Prim's Implementation"<<endl
        <<"6. Sequential Dijkstra's Implementation"<<endl
        <<"7. Concurrent Dijkstra's Implementation" <<endl;
        VectorGraph * tempG;
        cout<<endl<<endl<<"Selection :: ";
        cin>>choice;
        switch(choice)
        {
            case 1:
               tempG = g.mst_seq_kruskal();
                break;
            case 2:
                tempG = g.mst_p_kruskal();
                break;
            case 3:
                tempG = g.mst_seq_boruvka();
                break;
            case 4:
                tempG = g.mst_p_boruvka_fibo();
                break;
            case 5:
                tempG = g.mst_seq_prim();
                break;
            case 6:
                g.dijkstras(0);
                break;
            case 7:
                g.p_dijkstras(0);
                break;
            default:
                cout<<"Incorrect Choice"<<endl;
                break;
        }
    }while(choice > 0 && choice < 8);
    return 0;
}
