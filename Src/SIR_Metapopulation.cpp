#include <fstream>
#include <sstream>
#include <iostream>
#include <RandomNumber.hpp>
#include <vector>
#include <algorithm>
#include <numeric>
#include <cmath>
#include <queue>

using namespace std;
using namespace Snu::Cnrc;
using Size = unsigned int;


class SIR_Metapopulation {
private:
    unsigned int number_of_nodes;
    unsigned int number_of_edges;
    double       kappa;
    double       p;
    unsigned int number_of_ensembles;
    unsigned int uid;
    double       alpha;
    Size         b;
    Size         l;
    
public:
    SIR_Metapopulation(const unsigned, const unsigned, const double, const double, const unsigned, const unsigned, const double, const Size, const Size);
    ~SIR_Metapopulation();
    void Run();
    Size find_subgraph_given_distance(vector<vector<Size>>& graph_of_graph, Size picked_subgraph_1, Size distance, const Size b, const Size l);
    void switch_node(Size picked_site_1, Size distance, Size picked_node_in_subgraph_2, vector<vector<vector<Size>>>& graph, vector<vector<Size>>& state, vector<vector<Size>>& graph_of_graph, const Size b, const Size l, vector<vector<double>>& deltaP, vector<Size>& SI, vector<Size>& SE, vector<Size>& E, vector<Size>& I, vector<Size>& R, const double kappa, const double omega, const double alpha);
    void SIR(const Size number_of_nodes, const Size b, const Size l, vector<vector<vector<Size>>>& graph, vector<vector<Size>>& graph_of_graph, const double kappa, const double alpha, const Size subpopulation, const double p, const Size uid);
    bool isLinked(Size x, Size y, vector<vector<Size>>& graph); 
    void robinhood(vector<double>& P, vector<Size>& Y, Size N, double gamma); 
    void makeERGraph(const Size number_of_nodes, const Size number_of_edges, vector<vector<Size>>& graph);
    void makeSFGraph(vector<double>& P, vector<Size>& Y, Size number_of_nodes, Size number_of_edges, vector<vector<Size>>& graph);
    void makeTreeGraph(vector<vector<Size>>& graph_of_graph, const Size b, const Size l) ;
};


SIR_Metapopulation::SIR_Metapopulation(const unsigned int _number_of_nodes_, const unsigned int _number_of_edges_, const double _kappa_, const double _p_, const unsigned int _number_of_ensembles_, const unsigned int _uid_, const double _alpha_, const Size _b_, const Size _l_) 
    : number_of_nodes(_number_of_nodes_), number_of_edges(_number_of_edges_), kappa(_kappa_), p(_p_), number_of_ensembles(_number_of_ensembles_), uid(_uid_), alpha(_alpha_), b(_b_), l(_l_) {
}

SIR_Metapopulation::~SIR_Metapopulation() {
}

void SIR_Metapopulation::Run() {
    Size subpopulation = number_of_nodes/pow(b, l);
    Size subedges = number_of_edges/pow(b, l);
    vector<vector<vector<Size>>> graph(pow(b, l), vector<vector<Size>>(subpopulation));
    vector<vector<Size>> graph_of_graph((pow(b, l+1)-1)/(b-1));

    vector<double> P(subpopulation, 0);
    vector<Size> Y(subpopulation, 0);
    const double gamma = 2.5;
    robinhood(P, Y, subpopulation, gamma);
    for(Size n=0; n<pow(b, l); ++n) {
        makeSFGraph(P, Y, subpopulation, subedges, graph[n]);
        // makeERGraph(subpopulation, subedges, graph[n]);
    }
    

    makeTreeGraph(graph_of_graph, b, l);
    for(Size ensemble=0; ensemble<number_of_ensembles; ++ensemble) { // loop over ensemble
        SIR(number_of_nodes, b, l, graph, graph_of_graph, kappa, alpha, subpopulation, p, uid);
    }
}

Size SIR_Metapopulation::find_subgraph_given_distance(vector<vector<Size>>& graph_of_graph, Size picked_subgraph_1, 
    Size distance, const Size b, const Size l) {
    
    vector<Size> visited(graph_of_graph.size(), 0); // default values of vector are 0 
    queue<Size> bfs_queue;
    vector<Size> distance_from_seed(graph_of_graph.size(), 0);
    Size order_in_distance, picked_subgraph_2 = 0;
    bfs_queue.push(picked_subgraph_1);
    visited[picked_subgraph_1] = true;
    distance_from_seed[picked_subgraph_1] = 0;
    while(!bfs_queue.empty()) {
        Size u = bfs_queue.front();
        bfs_queue.pop();
        for(auto j : graph_of_graph[u]){
            if(!visited[j]){
                visited[j] = true;
                bfs_queue.push(j);
                distance_from_seed[j] = distance_from_seed[u]+1;
            }
            else if (graph_of_graph[j].size() == 1 && j!=picked_subgraph_1){
                distance_from_seed[j] = distance_from_seed[u]+1;
            }
        }
    }
    if(distance==1) {
        RandomIntGenerator rndInt_distance(0, b-2);
        order_in_distance=rndInt_distance();
    }
    else {
        RandomIntGenerator rndInt_distance(0, pow(b, distance-1)*(b-1)-1);
        order_in_distance = rndInt_distance();
    }
    Size count=0;
    for(Size i=(pow(b, l)-1)/(b-1); i<graph_of_graph.size(); ++i) {
        if(distance_from_seed[i]==2*distance && i!= picked_subgraph_1) {
            if(count==order_in_distance) {
                picked_subgraph_2=i-((pow(b, l)-1)/(b-1));
                break;
            }
            ++count;
        }
    }
    return picked_subgraph_2;
}

void SIR_Metapopulation::switch_node(Size picked_site_1, Size distance, Size picked_node_in_subgraph_2,
    vector<vector<vector<Size>>>& graph, vector<vector<Size>>& state, vector<vector<Size>>& graph_of_graph, 
    const Size b, const Size l, vector<vector<double>>& deltaP,
    vector<Size>& SI, vector<Size>& SE, vector<Size>& E, vector<Size>& I, vector<Size>& R,
    const double kappa, const double omega, const double alpha) {
    Size picked_subgraph_1 = picked_site_1%((Size)pow(b,l));
    Size picked_node_in_subgraph_1 = picked_site_1/(pow(b,l));
    Size picked_subgraph_2 = find_subgraph_given_distance(graph_of_graph, picked_subgraph_1+(pow(b, l)-1)/(b-1), distance, b, l);
    
    if(state[picked_subgraph_1][picked_node_in_subgraph_1] == state[picked_subgraph_2][picked_node_in_subgraph_2]) {
        return;
    }
    if(state[picked_subgraph_1][picked_node_in_subgraph_1] == 0) {
        if(state[picked_subgraph_2][picked_node_in_subgraph_2] == 1) {
            state[picked_subgraph_1][picked_node_in_subgraph_1] = 1; // S -> E
            deltaP[picked_subgraph_1][picked_node_in_subgraph_1] = alpha;
            ++E[picked_subgraph_1];
            for(auto nn : graph[picked_subgraph_1][picked_node_in_subgraph_1]){ 
                if( state[picked_subgraph_1][nn]==0 ) { // if neareast neighbor is S, 
                    deltaP[picked_subgraph_1][nn] += omega;
                    ++SE[picked_subgraph_1];
                }
                else if ( state[picked_subgraph_1][nn]==1 ) { // if neareast neighbor is E, 
                    --SE[picked_subgraph_1];
                }
                else if( state[picked_subgraph_1][nn]==2 ) { // if neareast neighbor is I, 
                    --SI[picked_subgraph_1];
                }
            }
            state[picked_subgraph_2][picked_node_in_subgraph_2] = 0; // E -> S
            --E[picked_subgraph_2];
            deltaP[picked_subgraph_2][picked_node_in_subgraph_2] = 0;
            for(auto nn : graph[picked_subgraph_2][picked_node_in_subgraph_2]){ 
                if( state[picked_subgraph_2][nn]==0 ) { // if neareast neighbor is S, 
                    deltaP[picked_subgraph_2][nn] -= omega;
                    --SE[picked_subgraph_2];
                }
                else if ( state[picked_subgraph_2][nn]==1 ) { // if neareast neighbor is E, 
                    ++SE[picked_subgraph_2];
                    deltaP[picked_subgraph_2][picked_node_in_subgraph_2] += omega;
                }
                else if( state[picked_subgraph_2][nn]==2 ) { // if neareast neighbor is I, 
                    ++SI[picked_subgraph_2];
                    deltaP[picked_subgraph_2][picked_node_in_subgraph_2] += kappa;
                }
            }
        }
        else if(state[picked_subgraph_2][picked_node_in_subgraph_2] == 2) {
            state[picked_subgraph_1][picked_node_in_subgraph_1] = 2; // S -> I
            deltaP[picked_subgraph_1][picked_node_in_subgraph_1] = 1;
            ++I[picked_subgraph_1];
            for(auto nn : graph[picked_subgraph_1][picked_node_in_subgraph_1]){
                if( state[picked_subgraph_1][nn]==0 ) { // if neareast neighbor is S, 
                    deltaP[picked_subgraph_1][nn] += kappa;
                    ++SI[picked_subgraph_1];
                }
                else if ( state[picked_subgraph_1][nn]==1 ) { // if neareast neighbor is E, 
                    --SE[picked_subgraph_1];
                }
                else if( state[picked_subgraph_1][nn]==2 ) { // if neareast neighbor is I, 
                    --SI[picked_subgraph_1];
                }
            }
            state[picked_subgraph_2][picked_node_in_subgraph_2] = 0; // I -> S
            deltaP[picked_subgraph_2][picked_node_in_subgraph_2] = 0;
            --I[picked_subgraph_2];
            for(auto nn : graph[picked_subgraph_2][picked_node_in_subgraph_2]){ 
                if( state[picked_subgraph_2][nn]==0 ) { // if neareast neighbor is S, 
                    deltaP[picked_subgraph_2][nn] -= kappa;
                    --SI[picked_subgraph_2];
                }
                else if ( state[picked_subgraph_2][nn]==1 ) { // if neareast neighbor is E, 
                    ++SE[picked_subgraph_2];
                    deltaP[picked_subgraph_2][picked_node_in_subgraph_2] += omega;
                }
                else if( state[picked_subgraph_2][nn]==2 ) { // if neareast neighbor is I, 
                    ++SI[picked_subgraph_2];
                    deltaP[picked_subgraph_2][picked_node_in_subgraph_2] += kappa;
                }
            }
        }
        else if(state[picked_subgraph_2][picked_node_in_subgraph_2] == 3) {
            state[picked_subgraph_1][picked_node_in_subgraph_1] = 3; // S -> R
            deltaP[picked_subgraph_1][picked_node_in_subgraph_1] = 0;
            ++R[picked_subgraph_1];
            for(auto nn : graph[picked_subgraph_1][picked_node_in_subgraph_1]){
                if ( state[picked_subgraph_1][nn]==1 ) { // if neareast neighbor is E, 
                    --SE[picked_subgraph_1];
                }
                else if( state[picked_subgraph_1][nn]==2 ) { // if neareast neighbor is I, 
                    --SI[picked_subgraph_1];
                }
            }
            state[picked_subgraph_2][picked_node_in_subgraph_2] = 0; // R -> S
            deltaP[picked_subgraph_2][picked_node_in_subgraph_2] = 0;
            --R[picked_subgraph_2];
            for(auto nn : graph[picked_subgraph_2][picked_node_in_subgraph_2]){ 
                if ( state[picked_subgraph_2][nn]==1 ) { // if neareast neighbor is E, 
                    ++SE[picked_subgraph_2];
                    deltaP[picked_subgraph_2][picked_node_in_subgraph_2] += omega;
                }
                else if( state[picked_subgraph_2][nn]==2 ) { // if neareast neighbor is I, 
                    ++SI[picked_subgraph_2];
                    deltaP[picked_subgraph_2][picked_node_in_subgraph_2] += kappa;
                }
            }
        }
    }
    else if(state[picked_subgraph_1][picked_node_in_subgraph_1] == 1) {
        if(state[picked_subgraph_2][picked_node_in_subgraph_2] == 0) {
            state[picked_subgraph_1][picked_node_in_subgraph_1] = 0; // E -> S
            --E[picked_subgraph_1];
            deltaP[picked_subgraph_1][picked_node_in_subgraph_1] = 0;
            for(auto nn : graph[picked_subgraph_1][picked_node_in_subgraph_1]){ 
                if( state[picked_subgraph_1][nn]==0 ) { // if neareast neighbor is S, 
                    deltaP[picked_subgraph_1][nn] -= omega;
                    --SE[picked_subgraph_1];
                }
                else if ( state[picked_subgraph_1][nn]==1 ) { // if neareast neighbor is E, 
                    ++SE[picked_subgraph_1];
                    deltaP[picked_subgraph_1][picked_node_in_subgraph_1] += omega;
                }
                else if( state[picked_subgraph_1][nn]==2 ) { // if neareast neighbor is I, 
                    ++SI[picked_subgraph_1];
                    deltaP[picked_subgraph_1][picked_node_in_subgraph_1] += kappa;
                }
            }
            state[picked_subgraph_2][picked_node_in_subgraph_2] = 1; // S -> E
            deltaP[picked_subgraph_2][picked_node_in_subgraph_2] = alpha;
            ++E[picked_subgraph_2];
            for(auto nn : graph[picked_subgraph_2][picked_node_in_subgraph_2]){ 
                if( state[picked_subgraph_2][nn]==0 ) { // if neareast neighbor is S, 
                    deltaP[picked_subgraph_2][nn] += omega;
                    ++SE[picked_subgraph_2];
                }
                else if ( state[picked_subgraph_2][nn]==1 ) { // if neareast neighbor is E, 
                    --SE[picked_subgraph_2];
                }
                else if( state[picked_subgraph_2][nn]==2 ) { // if neareast neighbor is I, 
                    --SI[picked_subgraph_2];
                }
            }
        }
        else if(state[picked_subgraph_2][picked_node_in_subgraph_2] == 2) {
            state[picked_subgraph_1][picked_node_in_subgraph_1] = 2; // E -> I
            deltaP[picked_subgraph_1][picked_node_in_subgraph_1] = 1;
            ++I[picked_subgraph_1];
            --E[picked_subgraph_1];
            for(auto nn : graph[picked_subgraph_1][picked_node_in_subgraph_1]){
                if( state[picked_subgraph_1][nn]==0 ) { // if neareast neighbor is S, 
                    deltaP[picked_subgraph_1][nn] += kappa;
                    deltaP[picked_subgraph_1][nn] -= omega;
                    --SE[picked_subgraph_1];
                    ++SI[picked_subgraph_1];
                }
            }
            state[picked_subgraph_2][picked_node_in_subgraph_2] = 1; // I -> E
            deltaP[picked_subgraph_2][picked_node_in_subgraph_2] = alpha;
            --I[picked_subgraph_2];
            ++E[picked_subgraph_2];
            for(auto nn : graph[picked_subgraph_2][picked_node_in_subgraph_2]){ 
                if( state[picked_subgraph_2][nn]==0 ) { // if neareast neighbor is S, 
                    deltaP[picked_subgraph_2][nn] -= kappa;
                    deltaP[picked_subgraph_2][nn] += omega;
                    ++SE[picked_subgraph_2];
                    --SI[picked_subgraph_2];
                }
            }
        }
        else if(state[picked_subgraph_2][picked_node_in_subgraph_2] == 3) {
            state[picked_subgraph_1][picked_node_in_subgraph_1] = 3; // E -> R
            deltaP[picked_subgraph_1][picked_node_in_subgraph_1] = 0;
            --E[picked_subgraph_1];
            ++R[picked_subgraph_1];
            for(auto nn : graph[picked_subgraph_1][picked_node_in_subgraph_1]){
                if( state[picked_subgraph_1][nn]==0 ) { // if neareast neighbor is S, 
                    deltaP[picked_subgraph_1][nn] -= omega;
                    --SE[picked_subgraph_1];
                }
            }
            state[picked_subgraph_2][picked_node_in_subgraph_2] = 1; // R -> E
            deltaP[picked_subgraph_2][picked_node_in_subgraph_2] = alpha;
            --R[picked_subgraph_2];
            ++E[picked_subgraph_2];
            for(auto nn : graph[picked_subgraph_2][picked_node_in_subgraph_2]){ 
                if( state[picked_subgraph_2][nn]==0 ) { // if neareast neighbor is S, 
                    deltaP[picked_subgraph_2][nn] += omega;
                    ++SE[picked_subgraph_2];
                }
            }
        }
    }
    else if(state[picked_subgraph_1][picked_node_in_subgraph_1] == 2) {
        if(state[picked_subgraph_2][picked_node_in_subgraph_2] == 0) {
            state[picked_subgraph_1][picked_node_in_subgraph_1] = 0; // I -> S
            --I[picked_subgraph_1];
            deltaP[picked_subgraph_1][picked_node_in_subgraph_1] = 0;
            for(auto nn : graph[picked_subgraph_1][picked_node_in_subgraph_1]){ 
                if( state[picked_subgraph_1][nn]==0 ) { // if neareast neighbor is S, 
                    deltaP[picked_subgraph_1][nn] -= kappa;
                    --SI[picked_subgraph_1];
                }
                else if ( state[picked_subgraph_1][nn]==1 ) { // if neareast neighbor is E, 
                    ++SE[picked_subgraph_1];
                    deltaP[picked_subgraph_1][picked_node_in_subgraph_1] += omega;
                }
                else if( state[picked_subgraph_1][nn]==2 ) { // if neareast neighbor is I, 
                    ++SI[picked_subgraph_1];
                    deltaP[picked_subgraph_1][picked_node_in_subgraph_1] += kappa;
                }
            }
            state[picked_subgraph_2][picked_node_in_subgraph_2] = 2; // S -> I
            deltaP[picked_subgraph_2][picked_node_in_subgraph_2] = 1;
            ++I[picked_subgraph_2];
            for(auto nn : graph[picked_subgraph_2][picked_node_in_subgraph_2]){ 
                if( state[picked_subgraph_2][nn]==0 ) { // if neareast neighbor is S, 
                    deltaP[picked_subgraph_2][nn] += kappa;
                    ++SI[picked_subgraph_2];
                }
                else if ( state[picked_subgraph_2][nn]==1 ) { // if neareast neighbor is E, 
                    --SE[picked_subgraph_2];
                }
                else if( state[picked_subgraph_2][nn]==2 ) { // if neareast neighbor is I, 
                    --SI[picked_subgraph_2];
                }
            }
        }
        else if(state[picked_subgraph_2][picked_node_in_subgraph_2] == 1) {
            state[picked_subgraph_1][picked_node_in_subgraph_1] = 1; // I -> E
            deltaP[picked_subgraph_1][picked_node_in_subgraph_1] = alpha;
            --I[picked_subgraph_1];
            ++E[picked_subgraph_1];
            for(auto nn : graph[picked_subgraph_1][picked_node_in_subgraph_1]){
                if( state[picked_subgraph_1][nn]==0 ) { // if neareast neighbor is S, 
                    deltaP[picked_subgraph_1][nn] -= kappa;
                    deltaP[picked_subgraph_1][nn] += omega;
                    ++SE[picked_subgraph_1];
                    --SI[picked_subgraph_1];
                }
            }
            state[picked_subgraph_2][picked_node_in_subgraph_2] = 2; // E -> I
            deltaP[picked_subgraph_2][picked_node_in_subgraph_2] = 1;
            --E[picked_subgraph_2];
            ++I[picked_subgraph_2];
            for(auto nn : graph[picked_subgraph_2][picked_node_in_subgraph_2]){ 
                if( state[picked_subgraph_2][nn]==0 ) { // if neareast neighbor is S, 
                    deltaP[picked_subgraph_2][nn] += kappa;
                    deltaP[picked_subgraph_2][nn] -= omega;
                    --SE[picked_subgraph_2];
                    ++SI[picked_subgraph_2];
                }
            }
        }
        else if(state[picked_subgraph_2][picked_node_in_subgraph_2] == 3) {
            state[picked_subgraph_1][picked_node_in_subgraph_1] = 3; // I -> R
            deltaP[picked_subgraph_1][picked_node_in_subgraph_1] = 0;
            --I[picked_subgraph_1];
            ++R[picked_subgraph_1];
            for(auto nn : graph[picked_subgraph_1][picked_node_in_subgraph_1]){
                if( state[picked_subgraph_1][nn]==0 ) { // if neareast neighbor is S, 
                    deltaP[picked_subgraph_1][nn] -= kappa;
                    --SI[picked_subgraph_1];
                }
            }
            state[picked_subgraph_2][picked_node_in_subgraph_2] = 2; // R -> I
            deltaP[picked_subgraph_2][picked_node_in_subgraph_2] = 1;
            --R[picked_subgraph_2];
            ++I[picked_subgraph_2];
            for(auto nn : graph[picked_subgraph_2][picked_node_in_subgraph_2]){ 
                if( state[picked_subgraph_2][nn]==0 ) { // if neareast neighbor is S, 
                    deltaP[picked_subgraph_2][nn] += kappa;
                    ++SI[picked_subgraph_2];
                }
            }
        }
    }
    else {
        if(state[picked_subgraph_2][picked_node_in_subgraph_2] == 0) {
            state[picked_subgraph_1][picked_node_in_subgraph_1] = 0; // R -> S
            --R[picked_subgraph_1];
            deltaP[picked_subgraph_1][picked_node_in_subgraph_1] = 0;
            for(auto nn : graph[picked_subgraph_1][picked_node_in_subgraph_1]){ 
                if ( state[picked_subgraph_1][nn]==1 ) { // if neareast neighbor is E, 
                    ++SE[picked_subgraph_1];
                    deltaP[picked_subgraph_1][picked_node_in_subgraph_1] += omega;
                }
                else if( state[picked_subgraph_1][nn]==2 ) { // if neareast neighbor is I, 
                    ++SI[picked_subgraph_1];
                    deltaP[picked_subgraph_1][picked_node_in_subgraph_1] += kappa;
                }
            }
            state[picked_subgraph_2][picked_node_in_subgraph_2] = 3; // S -> R
            deltaP[picked_subgraph_2][picked_node_in_subgraph_2] = 0;
            ++R[picked_subgraph_2];
            for(auto nn : graph[picked_subgraph_2][picked_node_in_subgraph_2]){ 
                if ( state[picked_subgraph_2][nn]==1 ) { // if neareast neighbor is E, 
                    --SE[picked_subgraph_2];
                }
                else if( state[picked_subgraph_2][nn]==2 ) { // if neareast neighbor is I, 
                    --SI[picked_subgraph_2];
                }
            }
        }
        else if(state[picked_subgraph_2][picked_node_in_subgraph_2] == 1) {
            state[picked_subgraph_1][picked_node_in_subgraph_1] = 1; // R -> E
            deltaP[picked_subgraph_1][picked_node_in_subgraph_1] = alpha;
            --R[picked_subgraph_1];
            ++E[picked_subgraph_1];
            for(auto nn : graph[picked_subgraph_1][picked_node_in_subgraph_1]){
                if( state[picked_subgraph_1][nn]==0 ) { // if neareast neighbor is S, 
                    deltaP[picked_subgraph_1][nn] += omega;
                    ++SE[picked_subgraph_1];
                }
            }
            state[picked_subgraph_2][picked_node_in_subgraph_2] = 3; // E -> R
            deltaP[picked_subgraph_2][picked_node_in_subgraph_2] = 0;
            --E[picked_subgraph_2];
            ++R[picked_subgraph_2];
            for(auto nn : graph[picked_subgraph_2][picked_node_in_subgraph_2]){ 
                if( state[picked_subgraph_2][nn]==0 ) { // if neareast neighbor is S, 
                    deltaP[picked_subgraph_2][nn] -= omega;
                    --SE[picked_subgraph_2];
                }
            }
        }
        else if(state[picked_subgraph_2][picked_node_in_subgraph_2] == 2) {
            state[picked_subgraph_1][picked_node_in_subgraph_1] = 2; // R -> I
            deltaP[picked_subgraph_1][picked_node_in_subgraph_1] = 1;
            --R[picked_subgraph_1];
            ++I[picked_subgraph_1];
            for(auto nn : graph[picked_subgraph_1][picked_node_in_subgraph_1]){ 
                if( state[picked_subgraph_1][nn]==0 ) { // if neareast neighbor is S, 
                    deltaP[picked_subgraph_1][nn] += kappa;
                    ++SI[picked_subgraph_1];
                }
            }
            state[picked_subgraph_2][picked_node_in_subgraph_2] = 3; // I -> R
            deltaP[picked_subgraph_2][picked_node_in_subgraph_2] = 0;
            --I[picked_subgraph_2];
            ++R[picked_subgraph_2];
            for(auto nn : graph[picked_subgraph_2][picked_node_in_subgraph_2]){
                if( state[picked_subgraph_2][nn]==0 ) { // if neareast neighbor is S, 
                    deltaP[picked_subgraph_2][nn] -= kappa;
                    --SI[picked_subgraph_2];
                }
            }
        }
    }
}

void SIR_Metapopulation::SIR(const Size number_of_nodes, const Size b, const Size l,
    vector<vector<vector<Size>>>& graph, vector<vector<Size>>& graph_of_graph, 
    const double kappa, const double alpha, const Size subpopulation, const double p, const Size uid) {
    
    RandomRealGenerator rnd(0.0, 1.0);
    RandomIntGenerator rndInt(0, number_of_nodes-1);
    double xi = 0.35;
    Size total_SI = 0, total_SE = 0, total_E = 0, total_I = 0, total_R = 0, picked_site;
    Size picked_site_1, picked_site_2, picked_distance, picked_subgraph, picked_node_in_subgraph;
    double t = 0, dt, r_1, r_2, r_subgraph, omega = 0;
    vector<vector<Size>> state(pow(b, l), vector<Size>(subpopulation, 0));
    vector<vector<double>> deltaP(pow(b, l), vector<double>(subpopulation, 0));
    vector<vector<Size>> order_of_site(pow(b, l), vector<Size>(number_of_nodes));
    vector<Size> SI(pow(b, l), 0), SE(pow(b, l), 0), E(pow(b, l), 0), I(pow(b, l), 0), R(pow(b, l), 0);
    vector<double> deltaP_subgraph(l);
    for(Size i=0; i<l; ++i) {
        deltaP_subgraph[i] = pow(b, i)*(b-1)*exp(-((double)i-1.)/xi);
    }
    partial_sum(deltaP_subgraph.begin(), deltaP_subgraph.end(), deltaP_subgraph.begin());

    // Single seed initial condition : I
    state[0][0] = 2;
    deltaP[0][0] = 1;
    ++I[0];
    for(auto nn : graph[0][0]) {
        if( state[0][nn]==0 ) {
            deltaP[0][nn] += kappa;
            ++SI[0];
        }
    }

    // pow(b, l) seed initial condition : I
    // for(Size n=0; n<pow(b, l); ++n) {
    //     state[n][0] = 2;
    //     deltaP[n][0] = 1;
    //     ++I[n];
    //     for(auto nn : graph[n][0]) {
    //         if( state[n][nn]==0 ) {
    //             deltaP[n][nn] += kappa;
    //             ++SI[n];
    //         }
    //     }
    // }


    // Single seed initial condition : E
    // state[0][0] = 1;
	// deltaP[0][0] = alpha;
	// ++E[0];
	// for(auto nn : graph[0][0]) {
	// 	if( state[0][nn]==0 ) {
	// 		deltaP[0][nn] += omega;
	// 		++SE[0];
	// 	}
	// }

    // pow(b, l) seed initial condition : E
    // for(Size n=0; n<pow(b, l); ++n) {
    //     state[n][0] = 1;
    //     deltaP[n][0] = alpha;
    //     ++E[n];
    //     for(auto nn : graph[n][0]) {
    //         if( state[n][nn]==0 ) {
    //             deltaP[n][nn] += omega;
    //             ++SE[n];
    //         }
    //     }
    // }


    ostringstream oss;
    oss << "SIR_" << number_of_nodes << "_" << kappa << "_"
    << p << "_" << uid << ".table";
    ofstream os(oss.str());
    double previous_total_I;
    while( true ) {
        for(Size n=0; n<pow(b, l); ++n) {
            total_SI += SI[n]; total_SE += SE[n]; total_E += E[n]; total_I += I[n]; total_R += R[n];
        }
        os << t << '\t' << (double)(number_of_nodes-total_E-total_I-total_R)/number_of_nodes
        << '\t' << (double)total_E/number_of_nodes 
        << '\t' << (double)total_I/number_of_nodes 
        << '\t' << (double)total_R/number_of_nodes << endl;
        if((total_I==0 && total_E==0)) break;
        r_1 = rnd();
		dt = -log(r_1)/(double)(total_SI*kappa+total_I*1.0);
        t += dt;
        if(rnd() < 1.0/total_I) {
            for(Size k=0; k<pow(b,l); ++k) {
                r_subgraph = rnd();
                picked_distance = distance(deltaP_subgraph.begin(), lower_bound(deltaP_subgraph.begin(), deltaP_subgraph.end(), deltaP_subgraph.back()*r_subgraph));
                picked_site_1 = rndInt();
                picked_site_2 = rndInt();
                switch_node(picked_site_1, picked_distance+1, picked_site_2/pow(b, l), 
                    graph, state, graph_of_graph, b, l, deltaP, SI, SE, E, I, R, kappa, omega, alpha);
            }
        }
        previous_total_I = total_I;
        total_SI = 0; total_SE = 0; total_E = 0; total_I = 0; total_R = 0;
        vector<double> deltaP_total;
        for(Size n=0; n<pow(b, l); ++n) {
            deltaP_total.insert(deltaP_total.end(), deltaP[n].begin(), deltaP[n].end());
        }
        r_2 = rnd();
        partial_sum(deltaP_total.begin(), deltaP_total.end(), deltaP_total.begin());
        picked_site = distance(deltaP_total.begin(), lower_bound(deltaP_total.begin(), deltaP_total.end(), deltaP_total.back()*r_2));
        picked_subgraph = picked_site/subpopulation;
        picked_node_in_subgraph = picked_site%subpopulation;
        if(state[picked_subgraph][picked_node_in_subgraph]==0) { // S -> I takes place
            state[picked_subgraph][picked_node_in_subgraph] = 2; // S -> I
            deltaP[picked_subgraph][picked_node_in_subgraph] = 1;
            ++I[picked_subgraph];
            for(auto nn : graph[picked_subgraph][picked_node_in_subgraph]){ 
                if( state[picked_subgraph][nn]==0 ) { // if neareast neighbor is S, 
                    deltaP[picked_subgraph][nn] += kappa;
                    ++SI[picked_subgraph];
                }
                else if ( state[picked_subgraph][nn]==2 ) { // if neareast neighbor is I, 
                    --SI[picked_subgraph];
                }
            }
        }
        else if(state[picked_subgraph][picked_node_in_subgraph]==2) { // I -> R takes place 
            state[picked_subgraph][picked_node_in_subgraph] = 3; // I -> R 
            deltaP[picked_subgraph][picked_node_in_subgraph] = 0;
            --I[picked_subgraph];
            ++R[picked_subgraph];
            for(auto nn : graph[picked_subgraph][picked_node_in_subgraph]){
                if( state[picked_subgraph][nn]==0 ) { // if neareast neighbor is S, 
                    deltaP[picked_subgraph][nn] -= kappa;
                    --SI[picked_subgraph];
                }
            }
        }
	}
}

bool SIR_Metapopulation::isLinked(Size x, Size y, vector<vector<Size>>& graph){ 
	return ((find(graph[x].begin(),graph[x].end(),y)) != graph[x].end());
}


void SIR_Metapopulation::robinhood(vector<double>& P, vector<Size>& Y, Size N, double gamma) { // builds robin hood table. used for dynamical scale-free network
    int n_poor, n_rich, n1, n2;
    double psum;
    double mu = 1.0/(gamma-1);

    vector<double> p(N);
    vector<Size> rich(N);
    vector<Size> poor(N);
    
    psum = 0;
    for(Size i=0; i<N; i++) {
        p[i] = pow(i+1,-mu); 
        psum += p[i];
    }
    for(Size i=0; i<N; i++) {
        p[i] /= psum;
    }
    
    n_poor=n_rich=0;
    for(Size i=0; i<N; i++) {
        if(p[i]>1/(double)N) rich[n_rich++] = i;
        else poor[n_poor++] = i;
    }

    for(Size i=0; i<N; i++) {
        if(n_rich==0) break;
        
        n1 = poor[n_poor-1]; // (3)
        n2 = rich[n_rich-1]; // (3)
        
        P[n1] = N*p[n1];
        Y[n1] = n2;
        p[n2] -= 1/(double)N - p[n1]; // (6)
        
        if(p[n2]>1/(double)N) //still rich
            n_poor--;
        else{ // becomes new poor
            poor[n_poor-1]=rich[n_rich-1];
            n_rich--;
        }
    }
    if(n_poor>0) {
        for(int i=0; i<n_poor; i++){
            n1 = poor[i];
            P[n1] = N*p[n1];
            Y[n1] = n1;
        }
    }
}

void SIR_Metapopulation::makeERGraph(const Size number_of_nodes, const Size number_of_edges, vector<vector<Size>>& graph) {
	RandomUnsignedIntGenerator rnd(0, number_of_nodes - 1);
	for(Size i = 0; i < number_of_edges; ++i) {
		Size n1 = rnd();
		Size n2 = rnd();
		while(n1 == n2 || isLinked(n1, n2, graph)) {
			n1 = rnd();
			n2 = rnd();
		}
		graph[n1].emplace_back(n2);
		graph[n2].emplace_back(n1);
	}
}

void SIR_Metapopulation::makeSFGraph(vector<double>& P, vector<Size>& Y, Size number_of_nodes, Size number_of_edges, vector<vector<Size>>& graph) {
    RandomIntGenerator rndInt(0, number_of_nodes-1);
    RandomRealGenerator rnd(0.0, 1.0);
    int n1, n2, K;
    for(Size i=0; i<number_of_edges; ++i){
        do {
            K = rndInt();
            if(rnd() < P[K]) n1 = K;
            else n1 = Y[K];
            K = rndInt();
            if(rnd() < P[K]) n2 = K;
            else n2 = Y[K];
        } while(n1==n2 || isLinked(n1, n2, graph));
        graph[n1].emplace_back(n2);
        graph[n2].emplace_back(n1);
    }
}

void SIR_Metapopulation::makeTreeGraph(vector<vector<Size>>& graph_of_graph, const Size b, const Size l) {
    for(Size i=0; i<l; ++i) {
        for(Size j=(pow(b, i)-1)/(b-1); j<(pow(b, i+1)-1)/(b-1); ++j) {
            for(Size k=1; k<=b; ++k) {
                graph_of_graph[j].emplace_back(j*b+k);
                graph_of_graph[j*b+k].emplace_back(j);
            }
        }
    }
}

 
int main(int argc, char* argv[]) {
	const unsigned int number_of_nodes = stoul(argv[1]);
	const unsigned int number_of_edges = stoul(argv[2]);
    const double       kappa = stof(argv[3]);
    const double       p = stof(argv[4]);
	const unsigned int number_of_ensembles = stoul(argv[5]);
    const unsigned int uid = stoul(argv[6]);
    const double       alpha = 0;
    const Size         b = 2;
    const Size         l = 3;
        
    std::ios_base::sync_with_stdio(false);
    cin.tie(nullptr); cout.tie(nullptr);

    SIR_Metapopulation* Model = new SIR_Metapopulation(number_of_nodes, number_of_edges, kappa, p, number_of_ensembles, uid, alpha, b, l);
    
    Model -> Run();
    delete Model;
}