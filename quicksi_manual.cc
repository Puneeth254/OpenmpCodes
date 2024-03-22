#include <bits/stdc++.h>
#include <omp.h>
#include "../../graph.hpp"
using namespace std;

#define int long long
#define N 100000
vector<int>H(N, -1);
vector<int>F(N, 0);

class R{
    public:
    int deg;
    vector<int>edges;
    R(){
        deg = -1;
    }
};


class T{
    public:
    int v, p;
    R extra;
    T(){
        v = -1;
        p = -1;
    }
};

struct DisjointSets 
{ 
    int *parent, *rnk; 
    int n; 
  
    DisjointSets(int n) 
    { 
        this->n = n; 
        parent = new int[n+1]; 
        rnk = new int[n+1]; 

        #pragma omp parallel for
        for (int i = 0; i <= n; i++) 
        { 
            rnk[i] = 0; 
            parent[i] = i; 
        } 
    }
  

    int find(int u) 
    { 
        if (u != parent[u]) 
            parent[u] = find(parent[u]); 
        return parent[u]; 
    } 
  
    void merge(int x, int y) 
    { 
        x = find(x), y = find(y); 
  
        if (rnk[x] > rnk[y]){
            parent[y] = x; 
        }
        else{
            parent[x] = y;
        }
  
        if (rnk[x] == rnk[y]){
            rnk[y]++;
        }
    } 
}; 


vector<T*> qi_seq(graph& q){
    int par = -1;
    vector<T*>res;
    vector<T*>local[omp_get_max_threads()];

    DisjointSets ds(q.num_nodes());

    map<int32_t, vector<edge>>edges = q.getEdges();
    vector<int>parent(q.num_nodes(), -1);

    // #pragma omp parallel for
    for(auto i : edges){
        T* temp = new T();
        temp->v = i.first;
        R extra;
        if(q.getOutDegree(i.first) >= 3){
            extra.deg = q.getOutDegree(i.first);
        }

        vector<int>local_extra[omp_get_max_threads()];
        #pragma omp parallel for
        for(int j = 0; j < i.second.size(); j++){
            int u = i.second[j].source;
            int v = i.second[j].destination;

            int set_u = ds.find(u);
            int set_v = ds.find(v); 
            if (set_u != set_v)
            {
                ds.merge(set_u, set_v);
                parent[v] = u;
            }
            else{
                local_extra[omp_get_thread_num()].push_back(v);
            }
        }

        // #pragma omp parallel for
        for(int j = 0; j < omp_get_max_threads(); j++){
            extra.edges.insert(extra.edges.end(), local_extra[j].begin(), local_extra[j].begin() + local_extra[j].size());
        }

        temp->extra = extra;
        local[omp_get_thread_num()].push_back(temp);
    }

    for(int i = 0; i < omp_get_max_threads(); i++){
        res.insert(res.end(), local[i].begin(), local[i].begin() + local[i].size());
    }
    cout<<"Seq Done\n";

    // #pragma omp parallel for
    for(auto i : res){
        i->p = parent[i->v];
        cout<<i->v<<" "<<i->p<<endl;
        cout<<i->extra.deg<<"\n"<<"HI ";
        for(auto j : i->extra.edges){
            cout<<j<<" ";
        }
        cout<<endl;
    }

    return res;
}

bool quicksi(vector<T*>&seq, graph& g, int d){
    int alpha = g.num_nodes();
    if(d >= seq.size()){
        return true;
    }
    cout<<"Enter "<<d<<" "<<seq.size()<<endl;
    T* temp = seq[d];
    vector<int>vertices;
    vector<int>local[omp_get_max_threads()];
    if(d == 0){
        
        #pragma omp parallel for
        for(int i = 0; i < alpha; i++){
            if(F[i] == 0){
                local[omp_get_thread_num()].push_back(i);
            }
        }
        
        // #pragma omp parallel for
        for(int i = 0; i < omp_get_max_threads(); i++){
            vertices.insert(vertices.end(), local[i].begin(), local[i].begin() + local[i].size());
        }
    }
    else{
        cout<<temp->p<<endl;

        #pragma omp parallel for
        for(int i = 0; i < alpha; i++){
            if(F[i] == 0 && (temp->p == -1 || (temp->p != -1 && g.check_if_nbr(H[temp->p], i)))){
                local[omp_get_thread_num()].push_back(i);
            }
        }

        // #pragma omp parallel for
        for(int i = 0; i < omp_get_max_threads(); i++){
            vertices.insert(vertices.end(), local[i].begin(), local[i].begin() + local[i].size());
        }
    }
    cout<<d<<" "<<vertices.size()<<endl;

    bool res = false;

    // #pragma omp parallel for
    for(auto i : vertices){
        if(res){
            continue;
        }
        cout<<"NEXT ONE "<<i<<endl;
        R extra =temp->extra;
        if(extra.deg != -1){
            if(g.getOutDegree(i) < extra.deg){
                cout<<"DEG FAIL\n";
                continue;
            }
        }
        bool flag = true;

        #pragma omp parallel for
        for(auto j : extra.edges){
            if(!flag){
                continue;
            }
            if(H[j] != -1 && !(g.check_if_nbr(i, H[j]))){
                cout<<i<<" "<<H[j]<<endl;
                flag = false;
            }
        }

        if(!flag){
            cout<<"EDGE FAIL\n";
            continue;
        }

        H[d] = i;
        F[i] = 1;
        cout<<"NEXT LEVEL\n";
        if(quicksi(seq, g, d + 1) == true){
            res =  true;
        }
        F[i] = 0;
    }
    return res;
}

int32_t main(){
    graph Q("dataset/input1.txt");
    Q.parseGraph();
    cout<<Q.num_nodes()<<" "<<Q.num_edges()<<endl;

    graph G("dataset/input2.txt");
    G.parseGraph();
    cout<<G.num_nodes()<<" "<<G.num_edges()<<endl;


    vector<T*>seq_q = qi_seq(Q);

    
    bool res = quicksi(seq_q, G, 0);
    if(res){
        cout<<"YES\n";
    }
    else{
        cout<<"NO\n";
    }
    return 0;
}
