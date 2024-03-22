#include"quicksi.h"

auto find(int u , int* parent)
{
  int par = parent[u];
  while (u != par ){
    u = par;
    par = parent[u];
  }
  return par;

}
void merge(int u , int v , int* parent , int* rnk
)
{
  u = find(u,parent);
  v = find(v,parent);
  if (rnk[u] > rnk[v] )
    {
    parent[v] = u;
  }
  else
  {
    parent[u] = v;
  }
  if (rnk[u] == rnk[v] )
    {
    int temp = rnk[v];
    rnk[v] = temp + 1;
  }

}
auto qi_seq(std::vector<int> degree , std::vector<int> edges , int n , int* parent , 
  int* rnk)
{
  std::vector<std::vector<int>> records;
  std::vector<int> par;
  int vertex = 0;
  int index = 0;
  while (vertex < n ){
    std::vector<int> temp;
    records.push_back(temp);

    par.push_back(-1);

    vertex++;
  }
  vertex = 0;
  while (vertex < n ){
    int deg = degree[vertex];
    while (deg > 0 ){
      int set_u = find(vertex,parent);
      int set_v = find(edges[index],parent);
      if (set_u != set_v )
        {
        merge(set_u,set_v,parent,rnk);

        int dst = edges[index];
        par[dst] = vertex;
      }
      else
      {
        int dst = edges[index];
        records[vertex].push_back(dst);

      }
      index++;
      deg--;
    }
    vertex++;
  }
  vertex = 0;
  while (vertex < n ){
    records[vertex].push_back(par[vertex]);

    vertex++;
  }
  return records;

}
bool res = false;
auto quicksi(std::vector<int> degree , std::vector<std::vector<int>> records , graph& g , int d , 
  int* H , int* F)
{
  if (d >= records.size() )
    {
    return true;
  }
  std::vector<int> temp;
  temp = records[d];
  int par = temp.back();
  // #pragma omp parallel for default(shared) private(H, F)
  for (int v = 0; v < g.num_nodes(); v ++) 
  {
    // std::cout<<"HEY\n";
    if (!res && F[v] == 0 && ((d == 0) || (d > 0 && (par == -1 || (par != -1 && g.check_if_nbr(H[par],v))))) )
      {
      if (g.getOutDegree(v) >= degree[d] )
        {
        bool flag = true;
        int index = 0;
        while (index < temp.size() - 1 ){
          int val = temp[index];
          if (flag )
            {
            if (H[val] != -1 && !g.check_if_nbr(v,H[val]) )
              {
              flag = false;
            }
          }
          index++;
        }
        // std::cout<<"HI\n";
        if (flag )
          {
            H[d] = v;
            F[v] = 1;
            if (quicksi(degree,records,g,d + 1,H,F) )
              {
              res = true;
            }
            F[v] = 0;
        }
      }
    }
  }
  return res;

}

int32_t main(int argc, char* argv[]){

    int v, e;
    std::cin >> v >> e;
    std::cout<<v<<" "<<e<<std::endl;
    std::vector<int>degree(v, 0);
    std::vector<int>edges;
    std::vector<std::vector<int>>all_edges(v);
    for(int i = 0; i < e; i++){
      int x, y;
      std::cin>>x>>y;
      all_edges[x].push_back(y);
      degree[x]++;
    }

    for(int i = 0; i < v; i++){
      for(auto j : all_edges[i]){
        edges.push_back(j);
      }
    }


    int* parent = (int*)malloc((v + 1) * sizeof(int));
    int* rnk = (int*)malloc((v + 1) * sizeof(int));
    
    for(int i = 0 ; i < v; i++){
      parent[i] = i;
      rnk[i] = 0;
    }
    std::vector<std::vector<int>>records = qi_seq(degree, edges, v, parent, rnk);

    std::cout<<"Seq Done\n";

    graph G("dataset/input2.txt");
    G.parseGraph();
    std::cout<<G.num_nodes()<<" "<<G.num_edges()<<std::endl;

    int* F = (int*)malloc(G.num_nodes() * sizeof(int));

    int* H = (int*)malloc(G.num_nodes() * sizeof(int));
    #pragma omp parallel for
    for(int i = 0 ; i < G.num_nodes(); i++){
      F[i] = 0;
      H[i] = -1;
    }

    std::cout<<"Start quicksi\n";
    double t1=omp_get_wtime();
    bool res = quicksi(degree, records, G, 0, H, F);
    double t2=omp_get_wtime();
    for(int i = 0 ; i < v; i++){
      std::cout<<H[i]<<" ";
    }
    std::cout<<"\n";
    if(res){
      std::cout<<"YES\n";
    }
    else{
      std::cout<<"NO\n";
    }
    std::cout<<t2 - t1<<std::endl;

    return 0;
}
