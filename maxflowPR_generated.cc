#include"maxflowPR.h"

void preflow(graph& g , int s , int* height , int* excess_flow , 
  bool* active , int* flow , int* capacity)
{
  height[s] = g.num_nodes();
  #pragma omp parallel for
  for (int edge = g.indexofNodes[s]; edge < g.indexofNodes[s+1]; edge ++) 
  {int v = g.edgeList[edge] ;
    int e = edge;
    flow[e] = capacity[e];
    excess_flow[v] = excess_flow[v] + flow[e];
    int rev = g.getEdge(v,s).id;
    flow[rev] = flow[rev] - flow[e];
    active[v] = true;
  }

}
void globalRelabel(graph& g , int s , int t , int* height , 
  int* capacity , int* flow)
{
  #pragma omp parallel for
  for (int v = 0; v < g.num_nodes(); v ++) 
  {
    height[v] = g.num_nodes();
  }
  height[t] = 0;
  int qsize = 1;
  int* inqueue=new int[g.num_nodes()];
  int* discovered=new int[g.num_nodes()];
  #pragma omp parallel for
  for (int t = 0; t < g.num_nodes(); t ++) 
  {
    inqueue[t] = 0;
    discovered[t] = 0;
  }
  inqueue[t] = 1;
  while (qsize ){
    #pragma omp parallel for
    for (int v = 0; v < g.num_nodes(); v ++) 
    {
      if (inqueue[v] )
        {
        for (int edge = g.indexofNodes[v]; edge < g.indexofNodes[v+1]; edge ++) 
        {int u = g.edgeList[edge] ;
          int e = edge;
          int rev = g.getEdge(u,v).id;
          if (capacity[rev] - flow[rev] > 0 && height[u] == g.num_nodes() && u != s )
            {
            discovered[u] = 1;
            if (u != t )
              {
              height[u] = height[v] + 1;
            }
          }
        }
      }
    }
    qsize = 0;
    #pragma omp parallel for reduction(+ : qsize)
    for (int v = 0; v < g.num_nodes(); v ++) 
    {
      inqueue[v] = discovered[v];
      discovered[v] = 0;
      if (inqueue[v] )
        {
        qsize = qsize + 1;
      }
    }
  }

}
void push_flow(graph& g , int v , int t , int* height , 
  int* capacity , int* flow , int* excess_flow , int* added_excess , 
  int* discovered)
{
  bool stop = false;
  for (int edge = g.indexofNodes[v]; edge < g.indexofNodes[v+1]; edge ++) 
  {int u = g.edgeList[edge] ;
    if (stop == false )
      {
      int e = edge;
      int rev = g.getEdge(u,v).id;
      if ((height[v] == height[u] + 1) && (capacity[e] - flow[e] > 0) )
        {
        int f = min(capacity[e] - flow[e],excess_flow[v]);
        excess_flow[v] = excess_flow[v] - f;
        added_excess[u] = added_excess[u] + f;
        flow[e] = flow[e] + f;
        flow[rev] = flow[rev] - f;
        if (u != t )
          {
          discovered[u] = true;
        }
        if (excess_flow[v] == 0 )
          {
          stop = true;
        }
      }
    }
  }

}
auto cal_new_height(graph& g , int v , int* capacity , int* flow , 
  int* height)
{
  int n = g.num_nodes() - 1;
  for (int edge = g.indexofNodes[v]; edge < g.indexofNodes[v+1]; edge ++) 
  {int u = g.edgeList[edge] ;
    int e = edge;
    if (capacity[e] != flow[e] )
      {
      n = min(n,height[u]);
    }
  }
  return n + 1;

}
void relabel(graph& g , int v , int* excess_flow , int* height , 
  int* new_height , int* discovered , int* capacity , int* flow
)
{
  if (excess_flow[v] > 0 || height[v] == g.num_nodes() )
    {
    new_height[v] = cal_new_height(g,v,capacity,flow,height);
    if (new_height[v] != g.num_nodes() && new_height[v] != height[v] )
      {
      discovered[v] = true;
    }
  }
  else
  {
    new_height[v] = height[v];
  }

}
auto maxflow(graph& g , int s , int t , int* weight
)
{
  int count = 0;
  int* height=new int[g.num_nodes()];
  int* new_height=new int[g.num_nodes()];
  int* excess_flow=new int[g.num_nodes()];
  bool* active=new bool[g.num_nodes()];
  int* discovered=new int[g.num_nodes()];
  int* added_excess=new int[g.num_nodes()];
  int* capacity=new int[g.num_edges()];
  int* flow=new int[g.num_edges()];
  #pragma omp parallel for
  for (int t = 0; t < g.num_nodes(); t ++) 
  {
    height[t] = 0;
    new_height[t] = 0;
    excess_flow[t] = 0;
    active[t] = false;
    discovered[t] = 0;
    added_excess[t] = 0;
  }
  #pragma omp parallel for
  for (int e = 0; e < g.num_edges(); e ++) 
  {
    capacity[e] = weight[e];
    flow[e] = 0;
  }
  preflow(g,s,height,excess_flow,active,flow,capacity);

  int active_count = (g.indexofNodes[s+1]-g.indexofNodes[s]);
  while (active_count ){
    if (count % 4 == 0 )
      {
      globalRelabel(g,s,t,height,capacity,flow);

    }
    count++;
    #pragma omp parallel for
    for (int v = 0; v < g.num_nodes(); v ++) 
    {
      if (active[v] )
        {
        push_flow(g,v,t,height,capacity,flow,excess_flow,added_excess,discovered);

      }
    }
    #pragma omp parallel for
    for (int v = 0; v < g.num_nodes(); v ++) 
    {
      if (active[v] )
        {
        relabel(g,v,excess_flow,height,new_height,discovered,capacity,flow);

      }
    }
    #pragma omp parallel for
    for (int v = 0; v < g.num_nodes(); v ++) 
    {
      if (active[v] )
        {
        height[v] = new_height[v];
        int num = added_excess[v];
        excess_flow[v] = excess_flow[v]+ num;
        added_excess[v] = 0;
      }
    }
    active_count = 0;
    #pragma omp parallel for reduction(+ : active_count)
    for (int v = 0; v < g.num_nodes(); v ++) 
    {
      if (height[v] < g.num_nodes() && discovered[v] )
        {
        active[v] = true;
        active_count = active_count + 1;
      }
      else
      {
        active[v] = false;
      }
      discovered[v] = false;
    }
    #pragma omp parallel for
    for (int v = 0; v < g.num_nodes(); v ++) 
    {
      if (active[v] )
        {
        excess_flow[v] = excess_flow[v]+ added_excess[v];
        added_excess[v] = 0;
      }
    }
  }
  return added_excess[t] + excess_flow[t];

}

int main(){
  char filename[50];
  for(char i = '1'; i <= '4'; i++){
    sprintf(filename, "dataset/input%c.txt", i);
    graph G(filename);
    G.parseGraph();
    std::cout<<G.num_nodes()<<" "<<G.num_edges()<<std::endl;

    int* weight = new int[G.num_edges()];
    std::map<int, std::vector<edge>>edges = G.getEdges();
    int idx = 0;
    for(auto i : edges){
      for(auto j : i.second){
        weight[idx] = j.weight;
        idx++;
      }
    }
    
    std::cout<<maxflow(G, 0, G.num_nodes() - 1, weight)<<"\n";
  }
  return 0;
}
