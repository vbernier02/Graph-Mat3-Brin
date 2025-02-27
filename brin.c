#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <float.h>
#include <malloc.h>

typedef unsigned short Shu;

struct strand {
  Shu node; 
  int next;
};

typedef struct strand strand;

struct strandgraph {
  Shu nbs;
  int nbstr;
  int *node;
  strand *nxt;
};

typedef struct strandgraph strgr;

strgr creegraphe_brin(int nbs, int *memoire){ 
  int max_arc, num_arc;
  Shu i, j, num;
  float taux, v;
  strgr g;
  struct mallinfo2 m_debut, m_fin;
  
  m_debut = mallinfo2();
  num_arc = 0;
  taux = 25.0; 
  num = nbs / 10;
  
  while (num > 1){
    num /= 5;
    taux /= 3.0;
  }
  taux /= 100.0;
  
  g.nbs = nbs;
  max_arc = (nbs * nbs);
  g.nbstr = 0;
  g.node = (int *)malloc(nbs * sizeof(int));
  g.nxt = (strand *)malloc(max_arc * sizeof(strand));
  srand(time(NULL)+ rand());
  
  for (i = 0; i < nbs; i++){
    g.node[i] = -1;
  } 
  for (i = 0; i < nbs; i++){
    for (j = 0; j < nbs; j++){
      v = (float)rand() / RAND_MAX;
      if (v < taux) {
	g.nxt[num_arc].node = j;
	g.nxt[num_arc].next = g.node[i];
	g.node[i] = num_arc;
	num_arc++;
      }
    }
  }
  g.nbstr = num_arc;
  
  m_fin = mallinfo2();
  *memoire = m_fin.uordblks - m_debut.uordblks;
  *memoire = 0;
  return g;
}

int find(int *parent, int x){
  
  if (parent[x] != x){
    parent[x] = find(parent, parent[x]);
  }
  
  return parent[x];
}

/*void uni(int *parent, int x, int y){
  
  int racine_x, racine_y;
  
  racine_x = find(parent, x);
  racine_y = find(parent, y);
  
  if (racine_x != racine_y){
    parent[racine_y] = racine_x;
  }
}*/

int composantes_connexes(strgr g){

  int i, j, arc, *parent, *composante, racine, racine_i, racine_j, nb_composante, memoire;
  struct mallinfo2 m_debut, m_fin;
  m_debut = mallinfo2();
  
  parent = malloc(g.nbs * sizeof(int));
  composante = calloc(g.nbs, sizeof(int)); 
  nb_composante = 0;
  
  for (i = 0; i < g.nbs; i++){
    parent[i] = i;
  }
  for (i = 0; i < g.nbs; i++){
    for (arc = g.node[i]; arc != -1; arc = g.nxt[arc].next){
      racine_i = find(parent, i);
      racine_j = find(parent, g.nxt[arc].node);
      if (racine_i != racine_j) {
	parent[racine_j] = racine_i;
      }
    }
  } 
  for (i = 0; i < g.nbs; i++){
    racine = find(parent, i);
    if (!composante[racine]){
      nb_composante++;
      composante[racine] = 1;
    }
  } 
  m_fin = mallinfo2();
  memoire = (m_fin.uordblks - m_debut.uordblks);
  memoire = 0;
  printf("Composantes connexes : %d\n", nb_composante);
  
  free(parent);
  free(composante);
  
  return memoire;
}

int dijkstra(strgr g, int source){

  int u, i, step, v, arc, *predecesseur, *vu, memoire;
  float *distance, poid;
  struct mallinfo2 m_debut, m_fin;
  m_debut = mallinfo2();
    
  distance = malloc(g.nbs * sizeof(float));
  predecesseur = malloc(g.nbs * sizeof(int));
  vu = calloc(g.nbs, sizeof(int));
  
  for (i = 0; i < g.nbs; i++){
    distance[i] = FLT_MAX;
    predecesseur[i] = -1;
  }
  distance[source] = 0.0;  
  for (step = 0; step < g.nbs; step++){
    u = -1;
    for (i = 0; i < g.nbs; i++){
      if (!vu[i] && (u == -1 || distance[i] < distance[u]))
	u = i;
    }
    if (u == -1 || distance[u] == FLT_MAX){
      break;
    }
    vu[u] = 1;   
    for (arc = g.node[u]; arc != -1; arc = g.nxt[arc].next){
      v = g.nxt[arc].node;
      poid = 1.0;
      if (!vu[v] && distance[u] + poid < distance[v]){
	distance[v] = distance[u] + poid;
	predecesseur[v] = u;
      }
    }
  }
  
  /*for (i = 0; i < g.nbs; i++){
    if (distance[i] == FLT_MAX){
    printf("Sommet %d : Aucun chemin trouvé\n", i);
    }
    else{
    printf("Sommet %d : %.2f (précédent : %d)\n", i, distance[i], predecesseur[i]);
    }
   }*/
  
  m_fin = mallinfo2();
  memoire = (m_fin.uordblks - m_debut.uordblks);
  memoire = 0;
  
  free(distance);
  free(predecesseur);
  free(vu);
  
  return memoire;
}

void liberer_graphe(strgr *g){
  
  if (g->node) free(g->node);
  if (g->nxt) free(g->nxt);
  g->node = NULL;
  g->nxt = NULL;
  g->nbs = g->nbstr = 0;
}

void mesure(int nbs, int depart, int tour){
  
  int i;
  clock_t start, end;
  double time, total_time;
  int memoire_cree, memoire_dijk, memoire_compo, total_memoire, total_memoire_cree;
  strgr g;
  total_memoire = 0;
  total_time = 0.0;
  
  for (i = 0; i < tour; i++){
    
    start = clock();
    g = creegraphe_brin(nbs, &memoire_cree);
    memoire_dijk = dijkstra(g, depart);
    memoire_compo = composantes_connexes(g);
    
    end = clock();
    time = (double)(end - start)/ CLOCKS_PER_SEC;
    
    total_memoire_cree += memoire_cree;
    total_time += time;
    printf("Temps : %.3f s\n", total_time);
    total_memoire += (memoire_cree + memoire_dijk + memoire_compo);
    liberer_graphe(&g);
  }

  total_memoire_cree = total_memoire_cree/tour;
  printf("Mémoire grap brin : %d octet\n", total_memoire_cree);
  printf("Temps moyen : %.3f s\n", total_time / tour);
  printf("Mémoire moyenne : %d octet\n", total_memoire / tour);
}

int main(){
  
  int nbs, depart, tour;

  nbs = 1000;
  depart = 0;
  tour = 1;
  
  mesure(nbs, depart, tour);
  
  return 0;
}
