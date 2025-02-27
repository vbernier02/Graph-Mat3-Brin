#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <float.h>
#include <malloc.h>

typedef struct triplet{
  int i;
  int j;
  float poids;
} triple_t;

typedef struct gramaco{
  int nbs;
  int nba;
  triple_t *ares;
} gramaco;

gramaco creegraphe_mat3(int nbs, int *memoire){
  int max, nba, num, i, j;
  float taux, v;
  gramaco g;
  struct mallinfo2 m_debut, m_fin;
  
  m_debut = mallinfo2();
  max = nbs;
  nba = 0;
  num = nbs / 10;
  taux = 3.0;
  
  while (num > 1){
    num /= 5;
    taux /= 3.0;
  }
  taux /= 100.0;
  g.nbs = nbs;
  g.ares = (triple_t *)malloc(max * sizeof(triple_t));
  srand(time(NULL) + rand());
  
  for (i = 0; i < nbs; i++){
    for (j = 0; j < nbs; j++){
	v = (float)rand() / RAND_MAX;
	if (v < taux){
	  if (nba >= max){
	    max *= 2;
	    g.ares = (triple_t *)realloc(g.ares, max * sizeof(triple_t));
	  }
	  g.ares[nba].i = i;
	  g.ares[nba].j = j;
	  g.ares[nba].poids = 1;
	  nba++;
	}
    }
  }
  
  g.nba = nba;
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

int composantes_connexes(gramaco g){
  
  int i, x, *parent, racine_j, *composante, nb_composante, racine, memoire;
  struct mallinfo2 m_debut, m_fin;
  m_debut = mallinfo2();
  
  parent = (int *)malloc(g.nbs * sizeof(int));
  composante = (int *)malloc(g.nbs * sizeof(int));
  memset(composante, 0, g.nbs * sizeof(int)); //initialiser toutes les cases du tableau composante à zéro
  nb_composante = 0;
  
  for (i = 0; i < g.nbs; i++){
    parent[i] = i;
  }  
  for (x = 0; x < g.nba; x++){
    racine = find(parent, g.ares[x].i);
    racine_j = find(parent, g.ares[x].j);
    if (racine !=racine_j){
      parent[racine_j] = racine;
    }
  }  
  for (i = 0; i < g.nbs; i++){
    racine = find(parent, i);
    if (composante[racine] == 0){
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

int dijkstra(gramaco g, int source){
  
  int u, x, i, step, v, nbs, *predecesseur, *vu, memoire;
  float *distance, poid;
  struct mallinfo2 m_debut, m_fin;
  m_debut = mallinfo2();
  
  nbs = g.nbs;
  distance = (float *)malloc(nbs * sizeof(float));
  predecesseur = (int *)malloc(nbs * sizeof(int));
  vu = (int *)calloc(nbs, sizeof(int));
  
  for (i = 0; i < nbs; i++){
    distance[i] = FLT_MAX;
    predecesseur[i] = -1;
  }
  distance[source] = 0.0;
  for (step = 0; step < nbs; step++){
    u = -1;
    for (i = 0; i < nbs; i++){
      if (!vu[i] && (u == -1 || distance[i] < distance[u])){
	u = i;
      }
    }
    if (distance[u] == FLT_MAX){
      break;
    }
    vu[u] = 1;
    for (x = 0; x < g.nba; x++){
      if (g.ares[x].i == u){
	v = g.ares[x].j;
	poid = g.ares[x].poids;
	if (!vu[v] && distance[u] + poid < distance[v]){
	  distance[v] = distance[u] + poid;
	  predecesseur[v] = u;
	}
      }
    }
  }
  m_fin = mallinfo2();
  
  /*for (i = 0; i < nbs; i++){
    if (distance[i] == FLT_MAX){
      printf("Sommet %d : Aucun chemin trouvé\n", i);
    }
    else{
      printf("Sommet %d : %.2f (précédent : %d)\n", i, distance[i], predecesseur[i]);
    }
  }*/
  memoire = (m_fin.uordblks - m_debut.uordblks);
  memoire = 0;
  
  free(distance);
  free(predecesseur);
  free(vu);

  return memoire;
}

void liberer_graphe(gramaco *g){
  
  if (g->ares != NULL){
    free(g->ares);
    g->ares = NULL;
  }
  g->nbs = g->nba = 0;
}

void affiche(gramaco g) {
    printf("Graphe : %d sommets, %d arêtes\n", g.nbs, g.nba);
}

void mesure(int nbs, int depart, int tour){
  
  int i;
  clock_t start, end;
  double time, total_time;
  int memoire_cree, memoire_dijk, memoire_compo, total_memoire, total_memoire_cree;
  gramaco g;
  total_time = 0.0;
  total_memoire = 0;
  
  for (i = 0; i < tour; i++){
    
    start = clock();
    g = creegraphe_mat3(nbs, &memoire_cree);
    
    memoire_dijk = dijkstra(g, depart);
    memoire_compo = composantes_connexes(g);
    end = clock();    
    time = (double)(end - start)/ CLOCKS_PER_SEC;

    total_memoire_cree += memoire_cree;
    printf("Mémoire grap en calcule : %d octet\n", total_memoire_cree);
    total_time += time;
    total_memoire += (memoire_cree + memoire_dijk + memoire_compo);  
    liberer_graphe(&g);
    }

  total_memoire_cree= total_memoire_cree/tour;
  total_time = total_time/tour;
  total_memoire = total_memoire/tour;

  printf("Mémoire grap mat3 : %d octet\n", total_memoire_cree);
  printf("Temps moyen : %.3f s\n", total_time);
  printf("Mémoire moyenne : %d octet\n", total_memoire);
}


int main(){
  
  int nbs, depart, tour;
  
  nbs = 1000;
  depart = 0;
  tour = 1;

  mesure(nbs, depart, tour);

  return 0;
}
