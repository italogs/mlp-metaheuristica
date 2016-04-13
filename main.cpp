#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <limits.h>
#include <float.h>
#include <inttypes.h>
#include <unistd.h>
#include <algorithm>
#include <vector>
#include <time.h>
#include "library/mtwist-1.5/mtwist.h"

using namespace std;

double **avaliacMatrix = NULL;

typedef struct Vertex{
    int index;
    double xCoord,yCoord;
}Vertex;

typedef struct TList{
	int origin, destiny;
	double cost;
}TList;
TList **lista = NULL;

double calculate_EUC_2D(Vertex i, Vertex j) {
	return sqrt(pow( i.xCoord - j.xCoord   ,2) + pow( i.yCoord - j.yCoord   ,2) ) ;
}
double costFromIndex(int *solution,int index){
	double cost = 0.0;
	int i=0 , k = index-1;
	for (i = 0; i < index-1; ++i) {
		cost = cost + (k * avaliacMatrix[solution[i]][solution[i+1]]);
		k--;
	}
	return cost;
}
double calculateCostSolution(int *solution, int lenSol) {
	double cost = 0.0;
	int i=0 , k = lenSol-1;
	for (i = 0; i < lenSol-1; ++i) {
		cost = cost + (k * avaliacMatrix[solution[i]][solution[i+1]]);
		k--;
	}
	return cost;
}
bool sortDistances(TList i, TList j) { return (i.cost <= j.cost);}
bool sortSavings(TList i, TList j) { return (i.cost >= j.cost);}

double Constructive_method2(int maxIterations,int *solution,int lenSol) {
	double mlp = DBL_MAX;
	unsigned int i,j,k;
	int m,n;
	vector<TList> savings;

	TList element;
	for (m = 1; m < lenSol-1; m++) {
		for (n = m+1; n < lenSol; n++) {
			element.origin = m;
			element.destiny = n;
			element.cost = avaliacMatrix[0][m] + avaliacMatrix[0][n] - avaliacMatrix[m][n];
			savings.push_back(element);
		}
	}
	sort(savings.begin(),savings.end(),sortSavings);
	printf("%d \n", (int) savings.size() );
	i=0;
	vector<TList>::iterator it,it2;
	vector<int> subroute;
	vector< vector < int > > subroutes;
	for (m = 1; m < lenSol; ++m){
		subroute.push_back(0);
		subroute.push_back(m);
		subroute.push_back(0);
		subroutes.push_back(subroute);
		subroute.clear();
	}
	
	solution[0] = 0;


	int index_origin = -1, index_destiny = -1;

	// printf("size subroutes %d \n", subroutes.size() );
	while(i < savings.size()){
		index_origin = index_destiny = -1;
		for (j = 0; j < subroutes.size(); j++) {
			for (k = 1; k < subroutes[j].size()-1; k++) {
				if(subroutes[j][k] == savings[i].origin)
					index_origin = j;
				if(subroutes[j][k] == savings[i].destiny)
					index_destiny = j;
			}
		}
		if(index_origin != index_destiny){
				subroutes[index_origin].insert(subroutes[index_origin].begin()+(subroutes[index_origin].size()-1) ,
												subroutes[index_destiny].begin()+1,subroutes[index_destiny].end());
				subroutes[index_origin].erase(subroutes[index_origin].begin()+subroutes[index_origin].size()-1);
				subroutes.erase(subroutes.begin()+index_destiny);
		}

		if(subroutes.size() == 1) break;
		i++;
	}
	for (m = 0; m < lenSol; ++m){
		solution[m] = subroutes[0][m];
	}

	mlp = calculateCostSolution(solution,lenSol);
	

	return mlp;
}

double constructive_method1(int maxIterations, int *solution, int lenSol){
    int i, j, k, l,iter=0;
    double bestMLP = DBL_MAX , newMLP,accMLP = 0.0;

    TList element;
    vector <TList> distances;
    vector <TList> distancesLRC;
    int *bestSolution = (int *)malloc(sizeof(int)*lenSol);
    unsigned int lengthLRC, it;

    int origin, destiny;
    
    for ( iter = 0; iter < maxIterations; iter++) {
	  	mt_seed32new(iter+1);
	    int numPointsSol = 0;
	    // printf("inicio\n");
	    double alpha = (  mt_lrand() % 26) / 100.0;
	    /* Inicializa lista com todas as arestas e seus respectivos comprimentos (custos) */
	    for(i = 0; i < lenSol; i++){
	        for(j = (i + 1); j < lenSol; j++){
	            element.origin = i;
	            element.destiny = j;
	            element.cost = avaliacMatrix[i][j];
	            distances.push_back(element);
	        }
	    }
	    /* Ordena lista conforme os comprimentos das arestas */
	    // sort(distances.begin(), distances.end(), sortDistances);

	     sort(distances.begin(), distances.begin()+(distances.size()-1), sortDistances);

	    /* Inicializa a LRC com os pontos mais próximos da origem (depósito), conforme o comprimento das arestas
	        entres esses pontos (destino) e a origem, e seleciona um ponto na LRC randomicamente */
	    for(it = 0; it < distances.size(); it++){
	        if((distances[it].origin == 0) || (distances[it].destiny == 0)){
	            distancesLRC.push_back(distances[it]);
	        }
	    }
	    if(alpha == 0){
	        l = 0;
	    }else{
	        lengthLRC = ceil(alpha * distancesLRC.size());
	        l = mt_lrand() % lengthLRC;
	        
	    }
	    
	    
	    origin = distancesLRC[l].origin;
	    destiny = distancesLRC[l].destiny;
	    
	    if(destiny == 0){
	        destiny = origin;
	        origin = 0;
	    }
	    
	    /* Adiciona os dois primeiros pontos a solução */
	    solution[0] = origin;
	    solution[1] = destiny;
	    numPointsSol = 2;
	    k = numPointsSol;
	    // printf("ates while\n");
	    while(k < lenSol){
	        distancesLRC.clear();
	        for(it = 0; it < distances.size(); it++){
	            if((distances[it].origin == solution[k - 1]) || (distances[it].destiny == solution[k - 1])){
	                distancesLRC.push_back(distances[it]);
	            }
	        }

	        while(distancesLRC.size() > 0){
	            if(alpha == 0){
	                l = 0;
	            }else{
	                lengthLRC = ceil(alpha * distancesLRC.size());
	                // l = rand() % lengthLRC;
	                l = mt_lrand() % lengthLRC;
	            }
	            //printf("lengthLRC = %d\n", lengthLRC);
	            //printf("l = %d\n", l);
	            origin = distancesLRC[l].origin;
	            destiny = distancesLRC[l].destiny;

	            if(destiny == solution[k - 1]){
	                destiny = origin;
	            }

	            j = 0;
	            while((j < numPointsSol) && (destiny != solution[j])){
	                j++;
	            }
	            
	            if(j == numPointsSol){
	                break;
	            }else{
	                distancesLRC.erase(distancesLRC.begin() + l);
	            }
	        }
	        
	        solution[k] = destiny;
	        k++;
	        numPointsSol++;    
	    }

	    /* Calcula o custo total da solução construida */
	    newMLP = calculateCostSolution(solution, numPointsSol);
	    distances.clear();
	    distancesLRC.clear();
	    // printf("newMLP %lf\n",newMLP );
	    if( newMLP < bestMLP) {
	    	bestMLP = newMLP;
	    	for(i = 0; i < lenSol; i++)
	    		bestSolution[i] = solution[i];
	    }
	    accMLP +=newMLP;
    }
	for(i = 0; i < lenSol; i++) {
		solution[i] = bestSolution[i];
	}
	printf("avg: %lf \n", accMLP/(1.0*maxIterations));
    free(bestSolution);
    return bestMLP;
}


int main(int argc, char *argv[]){

	FILE *inputFile = fopen("instancias/att48.tsp","r");
	if(inputFile==NULL) {
		printf("erro ao ler o arquivo\n");
		exit(0);
	}
	int lenSol = 0,i,j, *solution = NULL, maxIterations = 10; 
	char str1[256];
	Vertex *v_vector = NULL;
	clock_t cInit,cEnd;
	int n_read;
	while(!feof(inputFile)) {
		n_read = fscanf(inputFile, "%s",str1);

		printf("%s\n", str1 );
		if(n_read != 1) {printf("erro fscanf(inputFile, \"%%s\",str1)");exit(0);}
		if(strcmp(str1,"DIMENSION:") == 0){
			n_read = fscanf(inputFile, "%d",&lenSol);
			if(n_read != 1) {printf("erro fscanf(inputFile, \"%%d\",&lenSol);");exit(0);}
			lenSol = lenSol;
			solution = (int *)malloc(sizeof(int)*lenSol);
			v_vector = (Vertex *)malloc(sizeof(Vertex)*lenSol);
		} else if(strcmp(str1,"NODE_COORD_SECTION") == 0) {
			i = 0;
			while(i < lenSol){
				n_read = fscanf(inputFile, "%d %lf %lf", &v_vector[i].index, &v_vector[i].xCoord, &v_vector[i].yCoord);
				if(n_read != 3) {printf("erro fscanf(inputFile, \"%%d %%lf %%lf\", &v_vector[i].index, &v_vector[i].xCoord, &v_vector[i].yCoord);");exit(0);}
				i++;
			}
		}
	}
	avaliacMatrix = (double **)malloc(sizeof(double *) * lenSol);
	for (i = 0; i < lenSol; ++i)
		avaliacMatrix[i] = (double *)malloc(sizeof(double) * lenSol);

	for (i = 0; i < lenSol; ++i)
		for (j = 0; j < lenSol; ++j)
			avaliacMatrix[i][j] = calculate_EUC_2D(v_vector[i],v_vector[j]);

	double optLat = 603910.000;
	printf("Opt Lat %lf\n",optLat);


	printf("**** Constructive_method1 - START ****\n");

	cInit = clock();
	double lat = constructive_method1(maxIterations,solution,lenSol);
	printf("Latency: %lf gap: %.2lf%% = ((optLat - lat) / lat) \n", lat, ( optLat - lat )/ lat * 100.0 );
	cEnd = clock();
	printf("Tempo: %lf\n", (double)difftime(cEnd, cInit) / (CLOCKS_PER_SEC));

	printf("**** Constructive_method1 - END ****\n");




	printf("**** Constructive_method2 - START ****\n");

	cInit = clock();
	lat = Constructive_method2(maxIterations,solution,lenSol);
	printf("Constructive_method2: %lf gap: %.2lf%% = ((optLat - lat) / lat) \n", lat, ( optLat - lat )/ lat * 100.0 );
	cEnd = clock();
	printf("Tempo: %lf\n", (double)difftime(cEnd, cInit) / (CLOCKS_PER_SEC));

	printf("**** Constructive_method2 - END ****\n");


	fclose(inputFile);

	for (i = 0; i < lenSol; ++i)
		free(avaliacMatrix[i]);

	free(avaliacMatrix);
	free(solution);
	free(v_vector);

	return 0;


}