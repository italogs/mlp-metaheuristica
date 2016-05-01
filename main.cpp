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

#define TRUE 1
#define FALSE 0

using namespace std;

FILE *outputFile = NULL;
double **avaliacMatrix = NULL;
double constructive_avg = 0.0;
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

double constructive_method2(int *solution,int lenSol) {
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

double constructive_method1(int *solution, int lenSol){
    int i, j, k, l;
    double newMLP;

    TList element;
    vector <TList> distances;
    vector <TList> distancesLRC;

    unsigned int lengthLRC, it;

    int origin, destiny;
    
    // for ( iter = 0; iter < maxIterations; iter++) {
	  	
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

    return newMLP;
}

double swap(int *originalSolution,int *currentSolution,int lenSol,int pivot,int pivot2) {
	int k,aux;
    for(k = 0; k < lenSol; k++){
        currentSolution[k] = originalSolution[k];
    }
    
    aux = currentSolution[pivot];
    currentSolution[pivot] = currentSolution[pivot2];
    currentSolution[pivot2] = aux;

    return calculateCostSolution(currentSolution,lenSol);
}


double two_opt(int *originalSolution,int *currentSolution,int lenSol,int pivot,int pivot2) {
	int k;
	for(k = 0; k < pivot; k++)
        currentSolution[k] = originalSolution[k];

    for(k = pivot; k <= pivot2; k++)
        currentSolution[k] = originalSolution[(pivot2 - k) + pivot];

    for(k = (pivot2 + 1); k < lenSol; k++)
        currentSolution[k] = originalSolution[k];
	return calculateCostSolution(currentSolution,lenSol);
}

double reinsertion(int *originalSolution,int *currentSolution,int lenSol,int pivot,int pivot2) {
	int k,aux;
	if(pivot < pivot2){
        for(k = 0; k < pivot; k++){
            currentSolution[k] = originalSolution[k];
        }

        aux = originalSolution[pivot];
        k = pivot;
        while(k < pivot2){
            currentSolution[k] = originalSolution[k + 1];
            k++;
        }

        currentSolution[k] = aux;
        for(k = (pivot2 + 1); k < lenSol; k++){
            currentSolution[k] = originalSolution[k];
        }
    }
    else if(pivot > pivot2){
        currentSolution = new int[lenSol];

        for(k = 0; k < pivot2; k++){
            currentSolution[k] = originalSolution[k];
        }

        aux = originalSolution[pivot];
        k = pivot;
        while(k > pivot2){
            currentSolution[k] = originalSolution[k - 1];
            k--;
        }

        currentSolution[k] = aux;
        for(k = (pivot + 1); k < lenSol; k++){
            currentSolution[k] = originalSolution[k];
        }
    } else {
    	return DBL_MAX;
    }
	return calculateCostSolution(currentSolution,lenSol);
}


double ms_local_search_best_improving(int maxIterations, int *bestGlobalSolution, int lenSol){
	int i,j,k;
	bool improvement_2opt, improvement_reinsertion, improvement_swap;
	double currentLat, neighborLat,bestLat = DBL_MAX,bestGlobalLat = DBL_MAX;
	int *bestSolution = new int[lenSol], *currentSolution = new int[lenSol], *neighborSolution = new int[lenSol];
	for (int iter = 0; iter < maxIterations; iter++) {
		currentLat = constructive_method1(currentSolution,lenSol);

		improvement_2opt = improvement_reinsertion = improvement_swap = TRUE;
		while (improvement_2opt || improvement_reinsertion || improvement_swap) {
			if(improvement_swap) {
				improvement_swap = FALSE;
				for(i = 1; i < (lenSol - 2); i++){
		        	for(j = (i + 1); j < (lenSol - 1); j++){
						neighborLat = swap(currentSolution,neighborSolution,lenSol,i,j);
						if(neighborLat < bestLat) {
							for (k = 0; k < lenSol; k++)
								bestSolution[k] = neighborSolution[k];
							bestLat = neighborLat;
							improvement_2opt = improvement_reinsertion = improvement_swap = TRUE;
						}
					}
				}
			} else if(improvement_reinsertion) {
				improvement_reinsertion = FALSE;
				for (i = 1; i < (lenSol - 1); i++) {
					for (j = 1; j < (lenSol - 1); j++) {
						neighborLat = reinsertion(currentSolution,neighborSolution,lenSol,i,j);
						if(neighborLat < bestLat) {
							for (k = 0; k < lenSol; k++)
								bestSolution[k] = neighborSolution[k];
							bestLat = neighborLat;
							improvement_2opt = improvement_reinsertion = improvement_swap = TRUE;
						}
					}
				}
			} else if(improvement_2opt) {
				improvement_2opt = FALSE;
				for (i = 1; i < (lenSol - 2); i++) {
					for (j = (i+1); j < (lenSol - 1); j++) {
						neighborLat = two_opt(currentSolution,neighborSolution,lenSol,i,j);
						if(neighborLat < bestLat) {
							for (k = 0; k < lenSol; k++)
								bestSolution[k] = neighborSolution[k];
							bestLat = neighborLat;
							improvement_2opt = improvement_reinsertion = improvement_swap = TRUE;
						}
					}
				}
			}
			if(bestLat < currentLat) {
				for (i = 0; i < lenSol; i++)
					currentSolution[i] = bestSolution[i];
				currentLat = bestLat;
			}
		}

		if(currentLat < bestGlobalLat) {
			for (i = 0; i < lenSol; i++)
				bestGlobalSolution[i] = currentSolution[i];
			bestGlobalLat = currentLat;
		}
	}
	delete[] neighborSolution;
	delete[] bestSolution;
	delete[] currentSolution;
	return bestGlobalLat;
}

double ms_local_search_first_improving(int maxIterations, int *bestGlobalSolution, int lenSol){

	int i,j,k;
	double currentLat, neighborLat,bestGlobalLat = DBL_MAX;
	int *currentSolution = new int[lenSol], *neighborSolution = new int[lenSol];
	bool improvement_2opt, improvement_reinsertion, improvement_swap;
	for (int iter = 0; iter < maxIterations; iter++) {
		currentLat = constructive_method1(currentSolution,lenSol);

		improvement_2opt = improvement_reinsertion = improvement_swap = TRUE;
		while (improvement_2opt || improvement_reinsertion || improvement_swap) {
			if(improvement_swap) {
				improvement_swap = FALSE;
					for(i = 1; i < (lenSol - 2); i++){
			        	for(j = (i + 1); j < (lenSol - 1); j++){
							neighborLat = swap(currentSolution,neighborSolution,lenSol,i,j);
							if(neighborLat < currentLat) {
								for (k = 0; k < lenSol; k++)
									currentSolution[k] = neighborSolution[k];
								currentLat = neighborLat;
								improvement_2opt = improvement_reinsertion = improvement_swap = TRUE;
							}
						}
					}
			} else if (improvement_reinsertion) {
				improvement_reinsertion = FALSE;
				for (i = 1; i < (lenSol - 1); i++) {
					for (j = 1; j < (lenSol - 1); j++) {
						neighborLat = reinsertion(currentSolution,neighborSolution,lenSol,i,j);
						if(neighborLat < currentLat) {
							for (k = 0; k < lenSol; k++)
								currentSolution[k] = neighborSolution[k];
							currentLat = neighborLat;
							improvement_2opt = improvement_reinsertion = improvement_swap = TRUE;
						}
					}
				}
			} else if(improvement_2opt) {
				improvement_2opt = FALSE;
				for (i = 1; i < (lenSol - 2); i++) {
					for (j = (i+1); j < (lenSol - 1); j++) {
						neighborLat = two_opt(currentSolution,neighborSolution,lenSol,i,j);
						if(neighborLat < currentLat) {
							for (k = 0; k < lenSol; k++)
								currentSolution[k] = neighborSolution[k];
							currentLat = neighborLat;
							improvement_2opt = improvement_reinsertion = improvement_swap = TRUE;
						}
					}
				}
			}
		}

		if(currentLat < bestGlobalLat) {
			for (i = 0; i < lenSol; i++)
				bestGlobalSolution[i] = currentSolution[i];
			bestGlobalLat = currentLat;
		}

	}
	delete[] neighborSolution;
	delete[] currentSolution;
	return bestGlobalLat;
}

int main(int argc, char *argv[]){
	
	int exec = 1;
	if(argc < 2 || argc > 3) {
		printf("input: ./main instance/example.tsp |<exec>\n");
		return 0;
	}
	if(argc == 3) {
		exec = atoi(argv[2]);
	}
	if(exec <= 0) {
		printf("exec must be greater than 0\n");
		return 0;
	}

	// mt_seed32new(exec);
	// long seed = time(NULL);
	long seed = exec;


	
	char filenameInput[128] = "";
	strcpy(filenameInput,argv[1]);
	printf("Seed %ld \n",seed );
	FILE *inputFile = fopen(filenameInput,"r");
	if(inputFile == NULL) {
		printf("erro ao ler o arquivo\n");
		exit(0);
	}
	int lenSol = 0,i,j, *solution = NULL, maxIterations = 20; 
	char str1[256];
	Vertex *v_vector = NULL;
	
	double optLat = DBL_MAX;
	int n_read;

	while(!feof(inputFile)) {
		n_read = fscanf(inputFile, "%s",str1);

		// printf("%s\n", str1 );
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
		} else if(strcmp(str1,"OPT:") == 0){
			n_read = fscanf(inputFile, "%lf",&optLat);
			if(n_read != 1) {printf("erro fscanf(inputFile, \"%%lf\",&optLat);");exit(0);}
		}
	}
	avaliacMatrix = (double **)malloc(sizeof(double *) * lenSol);
	for (i = 0; i < lenSol; ++i)
		avaliacMatrix[i] = (double *)malloc(sizeof(double) * lenSol);

	for (i = 0; i < lenSol; ++i)
		for (j = 0; j < lenSol; ++j)
			avaliacMatrix[i][j] = calculate_EUC_2D(v_vector[i],v_vector[j]);

	// double optLat = 603910.000;
	// printf("Opt Lat %lf\n",optLat);




	char filenameOutput[128] = "outputs/";

	strcat(filenameOutput,basename(filenameInput));
	strcat(filenameOutput,".csv");
	printf("filenameOutput %s\n", filenameOutput);

	FILE *all_outputsFile = fopen("outputs/all_outputs.csv","a");
	if(exec == 1) {
		outputFile = fopen(filenameOutput,"w");
	} else {
		outputFile = fopen(filenameOutput,"a");
	}
	if(outputFile == NULL) {
		printf("outputFile is not able to write\n");
		return 0;
	}
	if(all_outputsFile == NULL) {
		printf("all_outputs is not able to write\n");
		return 0;
	}
	if(exec == 1) {
		fprintf(outputFile, "Best Improving,,First Improving\n");
		fprintf(outputFile, "Best,Time,Best,Time\n");
	}

	clock_t cInitBest_Improving,cEndBest_Improving;
	mt_seed32new(seed);
	printf("**** local_search_best_improving - START ****\n");
	cInitBest_Improving = clock();
	double best_improving_latency = ms_local_search_best_improving(maxIterations,solution,lenSol);
	printf("Latency: %lf \n",best_improving_latency );
	cEndBest_Improving = clock();
	double best_improving_time = (double)difftime(cEndBest_Improving, cInitBest_Improving) / (CLOCKS_PER_SEC);
	printf("Tempo: %lf\n",best_improving_time );
	for (i = 0; i < lenSol; ++i) {
		printf("%d ", solution[i] );
	}
	printf("\n");
	printf("**** local_search_best_improving - END ****\n");


	mt_seed32new(seed);
	clock_t cInitFirst_Improving,cEndFirst_Improving;
	printf("**** local_search_first_improving - START ****\n");
	cInitFirst_Improving = clock();
	double first_improving_latency = ms_local_search_first_improving(maxIterations,solution,lenSol);
	printf("Latency: %lf \n",first_improving_latency );
	cEndFirst_Improving = clock();
	double first_improving_time = (double)difftime(cEndFirst_Improving, cInitFirst_Improving) / (CLOCKS_PER_SEC);
	printf("Tempo: %lf\n",first_improving_time );
	for (i = 0; i < lenSol; ++i) {
		printf("%d ", solution[i] );
	}
	printf("\n");
	printf("**** local_search_first_improving - END ****\n");

	fprintf(outputFile, "%lf,%lf,%lf,%lf\n",
			best_improving_latency,best_improving_time,
			first_improving_latency,first_improving_time);

	fprintf(all_outputsFile, "%s,%lf,%lf,%lf,%lf\n",
			basename(filenameInput),
			best_improving_latency,best_improving_time,
			first_improving_latency,first_improving_time);



	printf("FIM\n");
	

	fclose(inputFile);
	fclose(outputFile);
	fclose(all_outputsFile);
	for (i = 0; i < lenSol; ++i)
		free(avaliacMatrix[i]);

	free(avaliacMatrix);
	free(solution);
	free(v_vector);

	return 0;


}