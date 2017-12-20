#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <time.h>

#define NS_PER_US 1000
typedef struct
{
   int box_id;
   int upper_left_y;
   int upper_left_x;
   int height;
   int width;
   int num_top_neighbors;
   int * top_neighbor_ids;
   int num_bottom_neighbors;
   int * bottom_neighbor_ids;
   int num_left_neighbors;
   int * left_neighbor_ids;
   int num_right_neighbors;
   int * right_neighbor_ids;
   double temperature; 
}Box;

bool isConverging();
void scanInput();
double getWeightedAverageTemperatures(int i);
int contactDistance(int s1, int s2, int s3, int s4);
void commitDSVs();
void computeDSV();


int numberOfGridBoxes=0;
Box * boxes;
double* weighted_avg_adjacent_temp;
double affectRate=0;
double epsilon=0;
double maxDSV;
double minDSV;
		
int main(int argc, char *argv[] ) {
	int iter = 0;
	if(argc > 3)
		printf("Too many arguments supplied.\n");
	affectRate = atof(argv[1]);
	epsilon = atof(argv[2]);
	scanInput();
	
	clock_t clock_time = 0;
	time_t time_time=0;
	clock_t begin, end;
    time_t start, finish;
    struct timespec start_c, end_c; 
	double diff; 
	
	clock_gettime(CLOCK_REALTIME,& start_c);
    time(&start);
    begin = clock();
	
	while(!isConverging()){
		computeDSV();		
		iter++;
	}
	
	end = clock();
    time(&finish);
    clock_gettime(CLOCK_REALTIME,& end_c); 
	
	diff = (double)( ((end_c.tv_sec - start_c.tv_sec)*CLOCKS_PER_SEC) + ((end_c.tv_nsec -start_c.tv_nsec)/NS_PER_US) );
    clock_time = end - begin;
    time_time = finish - start;
	
	printf("%s","***********************************************************************\n");
	printf("Dissipation converged in %d iterations,\n", iter);
	printf("with max DSV = %lf and min DSV = %lf\n", maxDSV, minDSV);
	printf("Affect rate = %lf , epsilon = %lf\n", affectRate, epsilon);
	printf("elapsed convergence loop time(clock) : %ld \n",clock_time);
	printf("elapsed convergence loop time(time) : %ld \n",time_time);
	printf("elapsed convergence loop time (chrono): %lf \n",diff);
	printf("***********************************************************************\n");
	
   return 0;
}

void scanInput(){
	int num_grid_rows=0;
	int num_grid_cols=0;
	int i;
	scanf("%d %d %d", & numberOfGridBoxes, & num_grid_rows, & num_grid_cols);
	
	boxes = (Box*)calloc(numberOfGridBoxes,sizeof(Box));
	weighted_avg_adjacent_temp = (double*)calloc(numberOfGridBoxes, sizeof(double));
	
	for(i=0;i<numberOfGridBoxes;i++)
	{
		scanf("%d", & (boxes[i].box_id));
		
		scanf("%d %d %d %d", & boxes[i].upper_left_y, & boxes[i].upper_left_x, & boxes[i].height, & boxes[i].width);
		int j;
		scanf("%d", & boxes[i].num_top_neighbors);
		boxes[i].top_neighbor_ids = (int *)calloc(boxes[i].num_top_neighbors, sizeof(int));	
		for(j=0;j< boxes[i].num_top_neighbors;j++)
		{
				scanf("%d", & boxes[i].top_neighbor_ids[j]);
		}
		
		scanf("%d", & boxes[i].num_bottom_neighbors);
		boxes[i].bottom_neighbor_ids = (int *)calloc(boxes[i].num_bottom_neighbors, sizeof(int));
		for(j=0;j<boxes[i].num_bottom_neighbors;j++)
		{
				scanf("%d", & boxes[i].bottom_neighbor_ids[j]);
		}
		
		scanf("%d", & boxes[i].num_left_neighbors);
		boxes[i].left_neighbor_ids = (int *)calloc(boxes[i].num_left_neighbors, sizeof(int));
		for(j=0;j<boxes[i].num_left_neighbors;j++)
		{
				scanf("%d", & boxes[i].left_neighbor_ids[j]);
		}
		
		scanf("%d", & boxes[i].num_right_neighbors);
		boxes[i].right_neighbor_ids = (int *)calloc(boxes[i].num_right_neighbors, sizeof(int));
		for(j=0;j<boxes[i].num_right_neighbors;j++)
		{
				scanf("%d", & boxes[i].right_neighbor_ids[j]);
				
		}		
		scanf("%lf", & boxes[i].temperature);
	}
}

double getWeightedAverageTemperatures(int i)
{	
		double weightedSumTemperature=0;
		int j;
		if(boxes[i].num_top_neighbors > 0)
			{
				for(j=0;j< boxes[i].num_top_neighbors;j++)
					{ 		
						Box topNeighbor = (Box)boxes[boxes[i].top_neighbor_ids[j]];						
						int distance = contactDistance(topNeighbor.upper_left_x, topNeighbor.upper_left_x + topNeighbor.width, boxes[i].upper_left_x, boxes[i].upper_left_x + boxes[i].width);			
						weightedSumTemperature=weightedSumTemperature + topNeighbor.temperature*distance;		
					}				
			}
		else
			weightedSumTemperature=weightedSumTemperature+ boxes[i].temperature*boxes[i].width;
			
		if(boxes[i].num_bottom_neighbors > 0)
			{
				for(j=0;j< boxes[i].num_bottom_neighbors;j++)
					{
						Box bottomNeighbor = (Box)boxes[boxes[i].bottom_neighbor_ids[j]];					
						int distance = contactDistance(bottomNeighbor.upper_left_x, bottomNeighbor.upper_left_x + bottomNeighbor.width, boxes[i].upper_left_x, boxes[i].upper_left_x + boxes[i].width);					
						weightedSumTemperature=weightedSumTemperature+ bottomNeighbor.temperature*distance;
												
					}
			}
		else
			weightedSumTemperature=weightedSumTemperature+ boxes[i].temperature*boxes[i].width;
		
		if(boxes[i].num_left_neighbors > 0)
			{
				for(j=0;j< boxes[i].num_left_neighbors;j++)
					{						
						Box leftNeighbor = (Box)boxes[boxes[i].left_neighbor_ids[j]];						
						int distance = contactDistance(leftNeighbor.upper_left_y, leftNeighbor.upper_left_y + leftNeighbor.height, boxes[i].upper_left_y, boxes[i].upper_left_y + boxes[i].height);						
						weightedSumTemperature=weightedSumTemperature+ leftNeighbor.temperature*distance;	
						
					}
			}
		else
			weightedSumTemperature=weightedSumTemperature+ boxes[i].temperature*boxes[i].height;
			
		if(boxes[i].num_right_neighbors > 0)
			{
				for(j=0;j< boxes[i].num_right_neighbors;j++)
					{
						Box rightNeighbor = (Box)boxes[boxes[i].right_neighbor_ids[j]];
						int distance = contactDistance(rightNeighbor.upper_left_y, rightNeighbor.upper_left_y + rightNeighbor.height, boxes[i].upper_left_y, boxes[i].upper_left_y + boxes[i].height);
						weightedSumTemperature=weightedSumTemperature+ rightNeighbor.temperature*distance;				
					}
			}
		else
			weightedSumTemperature=weightedSumTemperature+ boxes[i].temperature*boxes[i].height;
		
		double weightedAverageTemperature = weightedSumTemperature/(2 * (boxes[i].width + boxes[i].height));
		return weightedAverageTemperature;	
	}


void computeDSV()
{
	int i;
	for (i = 0; i < numberOfGridBoxes; i++)
	{
		weighted_avg_adjacent_temp[i] = getWeightedAverageTemperatures(i);
	}
	commitDSVs();
}

void commitDSVs()
{
	int i;
	for (i = 0; i < numberOfGridBoxes; i++)
	{
		if(boxes[i].temperature > weighted_avg_adjacent_temp[i])
			boxes[i].temperature= boxes[i].temperature - (boxes[i].temperature - weighted_avg_adjacent_temp[i]) * affectRate;
		else
			boxes[i].temperature= boxes[i].temperature + (weighted_avg_adjacent_temp[i] - boxes[i].temperature) * affectRate;
	}
}

int contactDistance(int s1, int s2, int s3, int s4)
	{
		int max = s3, min =s2;		
		if(s1>s3)
			max = s1;
		if(s2>s4)
			min = s4;
		
		if(max-min >0)
			return max-min;
		else
			return min-max;
	}

bool isConverging()
	{
		int i;
		 maxDSV = 0.0;
		 minDSV = boxes[0].temperature;
		for(i = 0; i < numberOfGridBoxes; i++)
		{
			if(maxDSV < boxes[i].temperature)
				maxDSV = boxes[i].temperature;
			if(minDSV > boxes[i].temperature)
				minDSV = boxes[i].temperature;
		}
		if(maxDSV-minDSV <= (epsilon * maxDSV))
			return true;
		return false;	
	}


