// MIT License
// 
// Copyright (c) 2023 Trevor Bakker 
// 
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:
// 
// The above copyright notice and this permission notice shall be included in all
// copies or substantial portions of the Software.
// 
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
// SOFTWARE.

#include <time.h>
#include <sys/time.h>
#include <ctype.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>
#include <getopt.h>
#include <stdint.h>
#include "utility.h"
#include "star.h"
#include "float.h"
#include <pthread.h>

#define NUM_STARS 30000 
#define MAX_LINE 1024
#define DELIMITER " \t\n"

struct Star star_array[ NUM_STARS ];
uint8_t   (*distance_calculated)[NUM_STARS];

double  min  = FLT_MAX;
double  max  = FLT_MIN;
double  mean = 0;
uint64_t count = 0;

pthread_mutex_t mutex;

int NUM_THREADS = 100;


void showHelp()
{
  printf("Use: findAngular [options]\n");
  printf("Where options are:\n");
  printf("-t          Number of threads to use\n");
  printf("-h          Show this help\n");
}


void* determineAverageAngularDistance(void* arg)
{
    uint32_t i, j;

    int* current_thread = (int*)arg;
    printf("Current thread number is : %d\n", *current_thread + 1);

    int start = *current_thread * NUM_STARS / NUM_THREADS;
    int end = (*current_thread + 1) * NUM_STARS / NUM_THREADS;

    printf("Here we will check %d to %d\n", start, end);

    for (i = 0; i < NUM_STARS; i++)
    {
      for (j = start; j < end; j++)
      {        
        if( i!=j && distance_calculated[i][j] == 0 && distance_calculated[j][i] == 0)
        {
          distance_calculated[i][j] = 1;

          if (distance_calculated[j][i] == 1)
            continue;

          double distance = calculateAngularDistance( star_array[i].RightAscension, star_array[j].Declination,
                                                      star_array[j].RightAscension, star_array[i].Declination ) ;

          if (distance_calculated[j][i] == 1)
            continue;

          if (min > distance)
            min = distance;

          if (max < distance)
            max = distance;

          distance_calculated[j][i] = 1;

          pthread_mutex_lock(&mutex);
          count++;
          mean = mean + (distance - mean) / count;
          pthread_mutex_unlock(&mutex);
        }
      }
    }
    printf("Done checking %d to %d\n", start, end);
    pthread_exit(NULL);
}


int main( int argc, char * argv[] )
{
  FILE *fp;
  uint32_t star_count = 0;

  uint32_t n;

  distance_calculated = malloc(sizeof(uint8_t[NUM_STARS][NUM_STARS]));

  if (distance_calculated == NULL)
  {
    uint64_t num_stars = NUM_STARS;
    uint64_t size = num_stars * num_stars * sizeof(uint8_t);
    printf("Could not allocate %ld bytes\n", size);
    exit( EXIT_FAILURE );
  }

  // default every thing to 0 so we calculated the distance.
  memset(distance_calculated, 0, sizeof(distance_calculated[0][0]) * NUM_STARS * NUM_STARS);

  for (n = 1; n < argc; n++)          
  {
    if (strcmp(argv[n], "-h" ) == 0)
    {
      showHelp();
      exit(0);
    }
    if (strcmp(argv[n], "-t" ) == 0)
    {
      if (argv[n + 1] == NULL)
      {
        printf("Invalid use of '-t' command.\nUse: -t Number of threads to use\n");
        exit(0);
      }
      else
      {
        if (atoi(argv[n + 1]) == 0)
        {
          printf("'%s' is not a number.\nUse: -t Number of threads to use\n", argv[n + 1]);
          exit(0);
        }

        for (int i = 0; i < strlen(argv[n + 1]); i++)
        {
          if(!isdigit(argv[n + 1][i]))
          {
            printf("'%s' is an invalid amount of threads\n", argv[n + 1]);
            exit(0);
          }
        }

        NUM_THREADS = atoi(argv[n + 1]);
      }
    }
  }

  fp = fopen("data/tycho-trimmed.csv", "r");

  if (fp == NULL)
  {
    printf("ERROR: Unable to open the file data/tycho-trimmed.csv\n");
    exit(1);
  }

  char line[MAX_LINE];
  while (fgets(line, 1024, fp))
  {
    uint32_t column = 0;

    char* tok;
    for (tok = strtok(line, " ");
            tok && *tok;
            tok = strtok(NULL, " "))
    {
       switch(column)
       {
          case 0:
              star_array[star_count].ID = atoi(tok);
              break;
       
          case 1:
              star_array[star_count].RightAscension = atof(tok);
              break;
       
          case 2:
              star_array[star_count].Declination = atof(tok);
              break;

          default: 
             printf("ERROR: line %d had more than 3 columns\n", star_count);
             exit(1);
             break;
       }
       column++;
    }
    star_count++;
  }
  printf("%d records read\n", star_count);

  pthread_mutex_init(&mutex, NULL);

  pthread_t tid[NUM_THREADS];
  int thread_num[NUM_THREADS];

  struct timeval start, end;
  gettimeofday(&start, 0);

  for (int i = 0; i < NUM_THREADS; i++)
  {
    thread_num[i] = i;
    if(pthread_create(&tid[i], NULL, determineAverageAngularDistance, (void*)&thread_num[i]))
    {
      printf("ERROR: Thread creation failed\n");
			exit(EXIT_FAILURE);
    }
  }

  for (int i = 0; i < NUM_THREADS; i++)
  {
    if(pthread_join(tid[i], NULL))
    {
      printf("ERROR: Thread joining failed\n");
			exit(EXIT_FAILURE);
    }
  }

  gettimeofday(&end, 0);
  long seconds = end.tv_sec - start.tv_sec;
  long microseconds = end.tv_usec - start.tv_usec;
  double elapsed = seconds + microseconds*1e-6;

  pthread_mutex_destroy(&mutex);

  // Print the values of the mean, min, and max that were found using threads.
  printf("\nAverage distance found is %lf\n", mean );
  printf("Minimum distance found is %lf\n", min );
  printf("Maximum distance found is %lf\n", max );

  printf("\nTime taken to calculate average angular distance: %lf seconds\n", elapsed);

  return 0;
}

