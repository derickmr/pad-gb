/*
  This program is an adaptation of the Mandelbrot program
  from the Programming Rosetta Stone, see
  http://rosettacode.org/wiki/Mandelbrot_set
  Compile the program with:
  gcc -o mandelbrot -O4 mandelbrot.c
  Usage:
 
  ./mandelbrot <xmin> <xmax> <ymin> <ymax> <maxiter> <xres> <out.ppm>
  Example:
  ./mandelbrot 0.27085 0.27100 0.004640 0.004810 1000 1024 pic.ppm
  The interior of Mandelbrot set is black, the levels are gray.
  If you have very many levels, the picture is likely going to be quite
  dark. You can postprocess it to fix the palette. For instance,
  with ImageMagick you can do (assuming the picture was saved to pic.ppm):
  convert -normalize pic.ppm pic.png
  The resulting pic.png is still gray, but the levels will be nicer. You
  can also add colors, for instance:
  convert -negate -normalize -fill blue -tint 100 pic.ppm pic.png
  See http://www.imagemagick.org/Usage/color_mods/ for what ImageMagick
  can do. It can do a lot.
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdint.h>
#include <pthread.h>
#define MAXROWS 10000000

//TODO allocate it dynamically
unsigned char colorsToBeWrittenOnFile[MAXROWS][6];

//TODO create struct and pass as parameter to threads
int arraySize;
int yres;
int xres;
int numThreads;
double xmin;
double xmax;
double ymin;
double ymax;
uint16_t maxiter;

void *calculate_mandelbrot(void *arg){
    
    int threadIndex = (intptr_t) arg;
    int counter = yres/numThreads * threadIndex;
    int counterEnd = yres/numThreads * (threadIndex + 1);

    int yStart = yres * threadIndex;
    int yEnd = yres * (threadIndex + 1);
        
    /* Precompute pixel width and height. */
     double dx=(xmax-xmin)/xres;
     double dy=(ymax-ymin)/yres;
    
    double x, y; /* Coordinates of the current point in the complex plane. */
    double u, v; /* Coordinates of the iterated point. */
    int i,j; /* Pixel counters */
    int k; /* Iteration counter */
    
    printf("yStart: %d, yEnd: %d, counter: %d, counterEnd: %d \n", yStart, yEnd, counter, counterEnd);
    
    for (j = yStart; j < yEnd && counter < counterEnd; j++) {
      y = ymax - j * dy;
      for(i = 0; i < xres; i++) {
      
        double u = 0.0;
        double v= 0.0;
        double u2 = u * u;
        double v2 = v*v;
        x = xmin + i * dx;
        /* iterate the point */
        for (k = 1; k < maxiter && (u2 + v2 < 4.0); k++) {
              v = 2 * u * v + y;
              u = u2 - v2 + x;
              u2 = u * u;
              v2 = v * v;
        };
        /* compute  pixel color and write it to file */
        if (k >= maxiter) {
          /* interior */
            int colorCounter;
            for (colorCounter = 0; colorCounter < 6; colorCounter++){
                colorsToBeWrittenOnFile[counter][colorCounter] = 0;
            }
          
        }
        else {
          /* exterior */
            colorsToBeWrittenOnFile[counter][0] = k >> 8;
            colorsToBeWrittenOnFile[counter][1] = k & 255;
            colorsToBeWrittenOnFile[counter][2] = k >> 8;
            colorsToBeWrittenOnFile[counter][3] = k & 255;
            colorsToBeWrittenOnFile[counter][4] = k >> 8;
            colorsToBeWrittenOnFile[counter][5] = k & 255;

        }
          counter++;
      }
    }
    return NULL;
  }
    

int main(int argc, char* argv[])
{
  /* Parse the command line arguments. */
  if (argc != 8) {
    printf("Usage:   %s <xmin> <xmax> <ymin> <ymax> <maxiter> <xres> <out.ppm>\n", argv[0]);
    printf("Example: %s 0.27085 0.27100 0.004640 0.004810 1000 1024 pic.ppm\n", argv[0]);
    exit(EXIT_FAILURE);
  }
    
    pthread_t threads[4];

  /* The window in the plane. */
    xmin = atof(argv[1]);
    xmax = atof(argv[2]);
    ymin = atof(argv[3]);
    ymax = atof(argv[4]);

  /* Maximum number of iterations, at most 65535. */
    maxiter = (unsigned short)atoi(argv[5]);

  /* Image size, width is given, height is computed. */
  xres = atoi(argv[6]);
  yres = (xres*(ymax-ymin))/(xmax-xmin);
    
    if (xres * yres > MAXROWS){
        printf("Image size not supported. Please lower it.");
        exit(EXIT_FAILURE);
    }
    
    arraySize = yres * xres;

  /* The output file name */
  const char* filename = argv[7];

  /* Open the file and write the header. */
  FILE * fp = fopen(filename,"wb");
  char *comment="# Mandelbrot set";/* comment should start with # */

  /*write ASCII header to the file*/
  fprintf(fp,
          "P6\n# Mandelbrot, xmin=%lf, xmax=%lf, ymin=%lf, ymax=%lf, maxiter=%d\n%d\n%d\n%d\n",
          xmin, xmax, ymin, ymax, maxiter, xres, yres, (maxiter < 256 ? 256 : maxiter));

    
    int i;
    
    numThreads = 4;

    for (i = 0; i < numThreads; i++){
        pthread_create(&(threads[i]), NULL, calculate_mandelbrot, (void *) (intptr_t) i);
    }
    
    for (i = 0; i < numThreads; i++){
        pthread_join(threads[i], NULL);
    }
        
    for (i = 0; i < arraySize; i++){
        fwrite(colorsToBeWrittenOnFile[i], 6, 1, fp);
    }
        
  fclose(fp);
  return 0;
}
