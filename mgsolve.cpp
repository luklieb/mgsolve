#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <errno.h>
#include <sys/time.h>
#include <math.h>
#include <vector>

#define GHOST false

// converts a timeval struct to seconds
double timevalToDouble( struct timeval *t ){
    return (double)(t->tv_sec) + (((double)(t->tv_usec))/1000000.0);
}

double residuum(int nx, int ny, std::vector<double> &grid, std::vector<double> &f_x_y, double h){
    double residuum = 0;
    double sum = 0;
    for(int i=1; i<ny-1; i++){
        for(int k=1; k<nx-1; k++){
                double temp = ((-1.0f/h)*(grid[i*nx+(k-1)]+grid[i*nx+(k+1)]) - (1.0f/h)*(grid[(i-1)*nx+k]+grid[(i+1)*nx+k]) + ((2.0f/h+2.0f/h)*grid[i*nx+k]) - f_x_y[i*nx+k]);
                sum += temp*temp;
        }
    }
    residuum = sqrt((1.0f/(nx*ny))*sum);
    return residuum;
}

void Red_Black_Gauss(int nx, int ny, std::vector<double> &grid, std::vector<double> &f_x_y, double h){
    
     double const faktor = 1.0f/(2.0f/h + 2.0f/h);
     double const beta = 1.0f/h;
     double const gamma = 1.0f/h;

     for(int iterations=0; iterations<nx*ny; iterations++){
         for(int m=1; m<ny-1; m++){
             for(int q=1; q<nx-1; q++){
                grid[m*nx+q] = (faktor) * (f_x_y[m*nx+q] + (beta)*(grid[m*nx+(q-1)]+grid[m*nx+(q+1)]) + (gamma)*(grid[(m-1)*nx+q]+grid[(m+1)*nx+q]));
             }
         }
     }
}

// test liefert das gleiche wie andere version
void Red_Black_Gauss_v(int nx, int ny, std::vector<double> &grid, std::vector<double> &f_x_y, double h){
    
     for(int iterations=0; iterations<nx*ny; iterations++){
         for(int m=1; m<ny-1; m++){
             for(int q=1; q<nx-1; q++){
                grid[m*nx+q] = (1.0/4.0) * (h*h*f_x_y[m*nx+q] + (grid[m*nx+(q-1)]+grid[m*nx+(q+1)]+grid[(m-1)*nx+q]+grid[(m+1)*nx+q]));
             }
         }
     }
}

int main(int argc, char **argv){

	// Ueberpruefung, ob Eingabeparamter passen
    if(argc != 3){
		fprintf(stderr, "Usage: ./mgsolve l n\n");
		exit(EXIT_SUCCESS);
	}
	
    // Einlesen der Eingabeparameter l, n
    int l, n;
    l = atoi(argv[1]);
    n = atoi(argv[2]);
    fprintf(stderr, "%d/n", n);
	
	// nx and ny are the total points in x and y direction
    int nx = (int)(pow(2,l)+1);
    int ny = (int)(pow(2,l)+1);
    double h = 1.0/(nx-1);

	// Bei Neumann Randwerte werden in x-Richtung ghost-Layers benoetigt
    if(GHOST){
        nx=nx+2;
    }
	
	// Anlegen von f_x_y:
    std::vector<double> f_x_y;
	
    // Speicher fuer das Gitter allokieren:
    std::vector<double> grid(nx*ny, 0.0);

    // Initialisierung des Gitters
    if(GHOST){
        f_x_y = std::vector<double>(nx*ny, 2.0);
		// Dirichlet RDB setzten aber nicht auf den Ghost-Layers 
        for(int i=0; i<nx-2; i++){
            grid[i+1] = i*h*(1-i*h);
            grid[nx*(ny-1)+i+1] = i*h*(1-i*h);
        }
        for(int i=0; i<ny; i++){
            grid[i*nx] = -h;
            grid[i*nx+(nx-1)] = -h;
        }
    }
    else{
        f_x_y = std::vector<double>(nx*ny, 0.0);
        for(int i=0; i<nx; i++){
            grid[nx*(ny-1)+i] = sin(M_PI*i*h)*sinh(M_PI);
        }
    }

	// TIMER RED_BLACK_GAUSS
    struct timeval t1;
    gettimeofday(&t1, NULL);

    Red_Black_Gauss(nx, ny, grid, f_x_y, h);

    struct timeval t2;
    gettimeofday(&t2, NULL);
    fprintf(stdout, "Timer Red_Black_Gauss: %lf\n", timevalToDouble(&t2)-timevalToDouble(&t1));

    double res = residuum(nx, ny, grid, f_x_y, h);

    // Oeffnen der Ausgabedatei
	FILE *out = fopen("./solution.txt", "w");
	if(out == NULL){
		perror("fopen");
		exit(EXIT_FAILURE);
    }

    // Ausgabe fuer solution.txt
    fprintf(out, "# x y u(x,y)\n");
    for(int j=0; j<=ny; j++){
        for(int t=0; t<=nx; t++){
            double f = (double)t/(double)nx;
            double q = (double)j/(double)ny;
// warum 2x f ??
            fprintf(out, "%.5lf %.5lf %.8lf\n", f*2, q, grid[j*ny+t]);
        }
		fprintf(out, "\n");
	}
	fclose(out);
    fprintf(stdout, "Residuum: %lf\n", res);
}
