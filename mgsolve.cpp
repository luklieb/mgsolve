#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <errno.h>
#include <math.h>
#include <vector>
#include <fstream>
#include <assert.h>
#include <sys/time.h>


#define GHOST true

typedef double type;
typedef std::vector<type> typeVec;
typedef std::vector<int> intVec;

// TEST _____________________________________________________________________________________________________________


template<typename T> class grid{

private:
	size_t lengthInX;
	size_t lengthInY;
	std::vector<T> data;

public:
	grid(){lengthInX = 0;lengthInY = 0;}
	grid(size_t xDim, size_t yDim){ data = std::vector<T>(xDim * yDim, (T)(0.0)); lengthInX = xDim; lengthInY = yDim; }		// Standart-constructor 
	grid(size_t xDim, size_t yDim, T value){ data = std::vector<T>(xDim * yDim, value); lengthInX = xDim; lengthInY = yDim; }	// initalisation-constructor
	~grid(){}		// Destructor 

	size_t lengthX(){return lengthInX;}	// returns number of elements in x direction
	size_t lengthY(){ return lengthInX; }	// returns number of elements in y direction

	T& operator()(size_t i, size_t j){ 
		assert(i < lengthInX); 
		assert(j < lengthInY);
		return data[j*lengthInX + i]; }	// return a element at [j*nx+i] so i=x and j=y 

};


//___________________________________________________________________________________________________________________

// converts a timeval struct to seconds
double timevalToDouble( struct timeval *t ){
    return (double)(t->tv_sec) + (((double)(t->tv_usec))/1000000.0);
}

// calculates the l2-norm of the residual
double residuum(int nx, int ny, grid<type> &u, grid<type> &f, double h){
	double residuum = 0;
	double sum = 0;
	for (int i = 1; i<ny - 1; i++){
		for (int k = 1; k<nx - 1; k++){
			double temp = ((u(k - 1, i) + u(k + 1, i)) + (u(k, i - 1) + u(k, i + 1)) - (4 * u(k, i))) / (h*h) + f(k, i);
			sum += temp*temp;
		}
	}
	residuum = sqrt((1.0 / ((nx - 1)*(ny - 1)))*sum);
	return residuum;
}

// calculates the residual
void residual(int nx, int ny, grid<type> &u, grid<type> &f, grid<type> &res, double h){
	for (int i = 1; i<ny - 1; i++){
		for (int k = 1; k<nx - 1; k++){
			res(k, i) = ((u(k - 1, i) + u(k + 1, i)) + (u(k, i - 1) + u(k, i + 1)) - (4 * u(k, i))) / (h*h) + f(k, i);
		}
	}
}

// restrict a grid (from) with the full weighting stencile to another grid (to). From has size nx[l] and to has size nx[l-1].
void coarsening(int l, grid<type>& from, grid<type>& to, intVec& nx, intVec& ny){

	if (GHOST){
		for (int i = 1; i < ny[l - 1] - 1; i++){
			for (int j = 1; j < nx[l - 1] - 1; j++){
				to(j, i) = (from(2 * j-2, 2 * i + 1) + 2 * from(2 * j-1, 2 * i + 1) + from(2 * j, 2 * i + 1)
					+ 2 * from(2 * j-2, 2 * i) + 4 * from(2 * j-1, 2 * i) + 2 * from(2 * j, 2 * i)
					+ from(2 * j-2, 2 * i - 1) + 2 * from(2 * j-1, 2 * i - 1) + from(2 * j, 2 * i - 1)) / 16.0; 
			}
		}
	}
	else{
		for (int i = 1; i < ny[l - 1] - 1; i++){
			for (int j = 1; j < nx[l - 1] - 1; j++){
				to(j, i) = (from(2 * j - 1, 2 * i + 1) + 2 * from(2 * j, 2 * i + 1) + from(2 * j + 1, 2 * i + 1)
					+ 2 * from(2 * j - 1, 2 * i) + 4 * from(2 * j, 2 * i) + 2 * from(2 * j + 1, 2 * i)
					+ from(2 * j - 1, 2 * i - 1) + 2 * from(2 * j, 2 * i - 1) + from(2 * j + 1, 2 * i - 1)) / 16.0; 
			}
		}
	}
}

void interpolation(int l, grid<type>& from, grid<type>& to, intVec& nx, intVec& ny){
	if (GHOST){
		for (int i = 1; i < ny[l] - 1; i++){
			for (int j = 1; j < nx[l] - 1; j++){
				if (i % 2 == 0 && j % 2 != 0){
					// wert uebernehmen
					to(j, i) = from(j / 2 +1, i / 2);
				}
				else if (i % 2 != 0 && j % 2 == 0){
					// kreuz
					to(j, i) = (type)(from(j / 2+1, i / 2) + from(j / 2 + 1, (i / 2) + 1) + from((j / 2)+2 , i / 2) + from((j / 2)+2, (i / 2) + 1)) / 4.0;
				}
				else if (i % 2 == 0 && j % 2 == 0){
					// vertikal
					to(j, i) = (double)(from(j / 2+1, i / 2) + from((j / 2)+2, i / 2)) / 2.0;
				}
				else{
					// horizontal 
					to(j, i) = (double)(from(j / 2 +1, i / 2) + from(j / 2 +1, (i / 2) + 1)) / 2.0;
				}
			}
		}
	}
	else{
		for (int i = 1; i < ny[l] - 1; i++){
			for (int j = 1; j < nx[l] - 1; j++){
				if (i % 2 == 0 && j % 2 == 0){
					// wert uebernehmen
					to(j, i) = from(j / 2, i / 2);
				}
				else if (i % 2 != 0 && j % 2 != 0){
					// kreuz
					to(j, i) = (type)(from(j / 2, i / 2) + from(j / 2, (i / 2) + 1) + from((j / 2) + 1, i / 2) + from((j / 2) + 1, (i / 2) + 1)) / 4.0;
				}
				else if (i % 2 == 0 && j % 2 != 0){
					// vertikal
					to(j, i) = (double)(from(j / 2, i / 2) + from((j / 2) + 1, i / 2)) / 2.0;
				}
				else{
					// horizontal 
					to(j, i) = (double)(from(j / 2, i / 2) + from(j / 2, (i / 2) + 1)) / 2.0;
				}
			}
		}
	}
}

void Red_Black_Gauss(int nx, int ny, grid<type> &u, grid<type> &f, double h, int finest, int numIterations){

	for (int iterations = 0; iterations<numIterations; iterations++){
		for (int m = 1; m<ny - 1; m++){
			int q = 1;
			if (m % 2 == 0){
				q++;
			}
			for (; q<nx - 1; q = q + 2){
				u(q, m) = (1.0 / 4.0) * (h*h*f(q, m) + (u(q - 1, m) + u(q + 1, m) + u(q, m - 1) + u(q, m + 1)));
			}
		}
		for (int m = 1; m<ny - 1; m++){
			int q = 1;
			if (m % 2 != 0){
				q++;
			}
			for (; q<nx - 1; q = q + 2){
				u(q, m) = (1.0 / 4.0) * (h*h*f(q, m) + (u(q - 1, m) + u(q + 1, m) + u(q, m - 1) + u(q, m + 1)));
			}
		}
	}

    if (finest==1){
        if (GHOST){
            for (int i = 0; i<ny; i++){
                u(0, i) = u(1, i) - h;
                u(nx - 1, i) = u(nx - 2, i) - h;
            }
        }
    }

}

void multigrid(int l, std::vector<grid<type>>& u, std::vector<grid<type>>& f, intVec& nx,
    intVec& ny, std::vector<type>& h, std::vector<grid<type>>& res, int finest, int v1 = 2, int v2 = 1, int gamma = 1){

    //Presmoothing
    Red_Black_Gauss(nx[l], ny[l], u[l], f[l], h[l], finest, v1);

	// Residuum
	residual(nx[l], ny[l], u[l], f[l], res[l], h[l]);

	// restrict residual
	coarsening(l, res[l], f[l - 1], nx, ny);

	if (l <= 1){
		// solve
        Red_Black_Gauss(nx[l - 1], ny[l - 1], u[l - 1], f[l - 1], h[l - 1], finest, 1);

	}
	else{
		for (int i = 1; i<ny[l - 1] - 1; i++){
			for (int j = 1; j<nx[l - 1] - 1; j++){
				u[l - 1](j, i) = 0.0;
			}
		}
		for (int i = 0; i<gamma; i++){
            multigrid(l - 1, u, f, nx, ny, h, res, 0, v1, v2, gamma);
		}


		// interpolation
		grid<type> correction(nx[l], ny[l], 0.0);
		interpolation(l, u[l - 1], correction, nx, ny);

		// corretion
		for (int i = 1; i<ny[l] - 1; i++){
			for (int j = 1; j<nx[l] - 1; j++){
				u[l](j, i) += correction(j, i);
			}
		}
	}


	//Postsmothing
    Red_Black_Gauss(nx[l], ny[l], u[l], f[l], h[l], finest, v2);
}

int main(int argc, char **argv){

// Ueberpruefung, ob Eingabeparamter passen
    if (argc != 3){
        fprintf(stderr, "Usage: ./mgsolve levels numberOfVCycles\n");
        exit(EXIT_SUCCESS);
    }

	// definitions
    int l = atoi(argv[1]);	// number of levels
    int n = atoi(argv[2]);	// number of v-cycles

    intVec  nx(l, 0);	// total number of gird points in x-direction
	intVec  ny(l, 0);	// total number of grid points in y-direction
	typeVec h(l, 0.0);	// mesh size of each levels, where h[0] is the mesh size of the coarsesed grid
	type convergence = 0.0;	// convergence factor ( res_i+1 / res_i )
	type old_residuum = 0.0;	// l2-norm of the residuum from the previous iteration
	type new_residuum = 0.0;	// l2-norm of the residuum from the current iteration
	std::vector<grid<type>> f(l);	// vector of grids for the right hand side
	std::vector<grid<type>> res(l);	// vector of grids for the residuums
	std::vector<grid<type>> u(l);		// vector of grids for the approximation u

	// initialisation -------------------------------------------------------------------------------------------
	
	// vectors for mesh size and number of grid points
	for (int i = l - 1; i >= 0; i--){
		nx[i] = (int)(pow(2, i + 1) + 1);
		ny[i] = (int)(pow(2, i + 1) + 1);
		h[i] = 1.0 / (nx[i] - 1);
	}
	// Bei Neumann Randwerte werden in x-Richtung ghost-Layers benoetigt
	if (GHOST){
		for (int i = 0; i<l; i++)
			nx[i] = nx[i] + 2;
	}

	// Anlegen von f:
	for (int i = l - 1; i >= 0; i--){
		if (GHOST){
			f[i] = grid<type>(nx[i] , ny[i], 2.0);
		}
		else{
			f[i] = grid<type>(nx[i] , ny[i], 0.0);
		}
	}

	// Speicher fuer das residual allokieren
	for (int i = l - 1; i >= 0; i--){
		res[i] = grid<type>(nx[i] , ny[i], 0.0);
	}

	// Initialisierung des Gitters
	if (GHOST){
		for (int i = l - 1; i >= 0; i--){
			u[i] = grid<type>(nx[i] , ny[i], 0.0);
		}

		// Dirichlet RDB setzten aber nicht auf den Ghost-Layers 
		for (int i = 0; i<nx[l - 1] - 2; i++){
			u[l - 1](i + 1,	   		    0 ) = i*h[l - 1] * (1 - i*h[l - 1]);
			u[l - 1](i + 1, ny[l - 1] - 1 ) = i*h[l - 1] * (1 - i*h[l - 1]);
		}
		for (int i = 0; i<ny[l - 1]; i++){
			u[l - 1](0			  , i) = -h[l - 1];
			u[l - 1](nx[l - 1] - 1, i) = -h[l - 1];
		}

	}
	else{
		for (int i = l - 1; i >= 0; i--){
			u[i] = grid<type>(nx[i] , ny[i], 0.0);
		}
		for (int i = 0; i<nx[l - 1]; i++){
			u[l - 1](i, ny[l - 1] - 1) = sin(M_PI*i*h[l - 1])*sinh(M_PI);

		}
	}

    // TIMER
    struct timeval t1;
    gettimeofday(&t1, NULL);

	// Multigrid solver ------------------------------------------------------------------------------------------------
	for (int j = 0; j<n; j++){
        multigrid(l - 1, u, f, nx, ny, h, res, n, 1);

		// computing the convergence faktor in each cycle
		new_residuum = residuum(nx[l - 1], ny[l - 1], u[l - 1], f[l - 1], h[l - 1]);
		fprintf(stdout, "L2 Norm: %lf\n", new_residuum);
		if (j>0){
			convergence = new_residuum / old_residuum;
			fprintf(stdout, "Convergence: %lf\n", convergence);
		}
		old_residuum = new_residuum;
	}

    struct timeval t2;
    gettimeofday(&t2, NULL);
    fprintf(stdout, "Timer Multigrid: %lf\n", timevalToDouble(&t2)-timevalToDouble(&t1));

	// output --------------------------------------------------------------------------------------------------------------
	std::ofstream out;
	out.open("solution.txt", std::ios::out);

	double errorSum = 0.0;
	// Ausgabe fuer solution.txt

    if(GHOST){
        out << "# x y u(x,y)\n" << std::endl;
        for (int j = 0; j<ny[l - 1]; j++){
            for (int t = 0; t<nx[l - 1]-2; t++){
                    double x = (double)t / (double)(nx[l - 1] - 3);
                    double y = (double)j / (double)(ny[l - 1] - 1);
                    double temp = u[l - 1](t+1, j) - x*(1-x);
                    errorSum += temp * temp;

                    out <<  x << " " << y << " " << u[l - 1](t+1, j) << std::endl;
            }
            out << std::endl;
        }
    }else{
        out << "# x y u(x,y)\n" << std::endl;
        for (int j = 0; j<ny[l - 1]; j++){
            for (int t = 0; t<nx[l - 1]; t++){
                    double x = (double)t / (double)(nx[l - 1] - 1);
                    double y = (double)j / (double)(ny[l - 1] - 1);
                    double temp = u[l - 1](t, j) - sin(x*M_PI)*sinh(y*M_PI);
                    errorSum += temp * temp;

                    out <<  x << " " << y << " " << u[l - 1](t, j) << std::endl;
            }
            out << std::endl;
        }
    }

	out.close();
	

	double error = sqrt((1.0 / ((nx[l - 1] - 1)*(ny[l - 1] - 1)))*errorSum);
	fprintf(stdout, "Fehler zur korrekten Lsg: %lf\n", error);
	 
}
