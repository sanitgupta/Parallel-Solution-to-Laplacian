#include<iostream>
#include<fstream>
#include<execinfo.h>
#include<math.h>
#include<unistd.h>
#include<cstddef>
#include<metis.h>
#include<vector>
#include<algorithm>
#include<mpi.h>
#include<time.h>

using namespace std;

int main(int argc, char *argv[])
{
    int n = 64, iter = 1000, NT = 10, l = 0, minY, maxY;
    double s = 0, sk = 0, st1 = 0, st2 = 0, st3 = 0, st4 = 0;
    srand (static_cast <unsigned> (time(0)));
    clock_t t1, t2, t3, t4, t5, t6, t7, t8, t9, t10, t11, t12;

    double  *A;
    A = (double*)calloc(n*n*n, sizeof(double));

    A[n*n*int(n/2) + n*int(n/2) + int(n/2)] = 1000;

	int ierr = MPI_Init(&argc, &argv);
    int procid, numprocs;
    ierr = MPI_Comm_rank(MPI_COMM_WORLD, &procid);
    ierr = MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
	
	for (int T = 0; T < NT; T++) {

		if (procid == 0)
			t1 = clock();

		double  **XYZ;

        XYZ = (double**)calloc(n*n*n, sizeof(double *));  
        //temps=(double**)calloc(numprocs,sizeof(double *));  

        for(int i = 0; i < n*n*n; i++)
            XYZ[i] = (double*)calloc(3, sizeof(double));

        for(int i = 0; i < n; i++){
            for(int j = 0; j < n; j++){
                for(int k = 0; k < n; k++){
                XYZ[i*n*n+j*n+k][0] = i;
                XYZ[i*n*n+j*n+k][1] = j;
                XYZ[i*n*n+j*n+k][2] = k;
                }
            }
        }

       

        idx_t nVertices = n*n*n;
        idx_t nEdges    = ( 3*8 + 4*12*(n-2) + 5*6*(n-2)*(n-2) + 6*(n-2)*(n-2)*(n-2) )/ 2;//modify

        idx_t nWeights  = 1;
        idx_t nParts    = numprocs;

        idx_t objval;
        idx_t part[nVertices];
    
        idx_t xadj[nVertices + 1];

        idx_t options[METIS_NOPTIONS];      
        
        METIS_SetDefaultOptions(options);

		vector<int> a(nParts, 0);

		vector<vector<int> > v(nVertices, vector<int>(0));
	 
		vector<vector<vector<int> > > comm(nParts, vector<vector<int> >(nParts, vector<int>(0)));

		/* 
		Here, the graph of the square grid is formed,
		and stored in the CSR format so that it can be later used by METIS.
		*/

		xadj[0] = 0;

		for(int i = 0; i < n; i++){
            for(int j = 0; j < n; j++){
                for(int k = 0; k < n; k++){
                    if((i==0||i==n-1)&&(j==0||j==n-1)&&(k==0||k==n-1))
                        xadj[i*n*n+j*n+k+1]=xadj[i*n*n+j*n+k]+3;

                    else if(((i==0||i==n-1)&&(j==0||j==n-1))||((j==0||j==n-1)&&(k==0||k==n-1))||((k==0||k==n-1)&&(i==0||i==n-1)))     
                        xadj[i*n*n+j*n+k+1]=xadj[i*n*n+j*n+k]+4;

                    else if(i==0||i==n-1||j==0||j==n-1||k==0||k==n-1)     
                        xadj[i*n*n+j*n+k+1]=xadj[i*n*n+j*n+k]+5;
            
                    else
                        xadj[i*n*n+j*n+k+1]=xadj[i*n*n+j*n+k]+6;
                }
            }
        }

		idx_t adjncy[2 * nEdges];
		double adjnd[2 * nEdges];
        double adjnwgt[2 * nEdges];

        int ref;
        l=0;
        for(int i = 0; i < n; i++) {
            for(int j = 0; j < n; j++) {
                for(int k = 0; k < n; k++){

                    ref = i*n*n+j*n+k;
                    if(i==0 && j==0 && k==0)
                    {
                        adjncy[l] = (i)*n*n + (j)*n + (k+1);
                        adjnd[l] = pow(pow(XYZ[adjncy[l]][0]-XYZ[ref][0],2) + pow(XYZ[adjncy[l]][1]-XYZ[ref][1],2) + pow(XYZ[adjncy[l]][2]-XYZ[ref][2],2),0.5);
                        l++;

                        adjncy[l] = (i)*n*n + (j+1)*n + (k);
                        adjnd[l] = pow(pow(XYZ[adjncy[l]][0]-XYZ[ref][0],2) + pow(XYZ[adjncy[l]][1]-XYZ[ref][1],2) + pow(XYZ[adjncy[l]][2]-XYZ[ref][2],2),0.5);
                        l++;

                        adjncy[l] = (i+1)*n*n + (j)*n + (k);
                        adjnd[l] = pow(pow(XYZ[adjncy[l]][0]-XYZ[ref][0],2) + pow(XYZ[adjncy[l]][1]-XYZ[ref][1],2) + pow(XYZ[adjncy[l]][2]-XYZ[ref][2],2),0.5);
                        l++;
                    }

                    else if(i==0 && j==0 && k==n-1)
                    {
                        adjncy[l] = (i)*n*n + (j)*n + (k-1);
                        adjnd[l] = pow(pow(XYZ[adjncy[l]][0]-XYZ[ref][0],2) + pow(XYZ[adjncy[l]][1]-XYZ[ref][1],2) + pow(XYZ[adjncy[l]][2]-XYZ[ref][2],2),0.5);
                        l++;

                        adjncy[l] = (i)*n*n + (j+1)*n + (k);
                        adjnd[l] = pow(pow(XYZ[adjncy[l]][0]-XYZ[ref][0],2) + pow(XYZ[adjncy[l]][1]-XYZ[ref][1],2) + pow(XYZ[adjncy[l]][2]-XYZ[ref][2],2),0.5);
                        l++;

                        adjncy[l] = (i+1)*n*n + (j)*n + (k);
                        adjnd[l] = pow(pow(XYZ[adjncy[l]][0]-XYZ[ref][0],2) + pow(XYZ[adjncy[l]][1]-XYZ[ref][1],2) + pow(XYZ[adjncy[l]][2]-XYZ[ref][2],2),0.5);
                        l++;
                    }
                    
                    else if(i==0 && j==n-1 && k==0)
                    {
                        adjncy[l] = (i)*n*n + (j)*n + (k+1);
                        adjnd[l] = pow(pow(XYZ[adjncy[l]][0]-XYZ[ref][0],2) + pow(XYZ[adjncy[l]][1]-XYZ[ref][1],2) + pow(XYZ[adjncy[l]][2]-XYZ[ref][2],2),0.5);
                        l++;

                        adjncy[l] = (i)*n*n + (j-1)*n + (k);
                        adjnd[l] = pow(pow(XYZ[adjncy[l]][0]-XYZ[ref][0],2) + pow(XYZ[adjncy[l]][1]-XYZ[ref][1],2) + pow(XYZ[adjncy[l]][2]-XYZ[ref][2],2),0.5);
                        l++;

                        adjncy[l] = (i+1)*n*n + (j)*n + (k);
                        adjnd[l] = pow(pow(XYZ[adjncy[l]][0]-XYZ[ref][0],2) + pow(XYZ[adjncy[l]][1]-XYZ[ref][1],2) + pow(XYZ[adjncy[l]][2]-XYZ[ref][2],2),0.5);
                        l++;
                    }

                    else if(i==n-1 && j==0 && k==0)
                    {
                        adjncy[l] = (i)*n*n + (j)*n + (k+1);
                        adjnd[l] = pow(pow(XYZ[adjncy[l]][0]-XYZ[ref][0],2) + pow(XYZ[adjncy[l]][1]-XYZ[ref][1],2) + pow(XYZ[adjncy[l]][2]-XYZ[ref][2],2),0.5);
                        l++;

                        adjncy[l] = (i)*n*n + (j+1)*n + (k);
                        adjnd[l] = pow(pow(XYZ[adjncy[l]][0]-XYZ[ref][0],2) + pow(XYZ[adjncy[l]][1]-XYZ[ref][1],2) + pow(XYZ[adjncy[l]][2]-XYZ[ref][2],2),0.5);
                        l++;

                        adjncy[l] = (i-1)*n*n + (j)*n + (k);
                        adjnd[l] = pow(pow(XYZ[adjncy[l]][0]-XYZ[ref][0],2) + pow(XYZ[adjncy[l]][1]-XYZ[ref][1],2) + pow(XYZ[adjncy[l]][2]-XYZ[ref][2],2),0.5);
                        l++;
                    }

                    else if(i==0 && j==n-1 && k==n-1)
                    {
                        adjncy[l] = (i)*n*n + (j)*n + (k-1);
                        adjnd[l] = pow(pow(XYZ[adjncy[l]][0]-XYZ[ref][0],2) + pow(XYZ[adjncy[l]][1]-XYZ[ref][1],2) + pow(XYZ[adjncy[l]][2]-XYZ[ref][2],2),0.5);
                        l++;

                        adjncy[l] = (i)*n*n + (j-1)*n + (k);
                        adjnd[l] = pow(pow(XYZ[adjncy[l]][0]-XYZ[ref][0],2) + pow(XYZ[adjncy[l]][1]-XYZ[ref][1],2) + pow(XYZ[adjncy[l]][2]-XYZ[ref][2],2),0.5);
                        l++;

                        adjncy[l] = (i+1)*n*n + (j)*n + (k);
                        adjnd[l] = pow(pow(XYZ[adjncy[l]][0]-XYZ[ref][0],2) + pow(XYZ[adjncy[l]][1]-XYZ[ref][1],2) + pow(XYZ[adjncy[l]][2]-XYZ[ref][2],2),0.5);
                        l++;
                    }
                    else if(i==n-1 && j==0 && k==n-1)
                    {
                        adjncy[l] = (i)*n*n + (j)*n + (k-1);
                        adjnd[l] = pow(pow(XYZ[adjncy[l]][0]-XYZ[ref][0],2) + pow(XYZ[adjncy[l]][1]-XYZ[ref][1],2) + pow(XYZ[adjncy[l]][2]-XYZ[ref][2],2),0.5);
                        l++;

                        adjncy[l] = (i)*n*n + (j+1)*n + (k);
                        adjnd[l] = pow(pow(XYZ[adjncy[l]][0]-XYZ[ref][0],2) + pow(XYZ[adjncy[l]][1]-XYZ[ref][1],2) + pow(XYZ[adjncy[l]][2]-XYZ[ref][2],2),0.5);
                        l++;

                        adjncy[l] = (i+-1)*n*n + (j)*n + (k);
                        adjnd[l] = pow(pow(XYZ[adjncy[l]][0]-XYZ[ref][0],2) + pow(XYZ[adjncy[l]][1]-XYZ[ref][1],2) + pow(XYZ[adjncy[l]][2]-XYZ[ref][2],2),0.5);
                        l++;
                    }

                    else if(i==n-1 && j==n-1 && k==0)
                    {
                        adjncy[l] = (i)*n*n + (j)*n + (k+1);
                        adjnd[l] = pow(pow(XYZ[adjncy[l]][0]-XYZ[ref][0],2) + pow(XYZ[adjncy[l]][1]-XYZ[ref][1],2) + pow(XYZ[adjncy[l]][2]-XYZ[ref][2],2),0.5);
                        l++;

                        adjncy[l] = (i)*n*n + (j-1)*n + (k);
                        adjnd[l] = pow(pow(XYZ[adjncy[l]][0]-XYZ[ref][0],2) + pow(XYZ[adjncy[l]][1]-XYZ[ref][1],2) + pow(XYZ[adjncy[l]][2]-XYZ[ref][2],2),0.5);
                        l++;

                        adjncy[l] = (i-1)*n*n + (j)*n + (k);
                        adjnd[l] = pow(pow(XYZ[adjncy[l]][0]-XYZ[ref][0],2) + pow(XYZ[adjncy[l]][1]-XYZ[ref][1],2) + pow(XYZ[adjncy[l]][2]-XYZ[ref][2],2),0.5);
                        l++;
                    }

                    else if(i==n-1 && j==n-1 && k==n-1)
                    {
                        adjncy[l] = (i)*n*n + (j)*n + (k-1);
                        adjnd[l] = pow(pow(XYZ[adjncy[l]][0]-XYZ[ref][0],2) + pow(XYZ[adjncy[l]][1]-XYZ[ref][1],2) + pow(XYZ[adjncy[l]][2]-XYZ[ref][2],2),0.5);
                        l++;

                        adjncy[l] = (i)*n*n + (j-1)*n + (k);
                        adjnd[l] = pow(pow(XYZ[adjncy[l]][0]-XYZ[ref][0],2) + pow(XYZ[adjncy[l]][1]-XYZ[ref][1],2) + pow(XYZ[adjncy[l]][2]-XYZ[ref][2],2),0.5);
                        l++;

                        adjncy[l] = (i-1)*n*n + (j)*n + (k);
                        adjnd[l] = pow(pow(XYZ[adjncy[l]][0]-XYZ[ref][0],2) + pow(XYZ[adjncy[l]][1]-XYZ[ref][1],2) + pow(XYZ[adjncy[l]][2]-XYZ[ref][2],2),0.5);
                        l++;
                    }

                    else if(i==0 && j==0)
                    {
                        adjncy[l] = (i)*n*n + (j)*n + (k-1);
                        adjnd[l] = pow(pow(XYZ[adjncy[l]][0]-XYZ[ref][0],2) + pow(XYZ[adjncy[l]][1]-XYZ[ref][1],2) + pow(XYZ[adjncy[l]][2]-XYZ[ref][2],2),0.5);
                        l++;

                        adjncy[l] = (i)*n*n + (j+1)*n + (k);
                        adjnd[l] = pow(pow(XYZ[adjncy[l]][0]-XYZ[ref][0],2) + pow(XYZ[adjncy[l]][1]-XYZ[ref][1],2) + pow(XYZ[adjncy[l]][2]-XYZ[ref][2],2),0.5);
                        l++;

                        adjncy[l] = (i+1)*n*n + (j)*n + (k);
                        adjnd[l] = pow(pow(XYZ[adjncy[l]][0]-XYZ[ref][0],2) + pow(XYZ[adjncy[l]][1]-XYZ[ref][1],2) + pow(XYZ[adjncy[l]][2]-XYZ[ref][2],2),0.5);
                        l++;

                        adjncy[l] = (i)*n*n + (j)*n + (k+1);
                        adjnd[l] = pow(pow(XYZ[adjncy[l]][0]-XYZ[ref][0],2) + pow(XYZ[adjncy[l]][1]-XYZ[ref][1],2) + pow(XYZ[adjncy[l]][2]-XYZ[ref][2],2),0.5);
                        l++;
                    }

                    else if(i==0 && j==n-1)
                    {
                        adjncy[l] = (i)*n*n + (j)*n + (k-1);
                        adjnd[l] = pow(pow(XYZ[adjncy[l]][0]-XYZ[ref][0],2) + pow(XYZ[adjncy[l]][1]-XYZ[ref][1],2) + pow(XYZ[adjncy[l]][2]-XYZ[ref][2],2),0.5);
                        l++;

                        adjncy[l] = (i)*n*n + (j-1)*n + (k);
                        adjnd[l] = pow(pow(XYZ[adjncy[l]][0]-XYZ[ref][0],2) + pow(XYZ[adjncy[l]][1]-XYZ[ref][1],2) + pow(XYZ[adjncy[l]][2]-XYZ[ref][2],2),0.5);
                        l++;

                        adjncy[l] = (i+1)*n*n + (j)*n + (k);
                        adjnd[l] = pow(pow(XYZ[adjncy[l]][0]-XYZ[ref][0],2) + pow(XYZ[adjncy[l]][1]-XYZ[ref][1],2) + pow(XYZ[adjncy[l]][2]-XYZ[ref][2],2),0.5);
                        l++;

                        adjncy[l] = (i)*n*n + (j)*n + (k+1);
                        adjnd[l] = pow(pow(XYZ[adjncy[l]][0]-XYZ[ref][0],2) + pow(XYZ[adjncy[l]][1]-XYZ[ref][1],2) + pow(XYZ[adjncy[l]][2]-XYZ[ref][2],2),0.5);
                        l++;
                    }

                    else if(i==n-1 && j==0)
                    {
                        adjncy[l] = (i)*n*n + (j)*n + (k-1);
                        adjnd[l] = pow(pow(XYZ[adjncy[l]][0]-XYZ[ref][0],2) + pow(XYZ[adjncy[l]][1]-XYZ[ref][1],2) + pow(XYZ[adjncy[l]][2]-XYZ[ref][2],2),0.5);
                        l++;

                        adjncy[l] = (i)*n*n + (j+1)*n + (k);
                        adjnd[l] = pow(pow(XYZ[adjncy[l]][0]-XYZ[ref][0],2) + pow(XYZ[adjncy[l]][1]-XYZ[ref][1],2) + pow(XYZ[adjncy[l]][2]-XYZ[ref][2],2),0.5);
                        l++;

                        adjncy[l] = (i-1)*n*n + (j)*n + (k);
                        adjnd[l] = pow(pow(XYZ[adjncy[l]][0]-XYZ[ref][0],2) + pow(XYZ[adjncy[l]][1]-XYZ[ref][1],2) + pow(XYZ[adjncy[l]][2]-XYZ[ref][2],2),0.5);
                        l++;

                        adjncy[l] = (i)*n*n + (j)*n + (k+1);
                        adjnd[l] = pow(pow(XYZ[adjncy[l]][0]-XYZ[ref][0],2) + pow(XYZ[adjncy[l]][1]-XYZ[ref][1],2) + pow(XYZ[adjncy[l]][2]-XYZ[ref][2],2),0.5);
                        l++;
                    }

                    else if(i==n-1 && j==n-1)
                    {
                        adjncy[l] = (i)*n*n + (j)*n + (k-1);
                        adjnd[l] = pow(pow(XYZ[adjncy[l]][0]-XYZ[ref][0],2) + pow(XYZ[adjncy[l]][1]-XYZ[ref][1],2) + pow(XYZ[adjncy[l]][2]-XYZ[ref][2],2),0.5);
                        l++;

                        adjncy[l] = (i)*n*n + (j-1)*n + (k);
                        adjnd[l] = pow(pow(XYZ[adjncy[l]][0]-XYZ[ref][0],2) + pow(XYZ[adjncy[l]][1]-XYZ[ref][1],2) + pow(XYZ[adjncy[l]][2]-XYZ[ref][2],2),0.5);
                        l++;

                        adjncy[l] = (i-1)*n*n + (j)*n + (k);
                        adjnd[l] = pow(pow(XYZ[adjncy[l]][0]-XYZ[ref][0],2) + pow(XYZ[adjncy[l]][1]-XYZ[ref][1],2) + pow(XYZ[adjncy[l]][2]-XYZ[ref][2],2),0.5);
                        l++;

                        adjncy[l] = (i)*n*n + (j)*n + (k+1);
                        adjnd[l] = pow(pow(XYZ[adjncy[l]][0]-XYZ[ref][0],2) + pow(XYZ[adjncy[l]][1]-XYZ[ref][1],2) + pow(XYZ[adjncy[l]][2]-XYZ[ref][2],2),0.5);
                        l++;
                    }


                    else if(k==0 && j==0)
                    {
                        adjncy[l] = (i)*n*n + (j)*n + (k+1);
                        adjnd[l] = pow(pow(XYZ[adjncy[l]][0]-XYZ[ref][0],2) + pow(XYZ[adjncy[l]][1]-XYZ[ref][1],2) + pow(XYZ[adjncy[l]][2]-XYZ[ref][2],2),0.5);
                        l++;

                        adjncy[l] = (i)*n*n + (j+1)*n + (k);
                        adjnd[l] = pow(pow(XYZ[adjncy[l]][0]-XYZ[ref][0],2) + pow(XYZ[adjncy[l]][1]-XYZ[ref][1],2) + pow(XYZ[adjncy[l]][2]-XYZ[ref][2],2),0.5);
                        l++;

                        adjncy[l] = (i+1)*n*n + (j)*n + (k);
                        adjnd[l] = pow(pow(XYZ[adjncy[l]][0]-XYZ[ref][0],2) + pow(XYZ[adjncy[l]][1]-XYZ[ref][1],2) + pow(XYZ[adjncy[l]][2]-XYZ[ref][2],2),0.5);
                        l++;

                        adjncy[l] = (i-1)*n*n + (j)*n + (k);
                        adjnd[l] = pow(pow(XYZ[adjncy[l]][0]-XYZ[ref][0],2) + pow(XYZ[adjncy[l]][1]-XYZ[ref][1],2) + pow(XYZ[adjncy[l]][2]-XYZ[ref][2],2),0.5);
                        l++;
                    }

                    else if(k==0 && j==n-1)
                    {
                        adjncy[l] = (i)*n*n + (j)*n + (k+1);
                        adjnd[l] = pow(pow(XYZ[adjncy[l]][0]-XYZ[ref][0],2) + pow(XYZ[adjncy[l]][1]-XYZ[ref][1],2) + pow(XYZ[adjncy[l]][2]-XYZ[ref][2],2),0.5);
                        l++;

                        adjncy[l] = (i)*n*n + (j-1)*n + (k);
                        adjnd[l] = pow(pow(XYZ[adjncy[l]][0]-XYZ[ref][0],2) + pow(XYZ[adjncy[l]][1]-XYZ[ref][1],2) + pow(XYZ[adjncy[l]][2]-XYZ[ref][2],2),0.5);
                        l++;

                        adjncy[l] = (i+1)*n*n + (j)*n + (k);
                        adjnd[l] = pow(pow(XYZ[adjncy[l]][0]-XYZ[ref][0],2) + pow(XYZ[adjncy[l]][1]-XYZ[ref][1],2) + pow(XYZ[adjncy[l]][2]-XYZ[ref][2],2),0.5);
                        l++;

                        adjncy[l] = (i-1)*n*n + (j)*n + (k);
                        adjnd[l] = pow(pow(XYZ[adjncy[l]][0]-XYZ[ref][0],2) + pow(XYZ[adjncy[l]][1]-XYZ[ref][1],2) + pow(XYZ[adjncy[l]][2]-XYZ[ref][2],2),0.5);
                        l++;
                    }

                    else if(k==n-1 && j==0)
                    {
                        adjncy[l] = (i)*n*n + (j)*n + (k-1);
                        adjnd[l] = pow(pow(XYZ[adjncy[l]][0]-XYZ[ref][0],2) + pow(XYZ[adjncy[l]][1]-XYZ[ref][1],2) + pow(XYZ[adjncy[l]][2]-XYZ[ref][2],2),0.5);
                        l++;

                        adjncy[l] = (i)*n*n + (j+1)*n + (k);
                        adjnd[l] = pow(pow(XYZ[adjncy[l]][0]-XYZ[ref][0],2) + pow(XYZ[adjncy[l]][1]-XYZ[ref][1],2) + pow(XYZ[adjncy[l]][2]-XYZ[ref][2],2),0.5);
                        l++;

                        adjncy[l] = (i+1)*n*n + (j)*n + (k);
                        adjnd[l] = pow(pow(XYZ[adjncy[l]][0]-XYZ[ref][0],2) + pow(XYZ[adjncy[l]][1]-XYZ[ref][1],2) + pow(XYZ[adjncy[l]][2]-XYZ[ref][2],2),0.5);
                        l++;

                        adjncy[l] = (i-1)*n*n + (j)*n + (k+1);
                        adjnd[l] = pow(pow(XYZ[adjncy[l]][0]-XYZ[ref][0],2) + pow(XYZ[adjncy[l]][1]-XYZ[ref][1],2) + pow(XYZ[adjncy[l]][2]-XYZ[ref][2],2),0.5);
                        l++;
                    }

                    else if(k==n-1 && j==n-1)
                    {
                        adjncy[l] = (i)*n*n + (j)*n + (k-1);
                        adjnd[l] = pow(pow(XYZ[adjncy[l]][0]-XYZ[ref][0],2) + pow(XYZ[adjncy[l]][1]-XYZ[ref][1],2) + pow(XYZ[adjncy[l]][2]-XYZ[ref][2],2),0.5);
                        l++;

                        adjncy[l] = (i)*n*n + (j-1)*n + (k);
                        adjnd[l] = pow(pow(XYZ[adjncy[l]][0]-XYZ[ref][0],2) + pow(XYZ[adjncy[l]][1]-XYZ[ref][1],2) + pow(XYZ[adjncy[l]][2]-XYZ[ref][2],2),0.5);
                        l++;

                        adjncy[l] = (i+1)*n*n + (j)*n + (k);
                        adjnd[l] = pow(pow(XYZ[adjncy[l]][0]-XYZ[ref][0],2) + pow(XYZ[adjncy[l]][1]-XYZ[ref][1],2) + pow(XYZ[adjncy[l]][2]-XYZ[ref][2],2),0.5);
                        l++;

                        adjncy[l] = (i-1)*n*n + (j)*n + (k+1);
                        adjnd[l] = pow(pow(XYZ[adjncy[l]][0]-XYZ[ref][0],2) + pow(XYZ[adjncy[l]][1]-XYZ[ref][1],2) + pow(XYZ[adjncy[l]][2]-XYZ[ref][2],2),0.5);
                        l++;
                    }


                    else if(i==0 && k==0)
                    {
                        adjncy[l] = (i)*n*n + (j)*n + (k+1);
                        adjnd[l] = pow(pow(XYZ[adjncy[l]][0]-XYZ[ref][0],2) + pow(XYZ[adjncy[l]][1]-XYZ[ref][1],2) + pow(XYZ[adjncy[l]][2]-XYZ[ref][2],2),0.5);
                        l++;

                        adjncy[l] = (i)*n*n + (j+1)*n + (k);
                        adjnd[l] = pow(pow(XYZ[adjncy[l]][0]-XYZ[ref][0],2) + pow(XYZ[adjncy[l]][1]-XYZ[ref][1],2) + pow(XYZ[adjncy[l]][2]-XYZ[ref][2],2),0.5);
                        l++;

                        adjncy[l] = (i)*n*n + (j-1)*n + (k);
                        adjnd[l] = pow(pow(XYZ[adjncy[l]][0]-XYZ[ref][0],2) + pow(XYZ[adjncy[l]][1]-XYZ[ref][1],2) + pow(XYZ[adjncy[l]][2]-XYZ[ref][2],2),0.5);
                        l++;

                        adjncy[l] = (i+1)*n*n + (j)*n + (k);
                        adjnd[l] = pow(pow(XYZ[adjncy[l]][0]-XYZ[ref][0],2) + pow(XYZ[adjncy[l]][1]-XYZ[ref][1],2) + pow(XYZ[adjncy[l]][2]-XYZ[ref][2],2),0.5);
                        l++;
                    }

                    else if(i==0 && k==n-1)
                    {
                        adjncy[l] = (i)*n*n + (j)*n + (k-1);
                        adjnd[l] = pow(pow(XYZ[adjncy[l]][0]-XYZ[ref][0],2) + pow(XYZ[adjncy[l]][1]-XYZ[ref][1],2) + pow(XYZ[adjncy[l]][2]-XYZ[ref][2],2),0.5);
                        l++;

                        adjncy[l] = (i)*n*n + (j-1)*n + (k);
                        adjnd[l] = pow(pow(XYZ[adjncy[l]][0]-XYZ[ref][0],2) + pow(XYZ[adjncy[l]][1]-XYZ[ref][1],2) + pow(XYZ[adjncy[l]][2]-XYZ[ref][2],2),0.5);
                        l++;

                        adjncy[l] = (i)*n*n + (j+1)*n + (k);
                        adjnd[l] = pow(pow(XYZ[adjncy[l]][0]-XYZ[ref][0],2) + pow(XYZ[adjncy[l]][1]-XYZ[ref][1],2) + pow(XYZ[adjncy[l]][2]-XYZ[ref][2],2),0.5);
                        l++;

                        adjncy[l] = (i+1)*n*n + (j)*n + (k);
                        adjnd[l] = pow(pow(XYZ[adjncy[l]][0]-XYZ[ref][0],2) + pow(XYZ[adjncy[l]][1]-XYZ[ref][1],2) + pow(XYZ[adjncy[l]][2]-XYZ[ref][2],2),0.5);
                        l++;
                    }

                    else if(i==n-1 && k==0)
                    {
                        adjncy[l] = (i)*n*n + (j)*n + (k+1);
                        adjnd[l] = pow(pow(XYZ[adjncy[l]][0]-XYZ[ref][0],2) + pow(XYZ[adjncy[l]][1]-XYZ[ref][1],2) + pow(XYZ[adjncy[l]][2]-XYZ[ref][2],2),0.5);
                        l++;

                        adjncy[l] = (i)*n*n + (j+1)*n + (k);
                        adjnd[l] = pow(pow(XYZ[adjncy[l]][0]-XYZ[ref][0],2) + pow(XYZ[adjncy[l]][1]-XYZ[ref][1],2) + pow(XYZ[adjncy[l]][2]-XYZ[ref][2],2),0.5);
                        l++;

                        adjncy[l] = (i)*n*n + (j-1)*n + (k);
                        adjnd[l] = pow(pow(XYZ[adjncy[l]][0]-XYZ[ref][0],2) + pow(XYZ[adjncy[l]][1]-XYZ[ref][1],2) + pow(XYZ[adjncy[l]][2]-XYZ[ref][2],2),0.5);
                        l++;

                        adjncy[l] = (i-1)*n*n + (j)*n + (k);
                        adjnd[l] = pow(pow(XYZ[adjncy[l]][0]-XYZ[ref][0],2) + pow(XYZ[adjncy[l]][1]-XYZ[ref][1],2) + pow(XYZ[adjncy[l]][2]-XYZ[ref][2],2),0.5);
                        l++;
                    }

                    else if(i==n-1 && k==n-1)
                    {
                        adjncy[l] = (i)*n*n + (j)*n + (k-1);
                        adjnd[l] = pow(pow(XYZ[adjncy[l]][0]-XYZ[ref][0],2) + pow(XYZ[adjncy[l]][1]-XYZ[ref][1],2) + pow(XYZ[adjncy[l]][2]-XYZ[ref][2],2),0.5);
                        l++;

                        adjncy[l] = (i)*n*n + (j-1)*n + (k);
                        adjnd[l] = pow(pow(XYZ[adjncy[l]][0]-XYZ[ref][0],2) + pow(XYZ[adjncy[l]][1]-XYZ[ref][1],2) + pow(XYZ[adjncy[l]][2]-XYZ[ref][2],2),0.5);
                        l++;

                        adjncy[l] = (i)*n*n + (j+1)*n + (k);
                        adjnd[l] = pow(pow(XYZ[adjncy[l]][0]-XYZ[ref][0],2) + pow(XYZ[adjncy[l]][1]-XYZ[ref][1],2) + pow(XYZ[adjncy[l]][2]-XYZ[ref][2],2),0.5);
                        l++;

                        adjncy[l] = (i-1)*n*n + (j)*n + (k);
                        adjnd[l] = pow(pow(XYZ[adjncy[l]][0]-XYZ[ref][0],2) + pow(XYZ[adjncy[l]][1]-XYZ[ref][1],2) + pow(XYZ[adjncy[l]][2]-XYZ[ref][2],2),0.5);
                        l++;
                    }

                    

                    else if(i==0)
                    {
                        adjncy[l] = (i)*n*n + (j)*n + (k-1);
                        adjnd[l] = pow(pow(XYZ[adjncy[l]][0]-XYZ[ref][0],2) + pow(XYZ[adjncy[l]][1]-XYZ[ref][1],2) + pow(XYZ[adjncy[l]][2]-XYZ[ref][2],2),0.5);
                        l++;

                        adjncy[l] = (i)*n*n + (j)*n + (k+1);
                        adjnd[l] = pow(pow(XYZ[adjncy[l]][0]-XYZ[ref][0],2) + pow(XYZ[adjncy[l]][1]-XYZ[ref][1],2) + pow(XYZ[adjncy[l]][2]-XYZ[ref][2],2),0.5);
                        l++;

                        adjncy[l] = (i)*n*n + (j-1)*n + (k);
                        adjnd[l] = pow(pow(XYZ[adjncy[l]][0]-XYZ[ref][0],2) + pow(XYZ[adjncy[l]][1]-XYZ[ref][1],2) + pow(XYZ[adjncy[l]][2]-XYZ[ref][2],2),0.5);
                        l++;

                        adjncy[l] = (i)*n*n + (j+1)*n + (k);
                        adjnd[l] = pow(pow(XYZ[adjncy[l]][0]-XYZ[ref][0],2) + pow(XYZ[adjncy[l]][1]-XYZ[ref][1],2) + pow(XYZ[adjncy[l]][2]-XYZ[ref][2],2),0.5);
                        l++;


                        adjncy[l] = (i+1)*n*n + (j)*n + (k);
                        adjnd[l] = pow(pow(XYZ[adjncy[l]][0]-XYZ[ref][0],2) + pow(XYZ[adjncy[l]][1]-XYZ[ref][1],2) + pow(XYZ[adjncy[l]][2]-XYZ[ref][2],2),0.5);
                        l++;
                    }

                    else if(i==n-1)
                    {
                        adjncy[l] = (i)*n*n + (j)*n + (k-1);
                        adjnd[l] = pow(pow(XYZ[adjncy[l]][0]-XYZ[ref][0],2) + pow(XYZ[adjncy[l]][1]-XYZ[ref][1],2) + pow(XYZ[adjncy[l]][2]-XYZ[ref][2],2),0.5);
                        l++;

                        adjncy[l] = (i)*n*n + (j)*n + (k+1);
                        adjnd[l] = pow(pow(XYZ[adjncy[l]][0]-XYZ[ref][0],2) + pow(XYZ[adjncy[l]][1]-XYZ[ref][1],2) + pow(XYZ[adjncy[l]][2]-XYZ[ref][2],2),0.5);
                        l++;

                        adjncy[l] = (i)*n*n + (j-1)*n + (k);
                        adjnd[l] = pow(pow(XYZ[adjncy[l]][0]-XYZ[ref][0],2) + pow(XYZ[adjncy[l]][1]-XYZ[ref][1],2) + pow(XYZ[adjncy[l]][2]-XYZ[ref][2],2),0.5);
                        l++;

                        adjncy[l] = (i)*n*n + (j+1)*n + (k);
                        adjnd[l] = pow(pow(XYZ[adjncy[l]][0]-XYZ[ref][0],2) + pow(XYZ[adjncy[l]][1]-XYZ[ref][1],2) + pow(XYZ[adjncy[l]][2]-XYZ[ref][2],2),0.5);
                        l++;

                        adjncy[l] = (i-1)*n*n + (j)*n + (k);
                        adjnd[l] = pow(pow(XYZ[adjncy[l]][0]-XYZ[ref][0],2) + pow(XYZ[adjncy[l]][1]-XYZ[ref][1],2) + pow(XYZ[adjncy[l]][2]-XYZ[ref][2],2),0.5);
                        l++;
                    }

                    else if(j==0)
                    {
                        adjncy[l] = (i)*n*n + (j)*n + (k-1);
                        adjnd[l] = pow(pow(XYZ[adjncy[l]][0]-XYZ[ref][0],2) + pow(XYZ[adjncy[l]][1]-XYZ[ref][1],2) + pow(XYZ[adjncy[l]][2]-XYZ[ref][2],2),0.5);
                        l++;

                        adjncy[l] = (i)*n*n + (j)*n + (k+1);
                        adjnd[l] = pow(pow(XYZ[adjncy[l]][0]-XYZ[ref][0],2) + pow(XYZ[adjncy[l]][1]-XYZ[ref][1],2) + pow(XYZ[adjncy[l]][2]-XYZ[ref][2],2),0.5);
                        l++;

                        adjncy[l] = (i)*n*n + (j+1)*n + (k);
                        adjnd[l] = pow(pow(XYZ[adjncy[l]][0]-XYZ[ref][0],2) + pow(XYZ[adjncy[l]][1]-XYZ[ref][1],2) + pow(XYZ[adjncy[l]][2]-XYZ[ref][2],2),0.5);
                        l++;

                        adjncy[l] = (i-1)*n*n + (j)*n + (k);
                        adjnd[l] = pow(pow(XYZ[adjncy[l]][0]-XYZ[ref][0],2) + pow(XYZ[adjncy[l]][1]-XYZ[ref][1],2) + pow(XYZ[adjncy[l]][2]-XYZ[ref][2],2),0.5);
                        l++;

                        adjncy[l] = (i+1)*n*n + (j)*n + (k);
                        adjnd[l] = pow(pow(XYZ[adjncy[l]][0]-XYZ[ref][0],2) + pow(XYZ[adjncy[l]][1]-XYZ[ref][1],2) + pow(XYZ[adjncy[l]][2]-XYZ[ref][2],2),0.5);
                        l++;
                    }

                    else if(j==n-1)
                    {
                        adjncy[l] = (i)*n*n + (j)*n + (k-1);
                        adjnd[l] = pow(pow(XYZ[adjncy[l]][0]-XYZ[ref][0],2) + pow(XYZ[adjncy[l]][1]-XYZ[ref][1],2) + pow(XYZ[adjncy[l]][2]-XYZ[ref][2],2),0.5);
                        l++;

                        adjncy[l] = (i)*n*n + (j)*n + (k+1);
                        adjnd[l] = pow(pow(XYZ[adjncy[l]][0]-XYZ[ref][0],2) + pow(XYZ[adjncy[l]][1]-XYZ[ref][1],2) + pow(XYZ[adjncy[l]][2]-XYZ[ref][2],2),0.5);
                        l++;

                        adjncy[l] = (i)*n*n + (j-1)*n + (k);
                        adjnd[l] = pow(pow(XYZ[adjncy[l]][0]-XYZ[ref][0],2) + pow(XYZ[adjncy[l]][1]-XYZ[ref][1],2) + pow(XYZ[adjncy[l]][2]-XYZ[ref][2],2),0.5);
                        l++;

                        adjncy[l] = (i-1)*n*n + (j)*n + (k);
                        adjnd[l] = pow(pow(XYZ[adjncy[l]][0]-XYZ[ref][0],2) + pow(XYZ[adjncy[l]][1]-XYZ[ref][1],2) + pow(XYZ[adjncy[l]][2]-XYZ[ref][2],2),0.5);
                        l++;

                        adjncy[l] = (i+1)*n*n + (j)*n + (k);
                        adjnd[l] = pow(pow(XYZ[adjncy[l]][0]-XYZ[ref][0],2) + pow(XYZ[adjncy[l]][1]-XYZ[ref][1],2) + pow(XYZ[adjncy[l]][2]-XYZ[ref][2],2),0.5);
                        l++;
                    }

                    else if(k==0)
                    {
                        adjncy[l] = (i)*n*n + (j)*n + (k+1);
                        adjnd[l] = pow(pow(XYZ[adjncy[l]][0]-XYZ[ref][0],2) + pow(XYZ[adjncy[l]][1]-XYZ[ref][1],2) + pow(XYZ[adjncy[l]][2]-XYZ[ref][2],2),0.5);
                        l++;

                        adjncy[l] = (i)*n*n + (j-1)*n + (k);
                        adjnd[l] = pow(pow(XYZ[adjncy[l]][0]-XYZ[ref][0],2) + pow(XYZ[adjncy[l]][1]-XYZ[ref][1],2) + pow(XYZ[adjncy[l]][2]-XYZ[ref][2],2),0.5);
                        l++;

                        adjncy[l] = (i)*n*n + (j+1)*n + (k);
                        adjnd[l] = pow(pow(XYZ[adjncy[l]][0]-XYZ[ref][0],2) + pow(XYZ[adjncy[l]][1]-XYZ[ref][1],2) + pow(XYZ[adjncy[l]][2]-XYZ[ref][2],2),0.5);
                        l++;

                        adjncy[l] = (i-1)*n*n + (j)*n + (k);
                        adjnd[l] = pow(pow(XYZ[adjncy[l]][0]-XYZ[ref][0],2) + pow(XYZ[adjncy[l]][1]-XYZ[ref][1],2) + pow(XYZ[adjncy[l]][2]-XYZ[ref][2],2),0.5);
                        l++;

                        adjncy[l] = (i+1)*n*n + (j)*n + (k);
                        adjnd[l] = pow(pow(XYZ[adjncy[l]][0]-XYZ[ref][0],2) + pow(XYZ[adjncy[l]][1]-XYZ[ref][1],2) + pow(XYZ[adjncy[l]][2]-XYZ[ref][2],2),0.5);
                        l++;
                    }

                    else if(k==n-1)
                    {
                        adjncy[l] = (i)*n*n + (j)*n + (k-1);
                        adjnd[l] = pow(pow(XYZ[adjncy[l]][0]-XYZ[ref][0],2) + pow(XYZ[adjncy[l]][1]-XYZ[ref][1],2) + pow(XYZ[adjncy[l]][2]-XYZ[ref][2],2),0.5);
                        l++;

                        adjncy[l] = (i)*n*n + (j-1)*n + (k);
                        adjnd[l] = pow(pow(XYZ[adjncy[l]][0]-XYZ[ref][0],2) + pow(XYZ[adjncy[l]][1]-XYZ[ref][1],2) + pow(XYZ[adjncy[l]][2]-XYZ[ref][2],2),0.5);
                        l++;

                        adjncy[l] = (i)*n*n + (j+1)*n + (k);
                        adjnd[l] = pow(pow(XYZ[adjncy[l]][0]-XYZ[ref][0],2) + pow(XYZ[adjncy[l]][1]-XYZ[ref][1],2) + pow(XYZ[adjncy[l]][2]-XYZ[ref][2],2),0.5);
                        l++;

                        adjncy[l] = (i-1)*n*n + (j)*n + (k);
                        adjnd[l] = pow(pow(XYZ[adjncy[l]][0]-XYZ[ref][0],2) + pow(XYZ[adjncy[l]][1]-XYZ[ref][1],2) + pow(XYZ[adjncy[l]][2]-XYZ[ref][2],2),0.5);
                        l++;

                        adjncy[l] = (i+1)*n*n + (j)*n + (k);
                        adjnd[l] = pow(pow(XYZ[adjncy[l]][0]-XYZ[ref][0],2) + pow(XYZ[adjncy[l]][1]-XYZ[ref][1],2) + pow(XYZ[adjncy[l]][2]-XYZ[ref][2],2),0.5);
                        l++;
                    }

                    else 
                    {
                        adjncy[l] = (i)*n*n + (j)*n + (k-1);
                        adjnd[l] = pow(pow(XYZ[adjncy[l]][0]-XYZ[ref][0],2) + pow(XYZ[adjncy[l]][1]-XYZ[ref][1],2) + pow(XYZ[adjncy[l]][2]-XYZ[ref][2],2),0.5);
                        l++;

                        adjncy[l] = (i)*n*n + (j)*n + (k+1);
                        adjnd[l] = pow(pow(XYZ[adjncy[l]][0]-XYZ[ref][0],2) + pow(XYZ[adjncy[l]][1]-XYZ[ref][1],2) + pow(XYZ[adjncy[l]][2]-XYZ[ref][2],2),0.5);
                        l++;

                        adjncy[l] = (i)*n*n + (j-1)*n + (k);
                        adjnd[l] = pow(pow(XYZ[adjncy[l]][0]-XYZ[ref][0],2) + pow(XYZ[adjncy[l]][1]-XYZ[ref][1],2) + pow(XYZ[adjncy[l]][2]-XYZ[ref][2],2),0.5);
                        l++;

                        adjncy[l] = (i)*n*n + (j+1)*n + (k);
                        adjnd[l] = pow(pow(XYZ[adjncy[l]][0]-XYZ[ref][0],2) + pow(XYZ[adjncy[l]][1]-XYZ[ref][1],2) + pow(XYZ[adjncy[l]][2]-XYZ[ref][2],2),0.5);
                        l++;

                        adjncy[l] = (i-1)*n*n + (j)*n + (k);
                        adjnd[l] = pow(pow(XYZ[adjncy[l]][0]-XYZ[ref][0],2) + pow(XYZ[adjncy[l]][1]-XYZ[ref][1],2) + pow(XYZ[adjncy[l]][2]-XYZ[ref][2],2),0.5);
                        l++;

                        adjncy[l] = (i+1)*n*n + (j)*n + (k);
                        adjnd[l] = pow(pow(XYZ[adjncy[l]][0]-XYZ[ref][0],2) + pow(XYZ[adjncy[l]][1]-XYZ[ref][1],2) + pow(XYZ[adjncy[l]][2]-XYZ[ref][2],2),0.5);
                        l++;
                    }
                }       
            }
        }

        for (int i = 0; i < nVertices; i++) {
            double twgt = 0;
            for (int j = xadj[i]; j < xadj[i+1]; j++) {
                twgt += adjnd[j];
            }
            for (int j = xadj[i]; j < xadj[i+1]; j++) {
                adjnwgt[j] = adjnd[j]/twgt;
            }
        }

        // Weights of vertices
        // if all weights are equal then can be set to NULL
        idx_t vwgt[nVertices * nWeights];


        /*int METIS PartGraphKway(idx t *nvtxs, idx t *ncon, idx t *xadj, idx t *adjncy,
        idx t *vwgt, idx t *vsize, idx t *adjwgt, idx t *nparts, real t *tpwgts,
        real t ubvec, idx t *options, idx t *objval, idx t *part)           */

        int ret = METIS_PartGraphKway(&nVertices, &nWeights, xadj, adjncy,
                                       NULL, NULL, NULL, &nParts, NULL,
                                       NULL, options, &objval, part);


		vector<int> pp;   
		vector<vector<int> > t_pp(numprocs, vector<int>(0));

		for (int i = 0; i < n; i++) {
			for (int j = 0; j < n; j++) {
				for(int k = 0; k < n; k++) {
					t_pp[part[i*n*n + j*n + k]].push_back(i*n*n + j*n + k);
					if (i == 0 || j == 0 || k == 0 || i == n - 1 || j == n - 1 || k == n-1)
						;
					else if (i == int(n/2) && j == int(n/2) && k == int(n/2))
						;
					else
					{
						//pc[part[i * n + j]]++;
						if (procid == part[i*n*n + j*n + k])
							pp.push_back(i*n*n + j*n + k);
					}
				}
			}
		}

		/*  
		In the vector v[i*iter+j], the part/process every i*iter+j element has been allotted
		and the parts/processes all its neighbours have been allotted are stored.
		v[i*iter+j][0] is the part/process that the i*iter+j element has been allotted and
		the subsequent elements are the parts/processes its neighbours have been allotted.
		*/

		for (int i = 0; i < nVertices; i++) {
			a[part[i]]++;
			v[i].push_back(part[i]);
			for (int j = xadj[i]; j < xadj[i + 1]; j++) {
				if (find(v[i].begin(), v[i].end(), part[adjncy[j]])==v[i].end())
					v[i].push_back(part[adjncy[j]]);
			}
		}	


		/*
		In the vector comm[x][y], I store the indices of the square cells
		which need to be sent from process y to process x.
		*/


		for (int i = 0; i < nVertices; i++) {
			for  (int t = 1; t < v[i].size(); t++) {
				comm[v[i][t]][v[i][0]].push_back(i);
			}
		}

		/*
		if(procid==0)	
			{
			ofstream com;
			com.open("comm.txt");
	
			for (int temp=0;temp<nParts;temp++)
				{
				for (int l=0;l<nParts;l++)
					{
					com<<temp<<" "<<l<<" ";
					for  (int t=0; t<comm[temp][l].size();t++)
						{
						com<<comm[temp][l][t]<<" ";
		  				}
					com<<endl;
					}
				}
		   
			com.close();	
			}
		*/	
		
		/*ifstream div("sg.txt");	
		int a, b;

		while (div >> a >> b)
			{
			if(procid==b)
				A[a]=a;
			}
	

		div.close();*/

	
		/*for (int p=0;p<numprocs;p++)
			{	
			if(procid==p)
				{
				for (int i=0;i<n;i++)
					{
					for (int j=0;j<n;j++)
						{
						cout<<A[i*n+j]<<" ";
						}
					cout<<endl;
					}
				cout<<endl<<endl;
				}
			}*/
	
		/*
		temps and tempr store the actual values at i*iter+j (i.e. A[i*iter+j])
		to be sent and received from one process to the other. To identify which values need
		to be sent and map them to temps, it uses the vector comm.
		Later, MPI_Sendrecv is used to do the actual communication
		so that it doesn't get blocked.
		Then comm is used to map the values received in tempr back to A.
		*/
		
		for(int i = 0; i < nVertices; i++){
	    	free(XYZ[i]);
	    }
		
		free(XYZ);

		double  **tempr, **temps;
	
		tempr=(double**)calloc(numprocs, sizeof(double *));	
		temps=(double**)calloc(numprocs, sizeof(double *));	

		for (int temp = 0; temp < numprocs; temp++) {
			if(temp != procid)
			{
				temps[temp] = (double*)calloc(comm[temp][procid].size(), sizeof(double));
				tempr[temp] = (double*)calloc(comm[procid][temp].size(), sizeof(double));
			}
		}

		if(procid == 0)
			t3 = clock();

		for (int t = 0; t < iter; t++) {
	
			t5 = clock();

			for (int temp = 0; temp < numprocs; temp++) {
				if(temp != procid)
				{
					for (int p = 0; p < comm[temp][procid].size(); p++) {
						temps[temp][p] = A[comm[temp][procid][p]];
						//cout<<temps[temp][p]<<" "<<comm[temp][procid][p]<<" "<<procid<<endl;
					}
				}
			}
			
			t6 = clock();

			st1 += t6 - t5;

			t7 = clock();

			for (int temp = 0; temp < numprocs; temp++) {
				if(temp != procid)
				{
					MPI_Sendrecv(&temps[temp][0], comm[temp][procid].size(), MPI_DOUBLE, temp, 0, &tempr[temp][0], comm[procid][temp].size(), MPI_DOUBLE,temp, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);	
				}
			}

			t8 = clock();

			st2 += t8 - t7;

			t9 = clock();

			for (int temp = 0; temp < numprocs; temp++) {
				if(temp != procid)
				{
					for (int p = 0; p < comm[procid][temp].size(); p++) {
						A[comm[procid][temp][p]] = tempr[temp][p];
					}
				}		
			}

			t10 = clock();

			st3 += t10 - t9;
			
			t11 = clock();
			/*
			The actual computation over the whole square grid is performed.
			*/

			for (int i = 0; i < pp.size(); i++) {
				double temp = 0;
	           	for (int z = xadj[pp[i]]; z < xadj[pp[i] + 1]; z++) {
		           	temp += A[adjncy[z]] * adjnwgt[z];
		      	}
	         	A[pp[i]] = temp;		
			}

			t12 = clock();

			st4 += t12 - t11;
		}

		if(procid == 0)
			t4 = clock();

		for (int temp = 0; temp < numprocs; temp++) {
			if(temp != procid)
			{
				free(temps[temp]);
				free(tempr[temp]);
			}
		}
		
		free(temps);
		free(tempr);

		/*
		In the last iteration, all values are sent to process 0 
		so that we have the final result in a readable form.			
		*/

		vector<vector<double> > sv(numprocs,vector<double>(0));

		double **tempr2, *temps2;

		for (int i = 0; i < nVertices; i++) {
				sv[v[i][0]].push_back(A[i]);
		}
	
		if(procid != 0)
		{		
			temps2 = (double*)calloc(sv[procid].size(), sizeof(double));

			for (int temp = 0; temp < sv[procid].size(); temp++) {
				temps2[temp] = sv[procid][temp];
			}

			MPI_Send(&temps2[0], sv[procid].size(), MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
			free(temps2);
		}
	    
		//free(temps2);
	
		tempr2 = (double**)calloc(numprocs, sizeof(double *));	
	
		if (procid == 0)
		{
			for (int temp = 1; temp < numprocs; temp++) {
				tempr2[temp] = (double*)calloc(a[temp], sizeof(double));
				MPI_Recv(&tempr2[temp][0], a[temp], MPI_DOUBLE, temp, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
			}
		
			vector<int> z(nParts, 0);

			for (int i = 0; i < nVertices; i++) {
				if (v[i][0] != 0)
				{
					A[i] = tempr2[v[i][0]][z[v[i][0]]];
					z[v[i][0]]++;
				}
			}
		}
		
		for (int temp = 0; temp < a[temp]; temp++) {
			free(tempr2[temp]);
		}
		free(tempr2);
		
		
		if(procid == 0)
		{
			t2 = clock();
			double tK = ((double)t2 - (double)t1);
			double tk = ((double)t4 - (double)t3);
			s = s+tK;
			sk = sk+tk;
		}
	}	
	

	if(procid == 0)

	{
		cout << s * 1000000/(NT * CLOCKS_PER_SEC) << endl;
		cout << sk * 1000000/(NT * CLOCKS_PER_SEC) << endl;
		cout << st1 * 1000000/(NT * CLOCKS_PER_SEC) << endl;
		cout << st2 * 1000000/(NT * CLOCKS_PER_SEC) << endl;
		cout << st3 * 1000000/(NT * CLOCKS_PER_SEC) << endl;
		cout << st4 * 1000000/(NT * CLOCKS_PER_SEC) << endl;
	}

/*	
	if(procid==0)
	{
	    for(int i = 0; i < n; i++){
	        for(int j = 0; j < n; j++){
	            for(int k = 0; k < n; k++){
	                cout<<A[i*n*n+j*n+k]<<"  ";
	            }
	            cout<<endl<<endl;
	        }
	        cout<<endl<<endl<<endl<<endl;
	    }
	}
*/	

	free(A);

	MPI_Finalize();
	 
	return 0;	
}

