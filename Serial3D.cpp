#include<iostream>
#include<stdlib.h>
#include<math.h>
#include<cstdlib>
#include<ctime>
#include<time.h>
#include<vector>


using namespace std;

int main()
    {

    int n = 64, iter = 1000, NT = 10, o = 0, l = 0;
    double s=0,sk=0;
    srand (static_cast <unsigned> (time(0)));
    clock_t t1,t2,t3,t4;

    double  *A;
    A = (double*)calloc(n*n*n, sizeof(double));

    A[n*n*int(n/2) + n*int(n/2) + int(n/2)] = 1000;
    
    for(int T = 0; T < NT; T++) {
        t1=clock();
        double **XYZ;
        
        XYZ = (double**)calloc(n*n*n,sizeof(double *));  
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

       

        int nVertices = n*n*n;
        int nEdges    = ( 3*8 + 4*12*(n-2) + 5*6*(n-2)*(n-2) + 6*(n-2)*(n-2)*(n-2) )/ 2;//modify
        
        int xadj[nVertices+1];
        xadj[0]=0;

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

        int adjncy[2*nEdges];
        double adjnd[2*nEdges];
        double adjnwgt[2*nEdges];


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

        for(int i = 0; i < nVertices; i++){
            double twgt=0;
            for(int j=xadj[i];j<xadj[i+1];j++){
                twgt+=adjnd[j];
            }
            for(int j=xadj[i];j<xadj[i+1];j++){
                adjnwgt[j]=adjnd[j]/twgt;
            }
        }

        vector<int> pp;   

        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                for(int k = 0; k < n; k++) {
                    if (i == 0 || j == 0 || k == 0 || i == n - 1 || j == n - 1 || k == n-1)
                        ;
                    else if (i == int(n/2) && j == int(n/2) && k == int(n/2))
                        ;
                    else
                    {
                        pp.push_back(i*n*n + j*n + k);
                    }
                }
            }
        }

        for(int i=0;i<nVertices;i++){
            free(XYZ[i]);
        }
        free(XYZ);

        t3=clock();
        
        for(int t=0;t<iter;t++){
            for (int i = 0; i < pp.size(); i++) {
                double temp = 0;
                for (int z = xadj[pp[i]]; z < xadj[pp[i] + 1]; z++) {
                    temp += A[adjncy[z]] * adjnwgt[z];
                }
                A[pp[i]] = temp;        
            }
        }
        
        t4=clock();

                
        t2=clock();

        double tk=((double)t4-(double)t3);
        
        double t=((double)t2-(double)t1);
        
        sk=sk+tk;
        s=s+t;
    }

    cout<<s*1000000/(NT*CLOCKS_PER_SEC)<<endl;
    cout<<sk*1000000/(NT*CLOCKS_PER_SEC);
    
/*    
    for(int i = 0; i < n; i++){
        for(int j = 0; j < n; j++){
            for(int k = 0; k < n; k++){
                cout<<A[i*n*n+j*n+k]<<"  ";
            }
            cout<<endl<<endl;
        }
        cout<<endl<<endl<<endl<<endl;
    }
*/

    free(A);
}
