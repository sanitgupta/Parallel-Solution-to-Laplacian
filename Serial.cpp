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
    int n=256,k=1000,NT=10,o=0,l=0;
    double s=0,sk=0;
    srand (static_cast <unsigned> (time(0)));
    clock_t t1,t2,t3,t4;

    double  *A;
    A = (double*)calloc(n*n, sizeof(double));

    A[n*n/2+n/2]=100;
    
    for(int T=0;T<NT;T++) {
        t1=clock();
        double  **XY;
        
        XY=(double**)calloc(n*n,sizeof(double *));  
        //temps=(double**)calloc(numprocs,sizeof(double *));  
        
        for(int i=0;i<n*n;i++)
            XY[i]=(double*)calloc(2,sizeof(double));

        for(int i = 0; i < n; i++){
            for(int j = 0; j < n; j++){
                XY[i*n+j][0]=j;
                XY[i*n+j][1]=i;
            }
        }

       

        int nVertices = n*n;
        int nEdges    = 2*n*(n-1);
        
        int xadj[nVertices+1];
        xadj[0]=0;

        for(int i=0;i<n;i++){
            for(int j=0;j<n;j++){
                if((i==0||i==n-1)&&(j==0||j==n-1))
                    xadj[i*n+j+1]=xadj[i*n+j]+2;

                else if(i==0||i==n-1||j==0||j==n-1)     
                    xadj[i*n+j+1]=xadj[i*n+j]+3;
        
                else
                    xadj[i*n+j+1]=xadj[i*n+j]+4;
            }
        }

        int adjncy[2*nEdges];
        double adjnd[2*nEdges];
        double adjnwgt[2*nEdges];

        l=0;
        for(int i=0;i<n;i++) {
            for(int j=0;j<n;j++) {
                if(i==0&&j==0)
                {
                    adjncy[l]=i*n+j+1;
                    adjnd[l] = pow(pow(XY[adjncy[l]][0]-XY[i*n+j][0],2) + pow(XY[adjncy[l]][1]-XY[i*n+j][1],2),0.5);
                    l++;

                    adjncy[l]=(i+1)*n+j;
                    adjnd[l] = pow(pow(XY[adjncy[l]][0]-XY[i*n+j][0],2) + pow(XY[adjncy[l]][1]-XY[i*n+j][1],2),0.5);
                    l++;
                }

                else if(i==0&&j==n-1)
                {
                    adjncy[l]=i*n+j-1;
                    adjnd[l] = pow(pow(XY[adjncy[l]][0]-XY[i*n+j][0],2) + pow(XY[adjncy[l]][1]-XY[i*n+j][1],2),0.5);
                    l++;    

                    adjncy[l]=(i+1)*n+j;
                    adjnd[l] = pow(pow(XY[adjncy[l]][0]-XY[i*n+j][0],2) + pow(XY[adjncy[l]][1]-XY[i*n+j][1],2),0.5);
                    l++;
                }

                else if(i==n-1&&j==0)
                {
                    adjncy[l]=i*n+j+1;
                    adjnd[l] = pow(pow(XY[adjncy[l]][0]-XY[i*n+j][0],2) + pow(XY[adjncy[l]][1]-XY[i*n+j][1],2),0.5);
                    l++;

                    adjncy[l]=(i-1)*n+j;
                    adjnd[l] = pow(pow(XY[adjncy[l]][0]-XY[i*n+j][0],2) + pow(XY[adjncy[l]][1]-XY[i*n+j][1],2),0.5);
                    l++;
                }

                else if(i==n-1&&j==n-1)
                {
                    adjncy[l]=i*n+j-1;
                    adjnd[l] = pow(pow(XY[adjncy[l]][0]-XY[i*n+j][0],2) + pow(XY[adjncy[l]][1]-XY[i*n+j][1],2),0.5);
                    l++;

                    adjncy[l]=(i-1)*n+j;
                    adjnd[l] = pow(pow(XY[adjncy[l]][0]-XY[i*n+j][0],2) + pow(XY[adjncy[l]][1]-XY[i*n+j][1],2),0.5);
                    l++;
                }

                else if(i==0)
                {
                    adjncy[l]=i*n+j+1;
                    adjnd[l] = pow(pow(XY[adjncy[l]][0]-XY[i*n+j][0],2) + pow(XY[adjncy[l]][1]-XY[i*n+j][1],2),0.5);
                    l++;

                    adjncy[l]=i*n+j-1;
                    adjnd[l] = pow(pow(XY[adjncy[l]][0]-XY[i*n+j][0],2) + pow(XY[adjncy[l]][1]-XY[i*n+j][1],2),0.5);
                    l++;    

                    adjncy[l]=(i+1)*n+j;
                    adjnd[l] = pow(pow(XY[adjncy[l]][0]-XY[i*n+j][0],2) + pow(XY[adjncy[l]][1]-XY[i*n+j][1],2),0.5);
                    l++;
                }

                else if(j==0)
                {
                    adjncy[l]=i*n+j+1;
                    adjnd[l] = pow(pow(XY[adjncy[l]][0]-XY[i*n+j][0],2) + pow(XY[adjncy[l]][1]-XY[i*n+j][1],2),0.5);
                    l++;
                   
                    adjncy[l]=(i+1)*n+j;
                    adjnd[l] = pow(pow(XY[adjncy[l]][0]-XY[i*n+j][0],2) + pow(XY[adjncy[l]][1]-XY[i*n+j][1],2),0.5);
                    l++;    
                    
                    adjncy[l]=(i-1)*n+j;
                    adjnd[l] = pow(pow(XY[adjncy[l]][0]-XY[i*n+j][0],2) + pow(XY[adjncy[l]][1]-XY[i*n+j][1],2),0.5);
                    l++;
                }

                else if(i==n-1)
                {
                    adjncy[l]=i*n+j+1;
                    adjnd[l] = pow(pow(XY[adjncy[l]][0]-XY[i*n+j][0],2) + pow(XY[adjncy[l]][1]-XY[i*n+j][1],2),0.5);
                    l++;
                    
                    adjncy[l]=i*n+j-1;
                    adjnd[l] = pow(pow(XY[adjncy[l]][0]-XY[i*n+j][0],2) + pow(XY[adjncy[l]][1]-XY[i*n+j][1],2),0.5);
                    l++;    
                    
                    adjncy[l]=(i-1)*n+j;
                    adjnd[l] = pow(pow(XY[adjncy[l]][0]-XY[i*n+j][0],2) + pow(XY[adjncy[l]][1]-XY[i*n+j][1],2),0.5);
                    l++;
                }

                else if(j==n-1)
                {
                    adjncy[l]=i*n+j-1;
                    adjnd[l] = pow(pow(XY[adjncy[l]][0]-XY[i*n+j][0],2) + pow(XY[adjncy[l]][1]-XY[i*n+j][1],2),0.5);
                    l++;
                    
                    adjncy[l]=(i+1)*n+j;
                    adjnd[l] = pow(pow(XY[adjncy[l]][0]-XY[i*n+j][0],2) + pow(XY[adjncy[l]][1]-XY[i*n+j][1],2),0.5);
                    l++;    
                    
                    adjncy[l]=(i-1)*n+j;
                    adjnd[l] = pow(pow(XY[adjncy[l]][0]-XY[i*n+j][0],2) + pow(XY[adjncy[l]][1]-XY[i*n+j][1],2),0.5);
                    l++;
                }

                else 
                {
                    adjncy[l]=i*n+j+1;
                    adjnd[l] = pow(pow(XY[adjncy[l]][0]-XY[i*n+j][0],2) + pow(XY[adjncy[l]][1]-XY[i*n+j][1],2),0.5);
                    l++;
                    
                    adjncy[l]=i*n+j-1;
                    adjnd[l] = pow(pow(XY[adjncy[l]][0]-XY[i*n+j][0],2) + pow(XY[adjncy[l]][1]-XY[i*n+j][1],2),0.5);
                    l++;
                    
                    adjncy[l]=(i+1)*n+j;
                    adjnd[l] = pow(pow(XY[adjncy[l]][0]-XY[i*n+j][0],2) + pow(XY[adjncy[l]][1]-XY[i*n+j][1],2),0.5);
                    l++;    
                    
                    adjncy[l]=(i-1)*n+j;
                    adjnd[l] = pow(pow(XY[adjncy[l]][0]-XY[i*n+j][0],2) + pow(XY[adjncy[l]][1]-XY[i*n+j][1],2),0.5);
                    l++;
                }           
            }
        }

        for(int i = 0; i < n * n; i++){
            double twgt=0;
            for(int j=xadj[i];j<xadj[i+1];j++){
                twgt+=adjnd[j];
            }
            for(int j=xadj[i];j<xadj[i+1];j++){
                adjnwgt[j]=adjnd[j]/twgt;
            }
        }

        /*for(int i=0;i<2*nEdges;i++){
            cout<<adjnwgt[i]<<" ";
        }*/

        for(int i=0;i<n*n;i++){
            free(XY[i]);
        }
        free(XY);

        t3=clock();
        
        for(int t=0;t<k;t++){
            for(int i=1;i<n-1;i++){
                for(int j=1;j<n-1;j++){
                    if(i==n/2&&j==n/2)
                        ;
                    else
                    {
                        //A[i*n+j]=A[i*n+j]+cx*(A[(i+1)*n+j]+A[(i-1)*n+j]-2*A[i*n+j])+cy*(A[i*n+j+1]+A[i*n+j-1]-2*A[i*n+j]);
                        double temp=0;
                        for(int z=xadj[i*n+j];z<xadj[i*n+j+1];z++)
                            temp+=A[adjncy[z]]*adjnwgt[z];
                        A[i*n+j]=temp;
                    }
                }
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
    for(int i=0;i<n;i++){
        for(int j=0;j<n;j++){
            cout<<A[i*n+j]<<"  ";
        }
        cout<<endl<<endl;
    }
    */
    free(A);
}
