 /* ******************* DESCRIBTION OF THE PROGRAMME ********************************** */
/*
 1) This program want to approximate the percolation threshold of the area interaction process, see for instance
 Dereudre & Houdebert,  Sharp phase transition for the continuum Widom-Rowlinson model
 for a description of this Gibbs Point process.
 2) The process is sampled using a FK representation with the generalized continuum random cluster model (grcm) , see
 Houdebert,  Phase transition of the non-symmetric Continuum Potts model
 for a description of the grcm
 3) The grcm is sampled using a MCMC birth and death algorithm, see the book
 Moller & Waagepetersen, Statistical Inference and Simulation for Spatial Point Processes, section 7.1.2
 for a description of the algorithm
 
 The program does:
 
 a) sampling of the grcm using MCMC algorith, with a number of step given by the parameter "Number_MCMC" line 588
 b) remove some cluster according to the FK representation to get an area-interaction configuration
 c) test the connectivity of the center of the box to the boundary
 d) do this a number of time to get a Law of Large Number approximation of the probability, given by the parameter Number_LLN line 589
 e) increment the activity parameter and do it again, in order to have beta -> percolation_threshold(beta)
 
 This program uses the openmp parallel for loop.
 
 
 Parameters
 
 - max_p, line 60: maximum size of the sampled configuration. If too small the program crashes, if too large the program is slower
 - max_n, line 62: maximum number of neighbors of a point. same as before
 - A and B, line 583: size of the box ([0,A]x[0,B])
 - beta, line 585: inverse temperature of area
 - zz, line 586: initial activity of area (we are incrementing the activity)
- zz_max, line 588: maximal activity;
- step, line 589: increment of the activity;
- Number_MCMC, line 592: number of MCMC iteration
- Number_LLN, line 593:   number of iteration Law Large Number
 
 Output:
 
 two txt files
- "parameters.txt" containing all the parameters;
- "results.txt" containing on each line, the acticvity, the percolation probability and the experimental intensity
 
To compile the code
 g++ -fopenmp -O1 -march=native grcm_area_perco.cpp  (-O1 may be replace by -O2 or -O3)
 
 To run  ./a.out
 
 ********************************************************************** */

#include<stdio.h>
#include<math.h>
#include<stdlib.h>
#include<time.h>
#include <algorithm>
#include <iostream>
#include <fstream>


/* ******************* CONSTANTS AND STRUCTURES ************************************ */
#define pi 3.1415922654
#define max_p 18000
// this is the maximum number of points in a configuration, if too small, one get a SegFault
#define max_n 43
// Number of max neighbors of a point. if too small, SegFault
#define kissingnumber 6 //this number depends on the dimension, for d=2 this is 6
#define MY_RAND() (float)rand()/RAND_MAX

typedef struct {float coor1; float coor2;} POINT;

typedef struct {POINT centre; int cluster;} BALL; //One could add a radius here is needed

typedef struct {int nb;				//number of discs
		BALL ball[max_p];		    //discs
		int neighbors[max_p][max_n];	//neighbors[i][j] is the j-the
                                    //neighhouring disc of the i-th disc
		int noneighb[max_p];		    //noneighb[i] is the number of discs
                                    //neighbouring with the i-th disc
		int ncc;			        //number of connected components
        double helppap;             // help value for the computation of papangelou
        int number_cluster;         // help value giving the number of
                                    // balls in the cluster we are trying to remove
	       } CONFIGURATION;

typedef struct {int v[max_p]; int cluster[kissingnumber];} VECTOR;
// help structure used in procedure calculating the number of connected components when adding/removing a point in the configuration

// **********************************distance of points********************
//This function returns distance of points x and y.

float dist_of_points(POINT x, POINT y)
{ double a, b, dist;

  a = (x.coor1 - y.coor1)*(x.coor1 - y.coor1);
  b = (x.coor2 - y.coor2)*(x.coor2 - y.coor2);
  dist = sqrt(a + b);
  if(a+b<0){dist = 0;}
  return (dist);
}

/* *******************log-gamma function********************** */
//This function is used for simulation of random variable from Poisson distribution.

float ln_gam_func(float xx)
{ double x, y, tmp, ser;
  static double
    cof[6]={76.18009172947146, -86.50532032941677, 24.01409824083091,
	    -1.231739572450155, 0.1208650973866179e-2, -0.5395239384953e-5};
  int j;

  y = x = xx;
  tmp = x + 5.5;
  tmp -= (x+0.5)*log(tmp);
  ser = 1.000000000190015;
  for (j=0; j<=5; j++) ser += cof[j]/++y;
  return -tmp+log(2.5066282746310005*ser/x);
}

/* *******************Poisson distribution******************** */
//This function returns the value randomly generated from Poisson distribution.

int sim_pois_var(float lambda)
{ static float sq,alxm,g,oldm=(-1.0);
  float em, t, y;

  if (lambda < 12.0){
    if (lambda != oldm){
      oldm = lambda;
      g=exp(-lambda);
    }
    em = -1;     t = 1.0;
    do{ ++em;
      t *= MY_RAND();
    } while (t > g);
  }
  else {
    if (lambda != oldm){
      oldm = lambda;
      sq = sqrt(2.0*lambda);
      alxm = log(lambda);
      g = lambda*alxm - ln_gam_func(lambda+1.0);
    }
    do {
      do {
	y = tan(pi*MY_RAND());
	em = sq*y+lambda;
      } while (em < 0.0);
      t = 0.9 * (1.0 + y*y) * exp(em*alxm-ln_gam_func(em+1.0)-g);
    } while (MY_RAND() > t);
  }
  return (int) em ;
}

/* ****************** FILLING TESSALLATION ******************* */
// At the beginning of running the programme, this procedure fills the whole structure of configuration with zeros.

CONFIGURATION *build_conf()
{ CONFIGURATION *conf = (CONFIGURATION*)malloc(sizeof(CONFIGURATION));
  int i,j;

  conf->ncc=0;    conf->nb=0;   conf->helppap=0;

  for (i=0;i<=max_p;i++)
    {conf->ball[i].centre.coor1=0;      conf->ball[i].centre.coor2=0;
     conf->ball[i].cluster=0;
     conf->noneighb[i]=0;
     conf->number_cluster =0;
        
     for (j=0;j<=max_n;j++)
       {conf->neighbors[i][j]=0;}
    }
  return (conf);
}

/* *******************Determining the connectivity******************** */

// The use of this function is only at the beginning, when sampling a PPP configuration, because it produce a SegFault when a cluster has too many points

// Help function for determining the connectivity of the structure.
// Each disc in the same cluster are given the same value
// Input of this function is a configuration "conf", help vector "w",
// where i-th item of this vector is positive integer denoting in which component the i-th ball lies,
// parameter "p" denoting the component, which is actually considered,
// parameter "n" denoting the disc which is actually considered.
// coordonnee is the number of clusters already considered
// Output of this function is changed vector "w".

VECTOR Connect(CONFIGURATION *conf,VECTOR w,int p,int n,int coordonnee)
{ int i;
  for (i=0;i<conf->noneighb[n];i++)
    {   if (w.v[conf->neighbors[n][i]]==0)
       { w.v[conf->neighbors[n][i]]=p; //this point is in cluster p
           w.cluster[coordonnee]++;      // this cluster has one more point
        w=Connect(conf,w,p,conf->neighbors[n][i],coordonnee);
       }
    }
  return (w);
}

/* ********************** MAKING MATRIX OF NEIGHBORING BALLS ************************ */
// Input of this procedure is the old configuration and position "n" of newly added point.
// Output of this procedure is the new configuration, where the following items are changed:
// a) the help matrix "neighbors" so that in i-th line, there are the positions of balls intersecting the i-th ball.
// b) the help vector "noneighb" in which the i-th item is the number of balls intersecting the i-th ball.

CONFIGURATION *build_neighbors_add_ball(CONFIGURATION *conf,int n)
{ int i;

  conf->noneighb[n]=0;
  for (i=0;i<n;i++)
    {    if(dist_of_points(conf->ball[i].centre,conf->ball[n].centre)<=1)
        // 1 = 0.5 + 0.5 because in wr, radius are 0.5
         {    conf->noneighb[i]=conf->noneighb[i]+1;
           conf->noneighb[n]=conf->noneighb[n]+1;
           conf->neighbors[i][conf->noneighb[i]-1]=n;
           conf->neighbors[n][conf->noneighb[n]-1]=i;
       }
    }
  return(conf);
}

//  ***************** Deleting Ball  ************************** ///
//  This function remove a point in a given configuration, determine the nex connectivity (cluster and neighbor matrix) and compute helppap needed to accept or rejet this removal.
// Input: configuration "conf", location "n" of the point to remove, and the parameters alpha1 alpha 2 used only to compute helppap
// Output : the new configuration, with in particular the value helppap
CONFIGURATION *Remove_Ball(CONFIGURATION *conf,int n, double alpha1, double alpha2)
{   int i,j,k,coordonee,test,nombreboule;
    VECTOR w;
    int monT[conf->nb], monT2[conf->nb];

// Step 1: removing the ball n from the list of neighbours of the other points
//         but keeping the "old" neighbors of the removed point to use later
    for (i=0;i<conf->noneighb[n];i++)
    {   for (j=0;(conf->neighbors[conf->neighbors[n][i]][j]!=n)&&(j<conf->noneighb[conf->neighbors[n][i]]);j++){}
//             this loop is there to access the collumn where n is
            conf->noneighb[conf->neighbors[n][i]]--;    //we substract one from the number of neighbors
            conf->neighbors[conf->neighbors[n][i]][j] = conf->neighbors[conf->neighbors[n][i]][conf->noneighb[conf->neighbors[n][i]]];
//             the removed neighbor is replaced by the last one on the list
    }
// Step 2: if the removed point has no neighbour
    if (conf->noneighb[n]==0)
    {   for (i=0;i<conf->nb;i++)
        {   if (conf->ball[i].cluster==conf->ncc)
            {   conf->ball[i].cluster=conf->ball[n].cluster;    }
        }
        conf->ncc--;  conf->helppap = 1;
    }
// Step 3: otherwise we start the exploration of the created clusters
    else
    {   coordonee=0;
//         coordonnee is counting the number of "new" clusters, and allow to count the number of point in each new              cluster
        for (i=0;i<conf->nb;i++) // initialisation of the help vector w
        { w.v[i]=0;  if (i<kissingnumber){w.cluster[i]=0;}  }
        
        for (i=0;i<conf->noneighb[n];i++)
        {   // loop of all neighbours of n which have not been explored yet
            if ( w.v[conf->neighbors[n][i]] == 0 ) // not explored yet
            {   // we have found a "new" cluster
                conf->ncc++;
                w.v[conf->neighbors[n][i]] = conf->ncc;
                w.cluster[coordonee]++;
                
                for (j=0;j<conf->nb;j++)  {   monT2[j] = 0; monT[j] = 0;}   // initialisation of the help table monT(2)
                
                monT[conf->neighbors[n][i]] = 1;
                // monT stores the points that were added to the new cluster in the last iteration (= while loop)
                // Now we explore this new cluster:
                test = 1;
                while (test == 1) //we discovered new points at the last step, so we muss keep going
                {   test = 0;
                    for (j=0;j<conf->nb;j++)
                    {   if (monT[j] == 1)
                        {   for (k=0;k<conf->noneighb[j];k++)
                            {   if (w.v[conf->neighbors[j][k]] == 0)
                                {   test = 1;
                                    monT2[conf->neighbors[j][k]] = 1;
                                    w.v[conf->neighbors[j][k]] = conf->ncc;
                                    w.cluster[coordonee]++;
                                }
                            }
                        }
                    }
                    for (j=0;j<conf->nb;j++)
                    {   monT[j] = monT2[j]; monT2[j] = 0;   }
                }
                coordonee++;
            }
        }
        // At this point we know all the new cluster created by the removal of "n", with the number of points per cluster, and everything. It is stored in the help vectore w.
        
//Step 4: putting all the information stored in w into the configuration
        for (i=0;i<conf->nb;i++)
        {   if( (w.v[i]>0) && (w.v[i] < conf->ncc) )
            {conf->ball[i].cluster=w.v[i];}
        }           conf->ncc--; //we initially increased it too much so we correct it
//Step 5: computation of helppap
        conf->helppap = 1; nombreboule = 0;
        for (i=0;i<kissingnumber;i++)
        {   if (w.cluster[i]>0)
            {   nombreboule=nombreboule+w.cluster[i];
                conf->helppap = (conf->helppap)*(1 + pow(alpha2/alpha1,w.cluster[i]) );
            }
        }
        conf->helppap=(conf->helppap)/( alpha1 * (1  + pow(alpha2/alpha1,1+nombreboule)) );
    } //end of else
    
//Step 6: Removing completelly the ball "n" from the conf
    conf->nb--;
    
    if (conf->nb > n)
    {   conf->ball[n]=conf->ball[conf->nb];
//    modifying matrix of neighboring discs so that the n-th disc is replacing by the last disc:
          for (i=0;i<conf->noneighb[conf->nb];i++)
          {  conf->neighbors[n][i]=conf->neighbors[conf->nb][i];   }
          conf->noneighb[n]=conf->noneighb[conf->nb];

//        Since we moved the last point at the location where "n" was, we have to uptade the neighbor matrix
          for (i=0;i<conf->noneighb[n];i++)
          {for (j=0;(conf->neighbors[conf->neighbors[n][i]][j] != conf->nb) && (j<conf->noneighb[conf->neighbors[n][i]]);j++){}
              
          conf->neighbors[conf->neighbors[n][i]][j]=n;
          }
    }
      return (conf);
}

/* ****************** ADDING BALL ******************************* */
// This procedure adds a given point to a given configuration, update the connectivity (clusters, neighbor matrix) and compute helppap used for accepting or rejecting this proposal.
//Input: the configuration  "conf" and coordinates of the new point, and the parameters alpha1 alpha 2 used in the computation of helppa.
//Output of this procedure is the changed configuration.

CONFIGURATION *add_new_ball (CONFIGURATION *conf,float ballcoor1, float ballcoor2, double alpha1, double alpha2)
{   int i,j,helpnumberingcc,NumberConnectingCluster,helpcluster,nombreboule,k;
    VECTOR w;
    int helptableau[kissingnumber]={0};
    FILE *fw;
    
// Step 1: adding the new disc
  conf->ball[conf->nb].centre.coor1=ballcoor1;
  conf->ball[conf->nb].centre.coor2=ballcoor2;

// Step 2: making matrix of its neighbouring discs:
  conf=build_neighbors_add_ball(conf,conf->nb);
    
// Step 3: connectivity finding all the balls connected to the new one
    if (conf->noneighb[conf->nb]==0) // If the new point is not connected to aanything
    {   conf->ncc++;
        conf->ball[conf->nb].cluster = conf->ncc;
        conf->helppap = 1;
    }
    else // if the new point has neighbor(s)
    {   for (i=0;i<conf->nb;i++) //initialisation of the help vector
        { w.v[i]=0;  if (i<kissingnumber){w.cluster[i]=0;}  }
        NumberConnectingCluster = 0; //nombre de cluster que l'on connecte and counting each point in each cluster
        
        for (i=0;i<conf->noneighb[conf->nb];i++)
        {   if (w.v[conf->neighbors[conf->nb][i]] == 0) // one new cluster was connected
            {   for (j=0;j<conf->nb;j++) //counting the number of points in this cluster, and saving the location of                                                the point
                {   if (conf->ball[j].cluster == conf->ball[conf->neighbors[conf->nb][i]].cluster)
                    {   w.v[j]=1; w.cluster[NumberConnectingCluster]++; }
                }
                
//           helpnumberingcc stores the smallest number which will be used as a label for all the clusters we found
                if (NumberConnectingCluster == 0) //we found the first new cluster
                {   helpnumberingcc=conf->ball[conf->neighbors[conf->nb][i]].cluster; }
                else
                {   if (helpnumberingcc > conf->ball[conf->neighbors[conf->nb][i]].cluster)
                    {   helptableau[NumberConnectingCluster-1] = helpnumberingcc;
                        helpnumberingcc = conf->ball[conf->neighbors[conf->nb][i]].cluster;
                    }
                    else
                    {   helptableau[NumberConnectingCluster-1] = conf->ball[conf->neighbors[conf->nb][i]].cluster;
                    }
                }
                NumberConnectingCluster++;
            }
        }   // at this stage the vector w contains all points which are connected together,
            //and the w.cluster part knows how many cluster with how many points each
        // helptableau contains the label of the "former" cluster connected together
        
// Step 4: computation helppap
        conf->helppap=1; nombreboule=0;
        for (i=0;i<kissingnumber;i++)
        {   if (w.cluster[i]>0)
            {   nombreboule=nombreboule+w.cluster[i];
                conf->helppap = (conf->helppap)*( 1 + pow(alpha2/alpha1,w.cluster[i]) );
            }
        }
        conf->helppap = alpha1 * ( 1 + (pow(alpha2/alpha1,1+nombreboule)) )/(conf->helppap);

// Step 5: relabeling of clusters
        for (i=0;i<=conf->nb;i++) //labeling all the ball connected to the new one
        {   if (w.v[i]==1)
            {   conf->ball[i].cluster=helpnumberingcc;    }
        }
        conf->ball[conf->nb].cluster=helpnumberingcc;
        
        if (NumberConnectingCluster>1)
        {   std::sort(helptableau, helptableau + NumberConnectingCluster-1);
//            we need to sort help tableau, otherwise it might not label correctly

            for (i=NumberConnectingCluster -2;i>-1;i--) // relabeling cluster to have continuity of the numbering
            {   if (helptableau[i] < conf->ncc)
                {   for (j=0;j<conf->nb;j++)
                    {   if (conf->ball[j].cluster==conf->ncc)
                        {   conf->ball[j].cluster = helptableau[i];}
                    }
                }
                conf->ncc--;
            }
        }
    } // end else
  conf->nb++;
  return (conf);
}

// ******************* MCMC birth and death step **********
// This function does one birt and death, by randomly choosing if we try to add or remove a point, then doing the appropriate function to get in particular helppap used in computing papangelou intensity, and finally accepting or rejection.
// Input : almost everything
// Output : the new configuration (which might be the same)

CONFIGURATION *MCMC_step(CONFIGURATION *conf,float rho, float A, float B,double alpha1, double alpha2,float p)
{   CONFIGURATION *helpconf = (CONFIGURATION*)malloc(sizeof(CONFIGURATION));
    float u1, u2;
    double h;
    int help;
    
    *helpconf=*conf;
    u1=MY_RAND();
    
    if (u1<p)    //proposal of adding a ball
    {   helpconf=add_new_ball(helpconf,A*MY_RAND(),B*MY_RAND(),alpha1,alpha2);
        h=(helpconf->helppap)*A*B*rho/(helpconf->nb);
        
        u2=MY_RAND();
        if (u2<h)        //proposal accepted
        {   *conf=*helpconf;    }
    }
    else        //proposal of deleting a ball
    {   if (conf->nb>0)
        {   help=(int)(MY_RAND()*helpconf->nb);
        helpconf=Remove_Ball(helpconf,help,alpha1,alpha2);
        h=(helpconf->helppap)*(conf->nb)/(A*B*rho);
        u2=MY_RAND();

        if (u2<h)    //proposal accepted
        {   *conf=*helpconf;    }
        }
    }
    free(helpconf);
return(conf);
}

// ****************** Removing a cluster **********************************
// This procedure is needed to sample the area interaction process
// This procedure removes a specific cluster from the configuration, WITHOUT updating the neighbouring matrix (because not needed)
// Input  is the configuration and the cluster to remove, parameters alpha 1 and alpha2
// Input case is there to specify the probability of removal, and it is there because in our code we have alpha1>alpha2
// output is the new configuration, with conf->nombre_cluster telling how many points were removed

CONFIGURATION *remove_cluster (CONFIGURATION *conf,int c, double alpha1, double alpha2, int cas)
{   int i;
    CONFIGURATION *helpconf = (CONFIGURATION*)malloc(sizeof(CONFIGURATION));
    double help_variable;
    float u;
    
    *helpconf = *conf;
//    Step 1: removing the cluster from helpconf
    helpconf->number_cluster = 0; // this variable tells you the number of points in the removed cluster
    for (i=helpconf->nb;i>0;i--)
    {   if (helpconf->ball[i-1].cluster == c)
        {   helpconf->number_cluster++;
            helpconf->ball[i-1] = helpconf->ball[helpconf->nb - 1];
            helpconf->nb--;
        }
    }
    
    if (helpconf->number_cluster>0)
    {   for (i=0;i<helpconf->nb;i++)
        {   if (helpconf->ball[i].cluster == helpconf->ncc)
            {   helpconf->ball[i].cluster = c;    }
        }
        helpconf->ncc--;
    }
    
//    Step 2 : testing if we accept the removal
    
    if (cas == 0) //beta = alpha2 rho
    {   help_variable = 1/(1 + pow( alpha1/alpha2, helpconf->number_cluster ));}
    else // beta = alpha1 rho
    {help_variable = 1/(1 + pow( alpha2/alpha1, helpconf->number_cluster ));}
    
    u = MY_RAND();
    
    if ( u < help_variable)
    {   *conf = *helpconf;
//         printf("removed a cluster \n");
    }
    free(helpconf);
return(conf);
}

// ****************** Testing percolation **********************************

// This procedure tests if the the configuration connects the center of the box (called as abuse of notation origine) to the boundary of the same box in the area configuration
// Input : the configuration and box parameters
// Output : integer (equal to 1 if perco occurs)

int test_percolation (CONFIGURATION *conf,float A, float B, float radius =1)
{   int i,j;
    int perco = 0; // if 1, we have percolation
    int temp = 0;
    POINT origin;
    int monTableau1[conf->nb]; //les objets deja connectés à l'origine
    int monTableau2[conf->nb];  // les nouveaux objets qui sont connectés à l'origine
    int monTableau3[conf->nb]; // il sert à updater monTableau2
    
    origin.coor1 = A/2;
    origin.coor2 = B/2;
    
// Step 1: initialisation of the table;
    
    for (i=0;i<conf->nb;i++)
    {   monTableau1[i]=0; monTableau2[i]=0; monTableau3[i]=0;}
    
// Step 2: finding the points of the configuration connected to the origin
    
    for (i=0;i<conf->nb;i++)
    {   if (   dist_of_points(conf->ball[i].centre,origin) <= radius )
        {   monTableau1[i]=1; monTableau2[i]=1; temp=1;
            if ( (conf->ball[i].centre.coor1 <= 1)
                || (A-conf->ball[i].centre.coor1 <= 1)
                || (conf->ball[i].centre.coor2 <= 1)
                || (B-conf->ball[i].centre.coor2 <= 1) )
            {perco = 1; break;}
        }
    }
    
// Step 3: exploring the cluster of the origin
    while (perco == 0 && temp == 1)
    {   temp = 0;
        for (i=0;i<conf->nb;i++)
        {   if (monTableau1[i] == 0 && perco == 0)
            {   for (j=0;j<conf->nb;j++)
                {   if (monTableau2[j] == 1 && perco == 0
                        && ( dist_of_points(conf->ball[i].centre , conf->ball[j].centre)
                        <= 2*radius ) )
                    {   monTableau1[i] = 1; temp = 1; monTableau3[i] = 1;
                        if ( (conf->ball[i].centre.coor1 <= 1)
                        || ((A - conf->ball[i].centre.coor1) <= 1)
                        || (conf->ball[i].centre.coor2 <= 1)
                        || ((B - conf->ball[i].centre.coor2) <= 1))
                        {   perco = 1;}
                    }
                }
            }
        }
        for (i=0;i<conf->nb;i++)
        {monTableau2[i] = monTableau3[i]; monTableau3[i] = 0;}
    }
return(perco);
}

//  *********   Simulation of Poisson Point Process  ******  //
// This procedure sample the initial poisson configuration, and compute the connectivity (clusters, neighbor matrix)

CONFIGURATION *sim_ppp(CONFIGURATION *conf, float rho, float A, float B)
{   int i,k;

    k =  sim_pois_var(rho*A*B);
    for (i=0; i<k; i++)
    {   conf = add_new_ball(conf,A*MY_RAND(),B*MY_RAND(),1,0);    }
return(conf);
}

// *********************************************
// *********************** MAIN ****************
// *********************************************
int main(void)
{   int i,Number_MCMC,Number_LLN, cas,perco_sum,k,intensity_sum;
    float rho,A, B,p,proba_perco,intensity_exp,step,zz_max;
    double alpha1, alpha2, beta, zz;
    char filename[30];
    
    A = 100;            //box
    B = A;
    beta = 1.723;
    zz = 1.725;
    
    zz_max =1.72501;
    step = 0.1;
    
    p = 0.5;        //probability of adding a disc
    Number_MCMC = 1500000;     //number of iterations MCMC
    Number_LLN = 100;     //number of iteration Law Large Number
    
// writting the parameters in a files
    std::ofstream outfile ("parameters.txt");
    
    outfile << "beta" " " <<  std::fixed << beta << std::endl;
    outfile << "z_min" " " <<  std::fixed << zz << std::endl;
    outfile << "z_max" " " <<  std::fixed << zz_max << std::endl;
    outfile << "Box" " " <<  std::fixed << A << "x" " " << B << std::endl;
    outfile << "Number MCMC" " " <<  std::fixed << Number_MCMC << std::endl;
    outfile << "Number LLN" " " <<  std::fixed << Number_LLN << std::endl;
    
    outfile.close();

   srand(time(NULL));
    
    outfile.open("results.txt"); // opening the results file

while (zz < zz_max) // loop to increase the intensity
{   perco_sum = 0;
    intensity_sum = 0;
    
    rho = beta + zz;
    if (beta <= zz) {cas = 0; alpha1= zz/rho; alpha2 = beta/rho;}
    else {cas = 2; alpha2= zz/rho; alpha1 = beta/rho;}

//--------------------------- starting structure - Poisson process   ------------------//
 #pragma omp parallel for
for (k=0;k<Number_LLN;k++)
{
    printf("! beta: %f,z: %f,lln: %d started ! \n",beta, zz,k);
    int number_loop_mcmc,number_loop_removal_cluster,perco;
    CONFIGURATION *conf = NULL;

conf=build_conf(); 				//filling the stucture by zeros

conf=sim_ppp(conf,1.1,A,B);    //filling the stucture by discs

//--------------------------- begin of MCMC simulation ----------------------------//
    for (number_loop_mcmc=0;number_loop_mcmc < Number_MCMC;number_loop_mcmc++)
    {   conf = MCMC_step(conf,rho,A,B,alpha1,alpha2,p);
    }

// --------------- begin of removal of cluster to get area process -----------
    for (number_loop_removal_cluster = conf->ncc; number_loop_removal_cluster > 0;number_loop_removal_cluster--)
    {   conf = remove_cluster (conf,number_loop_removal_cluster,alpha1,alpha2,cas);        }
    
//  ----------------- testing percolation ------------
    
    perco = test_percolation (conf, A, B);
    
    perco_sum += perco;
    intensity_sum += conf->nb;
    free(conf);
    printf("beta: %f, z: %f,  lln: %d finished \n",beta,zz,k);
} // end of the LLN loop
    
    proba_perco =(float)perco_sum;
    intensity_exp =(float)intensity_sum;
    proba_perco = proba_perco/Number_LLN;
    intensity_exp = intensity_exp/(Number_LLN*A*B);
    
    printf("z: %f, beta: %f, proba perco: %f, intensity exp: %f \n",zz,beta,proba_perco,intensity_exp);
    
    // Writting the results into a file
    outfile << " " << zz << " " << proba_perco << " " <<  intensity_exp << std::endl;

    zz += step;
}   // end while zz loop
    
    outfile.close(); //closing the results file

} // end of main
