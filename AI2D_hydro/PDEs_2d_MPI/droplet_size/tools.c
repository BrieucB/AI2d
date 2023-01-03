#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <float.h>
#include <string.h>
#include <time.h>
#include <mpi.h>

#include "tools.h"
#include "allocation.h"
#include "exit_if.h"
#include "erreurs.h"

void printdmatrix(double **mat,
		  int lx,
		  int ly)
{
  for(int j=0 ; j<ly ; j++)
    {
      for(int i=0 ; i<lx ; i++)
	{
	  printf("%.2f ", mat[i][ly-1-j]);
	}
      printf("\n");
    }
  printf("\n");
}

void printContigdmatrix(double **mat,
			int lx,
			int ly)
{
  for(int j=0 ; j<ly ; j++)
    {
      for(int i=0 ; i<lx ; i++)
	{
	  printf("%2.2f ", mat[0][j+i*ly]);
	}
      printf("\n");
    }
  printf("\n");   
}



void init_frame(int nbx_box,
		int nby_box,
		int *frame)
{
  // frame[4*rank] = rank of the left box
  // frame[4*rank+1] = rank of the upper box
  // frame[4*rank+2] = rank of the bottom box
  // frame[4*rank+3] = rank of the right box
  int box, nei;
  
  for(box=0 ; box<nbx_box*nby_box ; box++)
    {
      /* Left box */
      nei=box-1;
      if((nei+nbx_box)%nbx_box>box%nbx_box)
	nei+=nbx_box;
      frame[4*box]=nei;

      /* Top box*/
      nei=box+nbx_box;
      if(nei>nbx_box*nby_box-1)
	nei=box%nbx_box;
      frame[4*box+1]=nei;
      
      /* Bottom box */
      nei=box-nbx_box;
      if(nei<0)
	nei=box+(nby_box-1)*nbx_box;
      frame[4*box+2]=nei;
      
      /* Right box */
      nei=box+1;
      if(nei%nbx_box==0)
	nei-=nbx_box;
      frame[4*box+3]=nei;
      
    }
  return ;
}

void initFields(int rank,
		double **rho_tot,
		double **m_tot,
		double **m_rec,
		int Nx,
		int Ny,
		int Ny2,
		int Nx_box,
		int Ny_box,
		double ds,
		double rhof,
		double rhol,
		double ml,
		double r0
		)
{

  if(rank==0)
    {

      /* Initialize to liquid phase */
      for(int x=0 ; x<Nx ; x++)
        {
          for(int y=0 ; y<Ny ; y++)
            {
              rho_tot[x][y]=rhol;
              m_tot[x][y]=ml;
            }
	}

      /* Add a perturbation */
      int Nxfront=(int) floor((Nx)/2);

      for(int x=Nxfront-((int) floor(r0/ds))+1 ; x<Nxfront+((int) floor(r0/ds)); x++)
        {
          double dsx=((double) x- (double) Nxfront);

          for(int y=Ny2-((int) floor(r0/ds))+1 ; y<Ny2+((int) floor(r0/ds)) ; y++)
            {
              double dsy=((double) y - (double) Ny2);

              /* printf("x = %d, y = %d\n", x, y); */
              if(ds*ds*(dsx*dsx+dsy*dsy)<=r0*r0)
                {
                  rho_tot[x][y] += (rhof/0.148)*exp(-1./(1.-ds*ds*(dsx*dsx+dsy*dsy)/(r0*r0)));
                  m_tot[x][y] -= (rhof/0.148)*exp(-1./(1.-ds*ds*(dsx*dsx+dsy*dsy)/(r0*r0))); // ADD A BLOB OF containing rhof*pi*r0^2 particles
                }
            }
        }


      /* Debug */
      /* int cpt=1; */
      /* for(int x=0 ; x<Nx ; x++) */
      /*   { */
      /*     for(int y=0 ; y<Ny ; y++) */
      /*       { */
      /*         rho_tot[x][y]=cpt; */
      /*         m_tot[x][y]=cpt; */
      /* 	      cpt++; */
      /*       } */
      /* 	} */

    }

  return;
}

void scatterTotalMatrices(int rank,
			  int Nx_box, 
			  int Ny_box,
			  int Ny,
			  int nbx_box,
			  int nby_box,
			  int ncpu,
			  double **m_tot,
			  double **rho_tot,
			  double **rho_loc,
			  double **m_loc,
			  double **rho_rec,
			  double **m_rec
			  )
{

  int dest;

  /* Send the submatrices of *_tot matrices to each appropriate cpu */
  if(rank==0)
    {
      
      MPI_Datatype init_mat;
      // MPI_TYPE_VECTOR (count, blocklength, stride, oldtype, newtype)
      MPI_Type_vector(Nx_box, Ny_box, Ny, MPI_DOUBLE, &init_mat);
      MPI_Type_commit(&init_mat);
      
      for(dest=1 ; dest<ncpu ; dest++)
        {
          MPI_Send(&(m_tot[(int) ((dest%nbx_box)*Nx_box)][(int) ((dest/nbx_box)*Ny_box)]), 1, init_mat, dest, 1, MPI_COMM_WORLD);
	  /* printf("Sent m_tot from %d to %d\n", rank, dest); */
        }
    }

  if(rank!=0)
    {
      MPI_Recv(&m_rec[0][0], Nx_box*Ny_box, MPI_DOUBLE, 0, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      /* printf("Received m_tot on %d from %d\n", rank, 0); */
  
    }
  
  if(rank==0)
    {      
      MPI_Datatype init_mat;

      MPI_Type_vector(Nx_box, Ny_box, Ny, MPI_DOUBLE, &init_mat);
      MPI_Type_commit(&init_mat);

      for(dest=1 ; dest<ncpu ; dest++)
	{
	  MPI_Send(&(rho_tot[(dest%nbx_box)*Nx_box][(dest/nbx_box)*Ny_box]), 1, init_mat, dest, 1, MPI_COMM_WORLD);
	  /* printf("Sent rho_tot from %d to %d\n", rank, dest); */
	}

    }

  if(rank!=0)
    {
      MPI_Recv(&rho_rec[0][0], Nx_box*Ny_box, MPI_DOUBLE, 0, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      /* printf("Received rho_tot on %d from %d\n", rank, 0); */

      for(int x=0 ; x<Nx_box ; x++)
        {
          for(int y=0 ; y<Ny_box ; y++)
            {

              rho_loc[x+1][y+1]=rho_rec[x][y];
              m_loc[x+1][y+1]=m_rec[x][y];

            }
        }

    }
  
  /* Fill *_loc matrices on node n°0 */
  if(rank==0)
    {

      for(int x=0 ; x<Nx_box ; x++)
	{
	  for(int y=0 ; y<Ny_box ; y++)
	    {

	      rho_loc[x+1][y+1]=rho_tot[x][y];
	      m_loc[x+1][y+1]=m_tot[x][y];

	    }
	}

    }
  
  return;
  
}

void sendBorders(int rank,
		 double **m_loc,
		 double *rec_hborder,
		 int *frame,
		 int Nx_box,
		 int Ny_box,
		 MPI_Datatype sent_hborder,
		 int ncpu		 
		 )
{

  int i, x;

  for(i=0 ; i<ncpu; i+=2)
    {
      // Send to left
      if(i==rank)
	{
	  MPI_Send(&(m_loc[1][0]), Ny_box+2, MPI_DOUBLE, frame[4*rank], 1, MPI_COMM_WORLD);
  	}

      // Receive from right
      else if(i==frame[4*rank+3])
	{
	  MPI_Recv(&(m_loc[Nx_box+1][0]), Ny_box+2, MPI_DOUBLE, i, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	}

      // Send to right
      if(i==rank)
	{
	  MPI_Send(&(m_loc[Nx_box][0]), Ny_box+2, MPI_DOUBLE, frame[4*rank+3], 1, MPI_COMM_WORLD);
	}

      // Receive from left      
      else if(i==frame[4*rank])
	{
	  MPI_Recv(&(m_loc[0][0]), Ny_box+2, MPI_DOUBLE, i, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	}

      // Send to top
      if(i==rank)
	{
	  MPI_Type_vector(Nx_box+2, 1, Ny_box+2, MPI_DOUBLE, &sent_hborder);
	  MPI_Type_commit(&sent_hborder);
	  MPI_Send(&(m_loc[0][Ny_box]), 1, sent_hborder, frame[4*rank+1], 1, MPI_COMM_WORLD);
	}

      // Receive from bottom
      else if(i==frame[4*rank+2])
	{
	  MPI_Recv(&(rec_hborder[0]), Nx_box+2, MPI_DOUBLE, i, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	  for(x=0 ; x<Nx_box+2 ; x++)
	    m_loc[x][0]=rec_hborder[x];
	}

      
      // Send to bottom
      if(i==rank)
	{
	  MPI_Type_vector(Nx_box+2, 1, Ny_box+2, MPI_DOUBLE, &sent_hborder);
	  MPI_Type_commit(&sent_hborder);
	  MPI_Send(&(m_loc[0][1]), 1, sent_hborder, frame[4*rank+2], 1, MPI_COMM_WORLD);
	}
      
      // Receive from top
      else if(i==frame[4*rank+1])
	{
	  MPI_Recv(&(rec_hborder[0]), Nx_box+2, MPI_DOUBLE, i, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	  for(x=0 ; x<Nx_box+2 ; x++)
	    m_loc[x][Ny_box+1]=rec_hborder[x];

	}



      
    }

  for(i=1 ; i<ncpu; i+=2)
    {
      // Send to left
      if(i==rank)
	{
	  MPI_Send(&(m_loc[1][0]), Ny_box+2, MPI_DOUBLE, frame[4*rank], 1, MPI_COMM_WORLD);
  	}

      // Receive from right
      else if(i==frame[4*rank+3])
	{
	  MPI_Recv(&(m_loc[Nx_box+1][0]), Ny_box+2, MPI_DOUBLE, i, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	}

      // Send to right
      if(i==rank)
	{
	  MPI_Send(&(m_loc[Nx_box][0]), Ny_box+2, MPI_DOUBLE, frame[4*rank+3], 1, MPI_COMM_WORLD);
	}

      // Receive from left      
      else if(i==frame[4*rank])
	{
	  MPI_Recv(&(m_loc[0][0]), Ny_box+2, MPI_DOUBLE, i, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	}

      // Send to top
      if(i==rank)
	{
	  MPI_Type_vector(Nx_box+2, 1, Ny_box+2, MPI_DOUBLE, &sent_hborder);
	  MPI_Type_commit(&sent_hborder);
	  MPI_Send(&(m_loc[0][Ny_box]), 1, sent_hborder, frame[4*rank+1], 1, MPI_COMM_WORLD);
	}

      // Receive from bottom
      else if(i==frame[4*rank+2])
	{
	  MPI_Recv(&(rec_hborder[0]), Nx_box+2, MPI_DOUBLE, i, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	  for(x=0 ; x<Nx_box+2 ; x++)
	    m_loc[x][0]=rec_hborder[x];
	}

      
      // Send to bottom
      if(i==rank)
	{
	  MPI_Type_vector(Nx_box+2, 1, Ny_box+2, MPI_DOUBLE, &sent_hborder);
	  MPI_Type_commit(&sent_hborder);
	  MPI_Send(&(m_loc[0][1]), 1, sent_hborder, frame[4*rank+2], 1, MPI_COMM_WORLD);
	}
      
      // Receive from top
      else if(i==frame[4*rank+1])
	{
	  MPI_Recv(&(rec_hborder[0]), Nx_box+2, MPI_DOUBLE, i, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	  for(x=0 ; x<Nx_box+2 ; x++)
	    m_loc[x][Ny_box+1]=rec_hborder[x];

	}
    }

  return;
}

void sendBordersOld(int rank,
		 double **m_loc,
		 double *rec_hborder,
		 int *frame,
		 int Nx_box,
		 int Ny_box,
		 MPI_Datatype sent_hborder
		 )
{

  int x;

  /********************* SEND ***********************/
  /* Send the left vertical border to the left neighbouring box */
  printf("%d : %d\n",rank, frame[4*rank]);
  printf("%f\n",m_loc[1][0]);
  MPI_Send(&(m_loc[1][0]), Ny_box+2, MPI_DOUBLE, frame[4*rank], 1, MPI_COMM_WORLD);
  printf("OK!\n");


  /* Send the top horizontal border to the top neighbouring box */
  MPI_Type_vector(Nx_box+2, 1, Ny_box+2, MPI_DOUBLE, &sent_hborder);
  MPI_Type_commit(&sent_hborder);
  MPI_Send(&(m_loc[0][Ny_box]), 1, sent_hborder, frame[4*rank+1], 1, MPI_COMM_WORLD);

  /* Send the bottom horizontal border to the bottom neighbouring box */
  MPI_Type_vector(Nx_box+2, 1, Ny_box+2, MPI_DOUBLE, &sent_hborder);
  MPI_Type_commit(&sent_hborder);
  MPI_Send(&(m_loc[0][1]), 1, sent_hborder, frame[4*rank+2], 1, MPI_COMM_WORLD);
  
  /* Send the right vertical border to the right neighbouring box */
  MPI_Send(&(m_loc[Nx_box][0]), Ny_box+2, MPI_DOUBLE, frame[4*rank+3], 1, MPI_COMM_WORLD);
  
  /******************* RECEIVE **********************/
  /* Receive from the right neighbouring box */
  MPI_Recv(&(m_loc[Nx_box+1][0]), Ny_box+2, MPI_DOUBLE, frame[4*rank+3], 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
  
  /* Receive from the bottom neighbouring box */
  MPI_Recv(&(rec_hborder[0]), Nx_box+2, MPI_DOUBLE, frame[4*rank+2], 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      
  for(x=0 ; x<Nx_box+2 ; x++)
    m_loc[x][0]=rec_hborder[x];
    
  /* Receive from the top neighbouring box */
  MPI_Recv(&(rec_hborder[0]), Nx_box+2, MPI_DOUBLE, frame[4*rank+1], 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
  for(x=0 ; x<Nx_box+2 ; x++)
    m_loc[x][Ny_box+1]=rec_hborder[x];

  /* Receive from the left neighbouring box */
  MPI_Recv(&(m_loc[0][0]), Ny_box+2, MPI_DOUBLE, frame[4*rank], 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

  return;
}

void updateFields(int Nx_box,
		  int Ny_box,
		  double **rho_loc,
		  double **m_loc,
		  double **rho_loc_old,
		  double **m_loc_old,
		  double beta,
		  double gamma,
		  double c_adv,
		  double c_diff,
		  double dt
		  )
{

  
  memcpy(rho_loc_old[0], rho_loc[0], (Ny_box+2)*(Nx_box+2)*sizeof(double));
  memcpy(m_loc_old[0], m_loc[0], (Ny_box+2)*(Nx_box+2)*sizeof(double));
  
  for(int x=1 ; x<Nx_box+1 ; x++)
    {
      for(int y=1 ; y<Ny_box+1 ; y++)
	{
	  double betaphi = beta*m_loc_old[x][y]/rho_loc_old[x][y];
	      
	  rho_loc[x][y] = rho_loc_old[x][y] 
	    - c_adv*(m_loc_old[x+1][y]-m_loc_old[x-1][y]) // Advection
	    + c_diff*(rho_loc_old[x-1][y] // Diffusion
		      + rho_loc_old[x+1][y]
		      + rho_loc_old[x][y-1]
		      + rho_loc_old[x][y+1]
		      - 4*rho_loc_old[x][y]);

	  m_loc[x][y] = m_loc_old[x][y]
	    - c_adv*(rho_loc_old[x+1][y]-rho_loc_old[x-1][y])
	    + c_diff*(m_loc_old[x-1][y]
		      + m_loc_old[x+1][y]
		      + m_loc_old[x][y-1]
		      + m_loc_old[x][y+1]
		      - 4*m_loc_old[x][y])
	    + 2.*dt*gamma*(rho_loc_old[x][y]*sinh(betaphi)
			   - m_loc_old[x][y]*cosh(betaphi));
	}
    }

  return;
}


void updateFieldsNoTail(int Nx_box,
			 int Ny_box,
			 int nbx_box,
			 double **rho_loc,
			 double **m_loc,
			 double **rho_loc_old,
			 double **m_loc_old,
			 double beta,
			 double gamma,
			 double c_adv,
			 double c_diff,
			 double dt,
			 int rank,
			 int box0,
			 int x0,
			 double ml,
			 double rhol,
			 int rank_x,
			 int shift,
			 int x0_shifted,
			 int box_x,
			 int x0_rel,
			 int x01,
			 int x01_rel,
			 int box_x1
			 )
{
  double betaphi;
  int x,y;
  
  memcpy(rho_loc_old[0], rho_loc[0], (Ny_box+2)*(Nx_box+2)*sizeof(double));
  memcpy(m_loc_old[0], m_loc[0], (Ny_box+2)*(Nx_box+2)*sizeof(double));
  
  rank_x = rank%nbx_box;

  x0_shifted = ((x0+Nx_box*nbx_box-shift-1)%(Nx_box*nbx_box));
  box_x = (int) floor(x0_shifted/Nx_box);
  x0_rel=x0_shifted%Nx_box;
  
  x01 = ((x0+Nx_box*nbx_box-shift-2)%(Nx_box*nbx_box)); // Absolute position to the left of the fluctuation front
  x01_rel = x01%Nx_box; // Relative position on the left site of the fluctuation front
  box_x1 = (int) floor(x01/Nx_box);
  
  for(x=1 ; x<Nx_box+1 ; x++)
    {

      if((x==x01_rel+1)&&(rank_x==box_x1)) // Left column
	{
	  for(y=1 ; y<Ny_box+1 ; y++)
	    {
	      betaphi = beta*m_loc_old[x][y]/rho_loc_old[x][y];
	      
	      rho_loc[x][y] = rho_loc_old[x][y] 
		- c_adv*(m_loc_old[x][y]-m_loc_old[x-1][y]) // Advection
		+ c_diff*(rho_loc_old[x-1][y] // Diffusion
			  + rho_loc_old[x][y] //rhol
			  + rho_loc_old[x][y-1]
			  + rho_loc_old[x][y+1]
			  - 4*rho_loc_old[x][y]);

	      m_loc[x][y] = m_loc_old[x][y]
		- c_adv*(rho_loc_old[x][y]-rho_loc_old[x-1][y])
		+ c_diff*(m_loc_old[x-1][y]
			  + m_loc_old[x][y] //ml
			  + m_loc_old[x][y-1]
			  + m_loc_old[x][y+1]
			  - 4*m_loc_old[x][y])
		+ 2.*dt*gamma*(rho_loc_old[x][y]*sinh(betaphi)
			       - m_loc_old[x][y]*cosh(betaphi));
	    }

	}

      else if((x==x0_rel+1)&&(rank_x==box_x)) // Right col
	{
	  for(y=1 ; y<Ny_box+1 ; y++)
	    {
	      betaphi = beta*m_loc_old[x][y]/rho_loc_old[x][y];
	      
	      rho_loc[x][y] = rho_loc_old[x][y] 
		- c_adv*(m_loc_old[x+1][y]-ml) // Advection
		+ c_diff*(rhol                 // Diffusion
			  + rho_loc_old[x+1][y]
			  + rho_loc_old[x][y-1]
			  + rho_loc_old[x][y+1]
			  - 4*rho_loc_old[x][y]);

	      m_loc[x][y] = m_loc_old[x][y]
		- c_adv*(rho_loc_old[x+1][y]-rhol)
		+ c_diff*(ml
			  + m_loc_old[x+1][y]
			  + m_loc_old[x][y-1]
			  + m_loc_old[x][y+1]
			  - 4*m_loc_old[x][y])
		+ 2.*dt*gamma*(rho_loc_old[x][y]*sinh(betaphi)
			       - m_loc_old[x][y]*cosh(betaphi));
	    }

	}

      else
	{
	  for(y=1 ; y<Ny_box+1 ; y++)
	    {
	      betaphi = beta*m_loc_old[x][y]/rho_loc_old[x][y];
	      
	      rho_loc[x][y] = rho_loc_old[x][y] 
		- c_adv*(m_loc_old[x+1][y]-m_loc_old[x-1][y]) // Advection
		+ c_diff*(rho_loc_old[x-1][y] // Diffusion
			  + rho_loc_old[x+1][y]
			  + rho_loc_old[x][y-1]
			  + rho_loc_old[x][y+1]
			  - 4*rho_loc_old[x][y]);

	      m_loc[x][y] = m_loc_old[x][y]
		- c_adv*(rho_loc_old[x+1][y]-rho_loc_old[x-1][y])
		+ c_diff*(m_loc_old[x-1][y]
			  + m_loc_old[x+1][y]
			  + m_loc_old[x][y-1]
			  + m_loc_old[x][y+1]
			  - 4*m_loc_old[x][y])
		+ 2.*dt*gamma*(rho_loc_old[x][y]*sinh(betaphi)
			       - m_loc_old[x][y]*cosh(betaphi));
	    }
	}
    }

  return;
}



void updateFieldsNoTail2(int Nx_box,
			 int Ny_box,
			 int nbx_box,
			 double **rho_loc,
			 double **m_loc,
			 double **rho_loc_old,
			 double **m_loc_old,
			 double beta,
			 double gamma,
			 double c_adv,
			 double c_diff,
			 double dt,
			 int rank,
			 int box0,
			 int x0,
			 double ml,
			 double rhol,
			 int rank_x,
			 int shift,
			 int x0_shifted,
			 int box_x,
			 int x0_rel,
			 int x01,
			 int x01_rel,
			 int box_x1,
			 int x0_mid,
			 int x0_mid_rel,
			 int box_mid
			 )
{
  double betaphi;
  int x,y;
  
  rank_x = rank%nbx_box;

  x0_shifted = ((x0+Nx_box*nbx_box-shift)%(Nx_box*nbx_box));
  box_x = (int) floor(x0_shifted/Nx_box);
  x0_rel=x0_shifted%Nx_box;
  
  x01 = ((x0+Nx_box*nbx_box-shift-2)%(Nx_box*nbx_box)); // Absolute position to the left of the fluctuation front
  x01_rel = x01%Nx_box; // Relative position on the left site of the fluctuation front
  box_x1 = (int) floor(x01/Nx_box);

  x0_mid = ((x0+Nx_box*nbx_box-shift-1)%(Nx_box*nbx_box)); // Absolute position to the left of the fluctuation front
  x0_mid_rel = x0_mid%Nx_box; // Relative position on the left site of the fluctuation front
  box_mid = (int) floor(x0_mid/Nx_box);
  
  /* if(rank==3) */
  /*   printf("x0 = %d , x0-1 = %d, x0-2 = %d\n", x0_rel, x0_mid_rel, x01_rel); */

  if(rank_x==box_mid)
    {

      for(y=1 ; y<Ny_box+1 ; y++)
	{
	  rho_loc[x0_mid_rel+1][y]=rhol;
	  m_loc[x0_mid_rel+1][y]=ml;

	  
	}
    }

  memcpy(rho_loc_old[0], rho_loc[0], (Ny_box+2)*(Nx_box+2)*sizeof(double));
  memcpy(m_loc_old[0], m_loc[0], (Ny_box+2)*(Nx_box+2)*sizeof(double));
  
  for(x=1 ; x<Nx_box+1 ; x++)
    {
      for(y=1 ; y<Ny_box+1 ; y++)
	{
	  betaphi = beta*m_loc_old[x][y]/rho_loc_old[x][y];
	      
	  rho_loc[x][y] = rho_loc_old[x][y] 
	    - c_adv*(m_loc_old[x+1][y]-m_loc_old[x-1][y]) // Advection
	    + c_diff*(rho_loc_old[x-1][y] // Diffusion
		      + rho_loc_old[x+1][y]
		      + rho_loc_old[x][y-1]
		      + rho_loc_old[x][y+1]
		      - 4*rho_loc_old[x][y]);

	  m_loc[x][y] = m_loc_old[x][y]
	    - c_adv*(rho_loc_old[x+1][y]-rho_loc_old[x-1][y])
	    + c_diff*(m_loc_old[x-1][y]
		      + m_loc_old[x+1][y]
		      + m_loc_old[x][y-1]
		      + m_loc_old[x][y+1]
		      - 4*m_loc_old[x][y])
	    + 2.*dt*gamma*(rho_loc_old[x][y]*sinh(betaphi)
			   - m_loc_old[x][y]*cosh(betaphi));
	}
    }

  
  return;
}


void gatherFields(int ncpu,
		  int rank,
		  double **rho_loc,
		  double **m_loc,
		  double **rho_tot,
		  double **m_tot,
		  double **m_rec,
		  double **rho_rec,
		  int Nx_box,
		  int Ny_box,
		  int nbx_box,
		  int nby_box,
		  int Nx,
		  int Ny
		  )
{
  

  // SEND THE LOCAL MATRICES
  if(rank!=0)//((rank==1))//||(rank==2))
    {
      MPI_Send(&(m_loc[0][0]), (Nx_box+2)*(Ny_box+2), MPI_DOUBLE, 0, 1, MPI_COMM_WORLD);
      MPI_Send(&(rho_loc[0][0]), (Nx_box+2)*(Ny_box+2), MPI_DOUBLE, 0, 2, MPI_COMM_WORLD);
      
      /* printf("Sent from core %d data :\n", rank); */
      /* printdmatrix(m_loc, Nx_box+2, Ny_box+2); */
      /* printf("\n"); */
    }


  else if(rank==0)
    {
      /* printf("Sent from core %d data :\n", rank); */
      /* printdmatrix(m_loc, Nx_box+2, Ny_box+2); */
      /* printf("\n"); */
      
      /* Write the content of cpu 0 on the total matrices */
      for(int x=1 ; x<Nx_box+1 ; x++)
	{
	  for(int y=1 ; y<Ny_box+1 ; y++)
	    {

	      int case_x = x-1;
	      int case_y = y-1;
	      m_tot[case_x][case_y] = m_loc[x][y];
	      rho_tot[case_x][case_y] = rho_loc[x][y];
	    }
	}

      
      /* printf("Before reception :\n"); */
      /* printdmatrix(m_tot, Nx, Ny); */
      /* printf("\n"); */
      
      // RECEIVE ON PROC 0
      for(int emit=1 ; emit<ncpu ; emit++)
	{
	  MPI_Recv(&(m_rec[0][0]), (Nx_box+2)*(Ny_box+2), MPI_DOUBLE, emit, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

	  MPI_Recv(&(rho_rec[0][0]), (Nx_box+2)*(Ny_box+2), MPI_DOUBLE, emit, 2, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

	  int box_x=(int) (emit%(nbx_box));
	  int box_y=(int) (emit/(nbx_box));
	  
	  /* printf("Received from core %d. Writing on site (%d, %d). m_tot : \n", emit, box_x, box_y); */
	  /* printdmatrix(m_tot, Nx, Ny); */
	  /* printf("\n"); */

	  /* printf("Writing data :\n"); */
	  /* printdmatrix(m_rec, Nx_box, Ny_box); */
	  /* printf("\n"); */

	  // Store the received matrix at the correct location
	  for(int x=1 ; x<Nx_box+1 ; x++)
	    {
	      for(int y=1 ; y<Ny_box+1 ; y++)
		{

		  int case_x = box_x*(Nx_box) + x-1;
		  int case_y = box_y*(Ny_box) + y-1;
		  /* printf("Writing on site (%d, %d)\n", case_x, case_y); */
		  
		  m_tot[case_x][case_y] = m_rec[x][y];
		  rho_tot[case_x][case_y] = rho_rec[x][y];
		}
	    }
	}
    }
  return;
}


void printMatBorders(int ncpu,
		     int rank,
		     double **rho_loc,
		     double **m_loc,
		     int Nx_box,
		     int Ny_box,
		     int nbx_box,
		     int nby_box,
		     int Nx,
		     int Ny
		     )
{
  
  double **m_print = dmatrix((long) (Nx_box+2)*nbx_box, (long) (Ny_box+2)*nby_box);
  double **m_rec = dmatrix((long) Nx_box+2, (long) Ny_box+2);
    
  
  // SEND THE LOCAL MATRICES
  if(rank!=0)
    {
      MPI_Send(&(m_loc[0][0]), (Nx_box+2)*(Ny_box+2), MPI_DOUBLE, 0, 1, MPI_COMM_WORLD);

      /* printf("Sent from core %d data :\n", rank); */
      /* printdmatrix(m_loc, Nx_box+2, Ny_box+2); */
      /* printf("\n"); */
    }

  // RECEIVE ON PROC 0
  else if(rank==0)
    {
      /* printf("Sent from core %d data :\n", rank); */
      /* printdmatrix(m_loc, Nx_box+2, Ny_box+2); */
      /* printf("\n"); */
      
      /* Write the content of cpu 0 on the total matrices */
      for(int x=0 ; x<Nx_box+2 ; x++)
	{
	  for(int y=0 ; y<Ny_box+2 ; y++)
	    {

	      int case_x = x;
	      int case_y = y;
	      m_print[case_x][case_y] = m_loc[x][y];
	    }
	}

      
      /* printf("Before reception :\n"); */
      /* printdmatrix(m_tot, Nx, Ny); */
      /* printf("\n"); */

      for(int emit=1 ; emit<ncpu ; emit++)
	{
	  MPI_Recv(&(m_rec[0][0]), (Nx_box+2)*(Ny_box+2), MPI_DOUBLE, emit, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

	  int box_x=(int) (emit%(nbx_box));
	  int box_y=(int) (emit/(nbx_box));
	  
	  /* printf("Received from core %d. Writing on site (%d, %d). m_tot : \n", emit, box_x, box_y); */
	  /* printdmatrix(m_tot, Nx, Ny); */
	  /* printf("\n"); */

	  /* printf("Writing data :\n"); */
	  /* printdmatrix(m_rec, Nx_box, Ny_box); */
	  /* printf("\n"); */

	  // Store the received matrix at the correct location
	  for(int x=0 ; x<Nx_box+2 ; x++)
	    {
	      for(int y=0 ; y<Ny_box+2 ; y++)
		{

		  int case_x = box_x*(Nx_box+2) + x;
		  int case_y = box_y*(Ny_box+2) + y;
		  /* printf("Writing on site (%d, %d)\n", case_x, case_y); */
		  
		  m_print[case_x][case_y] = m_rec[x][y];
		}
	    }

	  /* printf("m_tot :\n"); */
	  /* printdmatrix(m_tot, Nx, Ny); */
	  /* printf("\n"); */
	}
    }

  if(rank==0)
    printdmatrix(m_print, (Nx_box+2)*nbx_box, (Ny_box+2)*nby_box);
  
  free_dmatrix(m_print);
  
  return;
}


void outputFields(int rank,
		  int t,
		  double dt,
		  int Nx,
		  int Ny,
		  FILE *f_rho,
		  FILE *f_m,
		  double **rho_tot,
		  double **m_tot
		  )

{
  if(rank==0)
    {
      char f_name_m[200];
      int pos_m = 0;
      pos_m += sprintf(&f_name_m[pos_m], "f_m_t%d.dat", (int) floor(t*dt));
      //pos_m += sprintf(&f_name_m[pos_m], "f_m_t%d.dat", t);
      f_m=fopen(f_name_m, "w");

      char f_name_rho[200];
      int pos_rho = 0;
      pos_rho += sprintf(&f_name_rho[pos_rho], "f_rho_t%d.dat", (int) floor(t*dt));
      //      pos_rho += sprintf(&f_name_rho[pos_rho], "f_rho_t%d.dat", t);
      f_rho=fopen(f_name_rho, "w");

      for(int y=0 ; y<Ny ; y++) //for(int x=(xmin-100)%Nx ; x<(xmax+100)%Nx ; x++)      
	{
	  /* fprintf(f_rho, "%d ", x+x0); */
	  /* fprintf(f_m, "%d ", x+x0); */

	  for(int x=0 ; x<Nx ; x++) //for(int y=(ymin-100)%Ny ; y<(ymax+100)%Ny ; y++)  
	    {
	      fprintf(f_rho, "%f ", rho_tot[x][Ny-1-y]);
	      fprintf(f_m, "%f ", m_tot[x][Ny-1-y]);
	    }

	  fprintf(f_rho, "\n");
	  fprintf(f_m, "\n");

	}

      fclose(f_rho);
      fclose(f_m);
    }
}
