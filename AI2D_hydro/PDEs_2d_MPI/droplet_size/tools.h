void printdmatrix(double **,
		  int,
		  int);

void printContigdmatrix(double **,
			int,
			int);

void init_frame(int,
		int,
		int *);

void initFields(int,
		double **,
		double **,
		double **,
		int ,
		int ,
		int,
		int ,
		int,
		double,
		double,
		double,
		double,
		double
		);


void scatterTotalMatrices(int ,
			  int , 
			  int ,
			  int ,
			  int ,
			  int ,
			  int ,
			  double **,
			  double **,
			  double **,
			  double **,
			  double **,
			  double **
			  );


void sendBorders(int,
		 double **,
		 double *,
		 int *,
		 int,
		 int,
		 MPI_Datatype,
		 int
		 );

void updateFields(int,
		  int,
		  double **,
		  double **,
		  double **,
		  double **,
		  double ,
		  double ,
		  double ,
		  double ,
		  double 
		  );

void updateFieldsNoTail(int,
			int,
			int,
			double **,
			double **,
			double **,
			double **,
			double,
			double,
			double,
			double,
			double,
			int,
			int,
			int,
			double,
			double,
			int,
			int,
			int,
			int,
			int,
			int,
			int,
			int);

void updateFieldsNoTail2(int,
			 int,
			 int,
			 double **,
			 double **,
			 double **,
			 double **,
			 double,
			 double,
			 double,
			 double,
			 double,
			 int,
			 int,
			 int,
			 double,
			 double,
			 int,
			 int,
			 int,
			 int,
			 int,
			 int,
			 int,
			 int,
			 int,
			 int,
			 int
			 );

void printMatBorders(int,
		     int,
		     double **,
		     double **,
		     int ,
		     int ,
		     int ,
		     int ,
		     int ,
		     int 
		     );

void gatherFields(int,
		  int,
		  double **,
		  double **,
		  double **,
		  double **,
		  double **,
		  double **,
		  int ,
		  int,
		  int,
		  int,
		  int,
		  int
		  );

void outputFields(int,
		  int,
		  double,
		  int,
		  int,
		  FILE *,
		  FILE *,
		  double **,
		  double **
		  );
