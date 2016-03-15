#ifndef DATAINIT_HPP
#define DATAINIT_HPP

/*
 *This is a simple data struct to
 *save the calculated data. 
 *
 *energy should have the Energy Griddim
 *Psi has NE*NX as size
 *
 */


typedef struct SimData{
	double* energy;
	double* psi;
	SimData(size_t  ne,size_t nx){
		energy=(double*) malloc(sizeof(double)*ne);
		psi=(double*) malloc(sizeof(double)*nx*ne);
	};
	void DataInit(size_t ne,size_t nx){
		energy=(double*)malloc(sizeof(double)*ne);
		psi=(double*)malloc(sizeof(double)*nx*ne);
	}
} SimData;

/*
 *Initialization and destruction function.
 *I know it's not the most elegant way, but
 *it'll do the job!
 *
 */ 

SimData* Allocation(size_t ne,size_t nx){
	SimData* d=(SimData*)malloc(sizeof(SimData)); 
	d->DataInit(ne,nx);
	return d;
};

void Destruction(SimData* data){
	free(data->energy);
	free(data->psi);
    free(data);
};
#endif
