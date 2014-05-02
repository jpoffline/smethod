
// masterstruct.h

#ifndef STRUCTMASTER_H
#define STRUCTMASTER_H

struct MASTER{

	DATA params;
	GRIDINFO grid;
	FIELDCONTAINER field;
	POISS poiss;
	COSM cosmology;
	
	void setupMaster(struct MASTER *master, struct DATA *params, struct GRIDINFO *grid, 
				struct FIELDCONTAINER *field, struct POISS *poiss, struct COSM *cosmology){
		
		master->params = *params;
		master->grid = *grid;
		master->field = *field;
		master->poiss = *poiss;
		master->cosmology = *cosmology;
		
	} // END setupMaster()
	
	void dostuff(struct MASTER *master){
	
		cout << master->cosmology.eta << " " ;		
		cout << master->cosmology.H << " " ;
		cout << master->cosmology.a <<endl;		
	}
	
};

#endif