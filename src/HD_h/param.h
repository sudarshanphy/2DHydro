#if 0
 How many primitive variables are there in the simulation
#endif

#define NVAR_BEGIN 1

#define DENS_VAR 1
#define VELX_VAR 2
#define VELY_VAR 3
#define VELZ_VAR 4
#define PRES_VAR 5
#define ENER_VAR 6

#define NVAR_END 6
#define NVAR_NUMBER (NVAR_END - NVAR_BEGIN + 1)

#if 0
 Number of conserved variable in the simulation:
 fluxes for the quantities that we get from riemann solver
#endif 
#define NCONSVAR_BEGIN 1

#define DENS_CONS 1
#define MOMX_CONS 2
#define MOMY_CONS 3
#define MOMZ_CONS 4
#define ENER_CONS 5

#define NCONSVAR_END 5
#define NCONSVAR_NUMBER (NCONSVAR_END - NCONSVAR_BEGIN + 1)
