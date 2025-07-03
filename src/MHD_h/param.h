
#define MHD

#define MASTER_PROC 0

#define NDIM 2
#define IAXIS 1
#define JAXIS 2

#define LEFT 121
#define RIGHT 122

#define FLOW 221
#define REFLECT 222
#define PERIODIC 223
#define CUSTOM 224

#if 0
 How many primitive variables are there in the simulation
#endif

#define NVAR_BEGIN 1

#define DENS_VAR 1
#define VELX_VAR 2
#define VELY_VAR 3
#define VELZ_VAR 4
#define PRES_VAR 5
#define BMFX_VAR 6
#define BMFY_VAR 7
#define BMFZ_VAR 8
#define BPSI_VAR 9
#define ENER_VAR 10

#define NVAR_END 10

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
#define BMFX_CONS 6
#define BMFY_CONS 7
#define BMFZ_CONS 8
#define BPSI_CONS 9

#define NCONSVAR_END 9
#define NCONSVAR_NUMBER (NCONSVAR_END - NCONSVAR_BEGIN + 1)
