#if 0
 How many primitive variables are there in the simulation
#endif

#define NVAR_BEGIN 1

#define DENS_VAR 1
#define VELX_VAR 2
#define VELY_VAR 3
#define VELZ_VAR 4
#define PRES_VAR 5

#ifdef MHD
#define BMFX_VAR 6
#define BMFY_VAR 7
#define BMFZ_VAR 8
#define BPSI_VAR 9
#define ENER_VAR 10

#define NVAR_END 10
#else
#define ENER_VAR 6

#define NVAR_END 6
#endif

#define NVAR_NUMBER (NVAR_END - NVAR_BEGIN + 1)
