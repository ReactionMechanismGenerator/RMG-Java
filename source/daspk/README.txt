MRH transferred the following files from /home/rmg3/fortran_software/fortran/ODESolver/DASPK:
   ddaspk.f
   dbanpre.f
   daux.f
   dilupre.f
   dlinpk.f
   drbdpre.f
   drbgpre.f
   dsensd.f
   dsparsk.f
   adf_dummy.f
   mpi_dummy.f

I do not think we can distribute these files (hence why they were initially absent from
this directory).
They have also been added to .gitignore so that they will not be commited to the git repository (in case we later make this public .. it's hard to remove them from history once committed)