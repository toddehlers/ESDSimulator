Todd Ehlers 19 May, 2015
This version of Pecube is compiled with MPI so that it can run on multicore machines.  This is not
normally needed if you are doing a single simulation.  However, if you are using the
Monte Carlo or Genetic search algorithm then it can run multiple jobs at once.

The implications of this are that:
1. You have to have a version of MPI installed on your machine.  We are using OpenMPI, but others will likely
work.

2. Compile the program using scons.  On our system this is done with scons --use-mpi, and then scons -c to clean
out the .o files.

3. To run the program you need to have the pecube.in file in the same directory as the executable, and then
start the job with mpi run.  For example, on our system we do:

mpirun -n pecube

Where -n is the number of cores you want to use.  If you are only running 1 job, then there is likely no speed
difference if N=1 or N=8.  So, for 1 job, you write
 mpirun -n 1 pecube.

** End of file
