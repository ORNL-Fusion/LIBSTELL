      MODULE mpi_params
      INTEGER, PARAMETER :: master=0
      INTEGER, PARAMETER :: WORKER_SPLIT_KEY = 3
      INTEGER :: myid=master, numprocs, ierr_mpi
      INTEGER :: MPI_COMM_WORKERS=-1, worker_id=-1                 !communicator for worker processors only
      INTEGER :: MPI_COMM_WORKERS_OK=-1, worker_id_ok=-1           !communicator subgroup, vmec ran ok
      END MODULE mpi_params
