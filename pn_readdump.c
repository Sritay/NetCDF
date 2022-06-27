
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *    This program reads a LAMMPS NetCDF dump in parrallel.
 *    Different timesteps assigned to different procs.
 *    % to compile: cc pnread.c -o pnread -lpnetcdf
 *    % to run: aprun -n xx pnread dumpfilename.nc
 *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#include <stdio.h>
#include <stdlib.h>
#include <string.h> 
#include <unistd.h> 
#include <mpi.h>
#include <pnetcdf.h>

#define CHECK_ERR                                                                            \
    {                                                                                        \
        if (err != NC_NOERR)                                                                 \
        {                                                                                    \
            printf("Error at line %d in %s: %s\n", __LINE__, __FILE__, ncmpi_strerror(err)); \
            nerrs++;                                                                         \
            goto fn_exit;                                                                    \
        }                                                                                    \
    }

static void usage(char *argv0)
{
    char *help =
        "Usage: %s [-h] | [-q] [file_name]\n"
        "       [-h] Print help\n"
        "       [filename] input netCDF file name\n";
    fprintf(stderr, help, argv0);
}

int main(int argc, char **argv)
{
    extern int optind;
    char filename[256], str_att[NC_MAX_NAME];
    int i, rank, nprocs, err, nerrs = 0, verbose = 1, ncid, typVar, velVar, coordVar, idVar, dimid[3], *buf;
    float sum_vel = 0;
    MPI_Offset len, NREC, NATOMS, local_ny, local_nx;
    MPI_Offset start[3], count[3], stride[3] = {1, 1, 1};

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

    /* get command-line arguments */
    while ((i = getopt(argc, argv, "hq")) != EOF)
        switch (i)
        {
        case 'h':
        default:
            if (rank == 0)
                usage(argv[0]);
            MPI_Finalize();
            return 1;
        }
    if (argv[optind] == NULL)
        strcpy(filename, "dump_vel.nc");
    else
        snprintf(filename, 256, "%s", argv[optind]);

    /* open an existing file for reading -------------------------------------*/
    err = ncmpi_open(MPI_COMM_WORLD, filename, NC_NOWRITE, MPI_INFO_NULL, &ncid);
    CHECK_ERR

    /* get dimension IDs */
    err = ncmpi_inq_dimid(ncid, "frame", &dimid[0]);
    CHECK_ERR
    err = ncmpi_inq_dimid(ncid, "atom", &dimid[1]);
    CHECK_ERR
    err = ncmpi_inq_dimid(ncid, "spatial", &dimid[2]);
    CHECK_ERR

    /* get lengths for dimensions frame and atom */
    err = ncmpi_inq_dimlen(ncid, dimid[0], &NREC); // number of timesteps
    CHECK_ERR
    err = ncmpi_inq_dimlen(ncid, dimid[1], &NATOMS); // number of atoms
    CHECK_ERR

    /* get IDs of variables in dump file */
    err = ncmpi_inq_varid(ncid, "velocities", &velVar);
    CHECK_ERR
    err = ncmpi_inq_varid(ncid, "coordinates", &coordVar);
    CHECK_ERR
    err = ncmpi_inq_varid(ncid, "type", &typVar);
    CHECK_ERR
    err = ncmpi_inq_varid(ncid, "id", &idVar);
    CHECK_ERR

    /* Create arrays to store data read from file */
    int local_size = (int)(NATOMS * NREC / nprocs) + 1;
    float *coordin = (float *)malloc(local_size * 3 * sizeof(float)); // atom co-ordinates
    float *velin = (float *)malloc(local_size * 3 * sizeof(float));   // atom velocities
    int *idin = (int *)malloc(local_size * sizeof(int));              // atom IDs
    float *typin = (float *)malloc(local_size * sizeof(float));       // atom types

    // store post-processed data
    float *avg_vel = (float *)malloc(NREC * sizeof(float));
    float *v_recv = NULL;
    if (rank == 0)
    {
        v_recv = (float *)malloc(nprocs * 2 * sizeof(float));
    }
    float v_send[2];

    /* local start and stop */
    start[1] = 0;
    start[0] = NREC / nprocs * rank;
    start[2] = 0;
    count[0] = 1;
    count[1] = NATOMS;
    count[2] = 3;

    MPI_Barrier(MPI_COMM_WORLD);
    double tstart = MPI_Wtime();

    /* read variables in collective mode */
    for (MPI_Offset mystart = NREC / nprocs * rank; mystart < NREC / nprocs * (rank + 1); mystart++)
    {
        start[0] = mystart;
        err = ncmpi_get_vara_float_all(ncid, velVar, start, count, velin);
        CHECK_ERR
        err = ncmpi_get_vara_int_all(ncid, idVar, start, count, idin);
        CHECK_ERR

        // do some post processing with the data, e.g., calc avg x vel
        sum_vel = 0;
        for (i = 0; i < local_size; i++)
        {
            sum_vel += velin[i * 3];
        }
        sum_vel = sum_vel / NATOMS;
        v_send[0] = mystart;
        v_send[1] = sum_vel;
        MPI_Request request;
        MPI_Igather(&v_send, 2, MPI_FLOAT, v_recv, 2, MPI_FLOAT, 0, MPI_COMM_WORLD, &request);
        MPI_Wait(&request, MPI_STATUS_IGNORE);

        // store/output post processed data
        if (rank == 0)
        {
            printf("curr tstep = %d of total %d [ %0.3f ]\n", mystart, NREC / nprocs, (float)(mystart + 1) / (NREC / nprocs) * 100);
            for (i = 0; i < nprocs; i++)
            {
                avg_vel[(int)v_send[i * 2]] = v_send[i * 2 + 1];
            }
        }
    }
    CHECK_ERR

    free(velin);
    err = ncmpi_close(ncid);
    CHECK_ERR

    MPI_Barrier(MPI_COMM_WORLD);

    double tend = MPI_Wtime();
    if (rank == 0)
        printf("total time taken = %0.3f \n", tend - tstart);

    /* check if there is any PnetCDF internal malloc residue */
    MPI_Offset malloc_size, sum_size;
    err = ncmpi_inq_malloc_size(&malloc_size);
    if (err == NC_NOERR)
    {
        MPI_Reduce(&malloc_size, &sum_size, 1, MPI_OFFSET, MPI_SUM, 0, MPI_COMM_WORLD);
        if (rank == 0 && sum_size > 0)
            printf("heap memory allocated by PnetCDF internally has %lld bytes yet to be freed\n",
                   sum_size);
    }

fn_exit:
    MPI_Finalize();
    return (nerrs > 0);
}
