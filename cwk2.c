
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include "cwk2_extra.h"

int main(int argc, char **argv)
{
    int i;
    int rank, numProcs;
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &numProcs);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    if ((numProcs & (numProcs - 1)) != 0 || numProcs > 256)
    {
        // Only display the error message from one processes, but finalise and quit all of them.
        if (rank == 0)
            printf("ERROR: Launch with a number of processes that is a power of 2 (i.e. 2, 4, 8, ...) and <=256.\n");

        MPI_Finalize();
        return EXIT_FAILURE;
    }

    float *globalData = NULL;
    int globalSize = 0;
    int localSize = 0;

    if (rank == 0)
    {
        globalData = readDataFromFile(&globalSize);
        if (globalData == NULL)
        {
            MPI_Finalize(); // Should really communicate to all other processes that they need to quit as well ...
            return EXIT_FAILURE;
        }
        printf("Rank 0: Read in data set with %d floats.\n", globalSize);
    }
    double startTime = MPI_Wtime();
    
    //Calculate and broadcast local size
    if (rank == 0)
        localSize = globalSize / numProcs;

    MPI_Bcast(&localSize, 1, MPI_FLOAT, 0, MPI_COMM_WORLD);

    float *localData = (float *)malloc(localSize * sizeof(float));
    if (!localData)
    {
        printf("Error creating array on rank %d.\n", rank);
        MPI_Finalize();
        return EXIT_FAILURE;
    }

    //Scatter to other processes
    MPI_Scatter(globalData, localSize, MPI_FLOAT, localData, localSize, MPI_FLOAT, 0, MPI_COMM_WORLD);

    /*
        Calculate mean values
    */
    float count = 0.0f;
    float localMean = 0.0f;
    float cumulativeMean = 0.0f;
    int localIndex;

    for (localIndex = 0; localIndex < localSize; localIndex++)
        count += localData[localIndex];

    localMean = count / localSize;

    MPI_Reduce(&localMean, &cumulativeMean, 1, MPI_FLOAT, MPI_SUM, 0, MPI_COMM_WORLD);

    /*
        Variance calculation using binary tree reduction
    */
    float localVariance = 0.0f;
    count = 0.0f;
    for (localIndex = 0; localIndex < localSize; localIndex++)
        count += pow((localData[localIndex] - localMean), 2);

    localVariance = count / localSize;

    //Taken from lecture slides to calculate the number of levels
    int lev = 1;
    while (1 << lev <= numProcs)
        lev++;

    int currentLevel;
    for (currentLevel = 1; currentLevel < lev; currentLevel++)
    {
        int nodesToPerformOperation = numProcs / (currentLevel * 2);
        if(rank < nodesToPerformOperation){ 
            float receivedVariance = 0.0f;
            int sourceRank = rank + (numProcs / (currentLevel * 2));
            MPI_Recv(&receivedVariance, 1, MPI_FLOAT, sourceRank, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            localVariance = (localVariance + receivedVariance) / 2;
        }else{
            int destRank = rank - (numProcs / (currentLevel * 2));
            MPI_Send(&localVariance, 1, MPI_FLOAT, destRank, 0, MPI_COMM_WORLD);
        }
    }

    //Reduce variance. This could be cleaner by putting this into the below if statement but I left it here for clarity
    float cumulativeVariance = 0.0f;
    if(rank == 0) cumulativeVariance = localVariance;

    /*
        Present results and clean up
    */
    if (rank == 0)
    {

        cumulativeMean = cumulativeMean / numProcs;

        printf("Total time taken: %g s\n", MPI_Wtime() - startTime);

        finalMeanAndVariance(cumulativeMean, cumulativeVariance);

        float sum = 0.0;
        for (i = 0; i < globalSize; i++)
            sum += globalData[i];
        float mean = sum / globalSize;

        float sumSqrd = 0.0;
        for (i = 0; i < globalSize; i++)
            sumSqrd += (globalData[i] - mean) * (globalData[i] - mean);
        float variance = sumSqrd / globalSize;

        printf("SERIAL CHECK: Mean=%g and Variance=%g.\n", mean, variance);
    }

    if (rank == 0)
        free(globalData);

    MPI_Finalize();

    return EXIT_SUCCESS;
}