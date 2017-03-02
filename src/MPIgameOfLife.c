//* CMDA 3634, Spring 2017, HW03 Reference code*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* Q2 a) code starts here to include the mpi header file */
#include <mpi.h>
/* Q2 a) code ends here */

// Function prototypes
void printBoard(int, int, int *);
void mpiPrintBoard(int, int, int, int, int *);
void updateBoard(int, int, int *, int *);
int *mpiGameSetup(int, int, int *, int *, FILE *);
void haloExchange(int, int, int, int, int *);

int getHowManyRows(int rank, int size, int N) {
  // this function computes the number of rows a proces with rank "rank" receives
  // N = number of rows in the entire board
  int chunk = N / size;
  int remainder = N - chunk * size;
  int Nlocal = chunk + (rank < remainder);
  return Nlocal;
}

int main(int argc, char **argv) {

  /* Q2 b): use MPI_Init to initialize MPI */
  MPI_Init(&argc, &argv);
  /* Q2 b) code ends here */

  int rank, size;
  /* Q2 c): use MPI_Comm_rank and MPI_Comm_size to get rank and size */
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  /* Q2 c) code ends here */

  int N; // width of local board
  int M; // height of local board
  // Game Parameters are read from the file
  // argv[0] is the name of the executable
  // argv[1] is the name of the file we are running
  // compile:  mpicc main.c -o gameOfLife
  // run:  mpirun -n NUMBER_OF_PROCESSES ./gameOfLife inputFile.tx

  // open file for reading
  FILE *fp = fopen(argv[1], "r");

  if (fp == NULL) {
    printf("Game Of Life: could not load input file %s\n", argv[1]);
    exit(0);
  }
  int *boardA = mpiGameSetup(rank, size, &N, &M, fp);
  int Nlocal = getHowManyRows(rank, size, N);
  int *boardB = (int *)calloc((Nlocal + 2) * (M + 2), sizeof(int));

  // Print Initial Board
  if (rank == 0) {
    printf("Initial condition\n\n");
  }
  mpiPrintBoard(Nlocal, M, rank, size, boardA);

  int t = 0;

  // Start Game
  int T = 60; // number of steps
  while (t < T) {

    ++t;
    /* Q4 f): call haloExchange on boardA data */
    haloExchange(N, M, rank, size, boardA);
    /* Q4 f) code ends here */
    updateBoard(Nlocal, M, boardA, boardB);


    if (rank == 0) {
      printf("After %d iterations\n", t);
    }
    

    /* Q5 c): call mpiPrintBoard to print the local board  */
    mpiPrintBoard(Nlocal, M, rank, size, boardB);
    /* Q5 c) code ends here */

    if (t == T)
      break;

    ++t;
    /* Q4 f): call haloExchange on boardB data */
    haloExchange(N, M, rank, size, boardB);
    /* Q4 f) code ends here */
    updateBoard(Nlocal, M, boardB, boardA);

    if (rank == 0) {
      printf("After %d iterations\n", t);
    }

    /* Q5 c): call mpiPrintBoard to print the local board*/
    mpiPrintBoard(Nlocal, M, rank, size, boardA);
    /* Q5 c) code ends here */
  }

  // Finish
  free(boardA);
  free(boardB);


  /* Q2 d): use MPI_Finalize to finalize MPI */
  MPI_Finalize();
  /* Q2 d) code ends here */
  return (0);
}

// HW03 functions
int *mpiGameSetup(int rank, int size, int *sizeX, int *sizeY, FILE *fp) {

  if (rank == 0) {
    // keep reading the file until you find $Size
    char buf[BUFSIZ];
    do {
      fgets(buf, BUFSIZ, fp);
    } while (!strstr(buf, "$Size"));

    // read the size
    fgets(buf, BUFSIZ, fp);
    int N, M;
    sscanf(buf, "%d %d", &N, &M);
    sizeX[0] = N;
    sizeY[0] = M;
    printf("just read N M %d %d \n", N, M);

    // Initialize board
    int *initialBoard = (int *)calloc((N + 2) * (M + 2), sizeof(int));

    // read number of updates
    int T;

    do {
      fgets(buf, BUFSIZ, fp);
    } while (!strstr(buf, "$Updates"));
    fgets(buf, BUFSIZ, fp);
    sscanf(buf, "%d", &T);


    printf("number of updates: %d\n", T);
    int numAlive;
    // next, scan for how many alive cells you have
    do {
      fgets(buf, BUFSIZ, fp);
    } while (!strstr(buf, "$Alive"));

    // read the number of alive cells
    fgets(buf, BUFSIZ, fp);
    sscanf(buf, "%d", &numAlive);
    printf("initial number of alive cells:  %d \n", numAlive);

    // allocate the alive list (one list per every dimension)
    int *LiveList_i = (int *)calloc(numAlive, sizeof(int));
    int *LiveList_j = (int *)calloc(numAlive, sizeof(int));

    for (int i = 0; i < numAlive; i++) {
      fgets(buf, BUFSIZ, fp);
      sscanf(buf, "%d %d", &LiveList_i[i], &LiveList_j[i]);
    }
    // Close file
    fclose(fp);
    // Spawn Cells
    for (int n = 0; n < numAlive; ++n) {
      int i = LiveList_i[n];
      int j = LiveList_j[n];
      initialBoard[j + i * (N + 2)] = 1;
    }
    free(LiveList_i);
    free(LiveList_j);

    /* Q3 a): create a for-loop and use MPI_Send inside the loop
       to send the size of initialBoard to an appropriate process*/
    for (int i = 1; i < size; i++) {
      MPI_Send(&M, 1, MPI_INT, i, 0, MPI_COMM_WORLD);
      MPI_Send(&N, 1, MPI_INT, i, 0, MPI_COMM_WORLD);
    }
    /* Q3 a) code ends here */

    /* Q3 c): create a for-loop and use MPI_Send inside the loop
       to send rows of initialBoard to an appropriate process*/
    for (int i = 1; i < size; i++) {
      int k;
      int sum = 0;
      for (k = 0; k < i; ++k) {
        sum += getHowManyRows(k, size, N);
      }
      int localNi = getHowManyRows(i, size, N);
      MPI_Send(initialBoard + sum*(M+2), (M+2)*localNi, MPI_INT, i, 0, MPI_COMM_WORLD);
    }
    /* Q3 c) code ends here */

    int localN = getHowManyRows(0, size, N);

    int *boardA = (int *)calloc((localN + 2) * (M + 2), sizeof(int));
    for (int i = M + 2; i < (localN + 1) * (M + 2); i++) {
      boardA[i] = initialBoard[i];
    }
    return boardA;
  } 
  
  else {

    int N, M;
    /* Q3 a): use MPI_Recv to receive the size of the board (N rows, M columns)
     */
    MPI_Status status1;
    MPI_Status status2;
    MPI_Recv(&M, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, &status1);
    MPI_Recv(&N, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, &status2);
    /* Q3 a) ends here */
    sizeX[0] = N;
    sizeY[0] = M;

    /* Q3 b): use appropriate algebraic operation to calculate how many rows the
     * process receives*/
    int localN = getHowManyRows(rank, size, N);
    /* Q3 b) ends here */

    /* Q3 d): allocate a local board boardA of size (localN+2)*(M+2) and use MPI_Recv to
      receive the data and place it in localBoard*/
    int *boardA = (int *)calloc((localN + 2) * (M + 2), sizeof(int));
    MPI_Status status3;
    MPI_Recv(boardA + M+2, (M+2)*localN, MPI_INT, 0, 0, MPI_COMM_WORLD, &status3);

    /* Q3 d) ends here */

    return boardA;
  }
}
void mpiPrintBoard(int N, int M, int rank, int size, int *board) {

#if 0 
  for (int RankIndex = 0; RankIndex < size; ++RankIndex) {
    if (rank == RankIndex) {
      /* Q5 a) call serial printBoard here */
      printBoard(N, M, board);
      /* Q5 a) code ends here */
    }

    /* Q5 b): use MPI_Barrier here */
    MPI_Barrier(MPI_COMM_WORLD);
    /* Q5 b) code ends here */
  }
  #else

int *tmp = (int*) calloc((N+2)*(M+2), sizeof(int));

if(rank==0){
  printf("-------------------------\n");
  printBoard(N, M, board);
}

int tag = 999;
for(int r=1;r<size;++r){
   if(rank==0){
    MPI_Status status;
    MPI_Recv(tmp, (N+2)*(M+2), MPI_INT, r, tag, MPI_COMM_WORLD, &status);
    printBoard(N,M,tmp);
   }
   else if(r==rank){

    MPI_Send(board, (N+2)*(M+2), MPI_INT, 0, tag, MPI_COMM_WORLD);
   }



}

  #endif
}

void haloExchange(int N, int M, int rank, int size, int *board) {

  if(size==1) return;

  int localN = getHowManyRows(rank, size, N);
  if (rank == 0) {

    /* Q4 b) use MPI_Send to send board data to process 1 and MPI_Recv to receive
     * board data from process 1*/
    MPI_Status status;
    MPI_Send(board + localN*(M+2),     M+2, MPI_INT, 1, 0, MPI_COMM_WORLD);
    MPI_Recv(board + (localN+1)*(M+2), M+2, MPI_INT, 1, 0, MPI_COMM_WORLD, &status);
    /* Q4 b) code ends here */
  } else {
    if (rank == (size - 1)) {

      /* Q4 c) use MPI_Recv to receive board data from second to last process and
       * use MPI_Send to send board data to second to last process*/ 
      MPI_Status status;
      MPI_Recv(board,         M+2, MPI_INT, size-2, 0, MPI_COMM_WORLD, &status);
      MPI_Send(board + (M+2), M+2, MPI_INT, size-2, 0, MPI_COMM_WORLD);
      /* Q4 c) code ends here */
    } else { // all other processes


      /* Q4 d) use MPI_Recv to send board data from process rank+1 and use
       * MPI_Send to send board data to process rank-1 */
      /* Q4 e)  use MPI_Send to send board data to process rank-1 and MPI_Recv to
       * receive board data from process rank+1 */
        //Sends data to process k+1
        //Receives data from process  k-1 
        //Sends data to process k-1
        //Receives data from process k+1
      MPI_Status status;
      MPI_Send(board + localN*(M+2),     M+2, MPI_INT, rank + 1, 0, MPI_COMM_WORLD);
      MPI_Recv(board,                    M+2, MPI_INT, rank - 1, 0, MPI_COMM_WORLD, &status);
      MPI_Send(board + (M+2),            M+2, MPI_INT, rank - 1, 0, MPI_COMM_WORLD);
      MPI_Recv(board + (localN+1)*(M+2), M+2, MPI_INT, rank + 1, 0, MPI_COMM_WORLD, &status);
      /* Q4 d) and Q4 e) code ends here */
    }
  }
}


// HW02 functions
void updateBoard(int N, int M, int *oldBoard, int *newBoard) {
  // Update the board
  for (int i = 1; i < N + 1; ++i) { // starting at 1 to skip boundary layer
    for (int j = 1; j < M + 1; ++j) { // starting at 1 to skip boundary layer

      // Make useful indices
      int cell = j + i * (M + 2);            // Current cell
      int cellBelow = j + (i - 1) * (M + 2); // Cell below it
      int cellAbove = j + (i + 1) * (M + 2); // Cell above it

      // Split the sum into the 3 above, 3 below, and 2 level neighbors
      int sumBelow = oldBoard[cellBelow] + oldBoard[cellBelow - 1] + oldBoard[cellBelow + 1];
      int sumLevel = oldBoard[cell - 1]                            + oldBoard[cell + 1];
      int sumAbove = oldBoard[cellAbove] + oldBoard[cellAbove - 1] + oldBoard[cellAbove + 1];

      // Compute the sum
      int sum = sumBelow + sumLevel + sumAbove;

      // Get the current state of the cell
      int oldState = oldBoard[cell];

      // Change the state of the cell
      int newState;

      // If the cell was alive
      if (oldState == 1) {

        // Exactly 2 or 3 neighbors
        if (sum == 2 || sum == 3) {
          newState = 1;
        }

        // More then 3 or less then 2 neighbors
        else {
          newState = 0;
        }
      }

      // If the cell was dead
      else {

        // Exactly 3 neighbors
        if (sum == 3) {
          newState = 1;
        }

        else {
          newState = 0;
        }
      }

      // Update new board
      newBoard[cell] = newState;
    }
  }
}

void printBoard(int N, int M, int *board) {
  // this is a reference serial version
  // Formatted to start in top left corner, moving across each row
  for (int i = 1; i < N + 1; ++i) {
    for (int j = 1; j < M + 1; ++j) { // starting at 1 to skip boundary layer

      // Cell number and state
      int cell = j + i * (M + 2);
      int state = board[cell];
      // printf("cell: %d\n",cell);
      if (state == 1) {
        printf("X ");
      }

      else {
        printf("0 ");
      }
    }
    printf("\n");
  }
}