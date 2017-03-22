#include <stdio.h>
#include <stdlib.h>
#include <string.h>
/* Q2 a) code starts here */
#include <mpi.h>

/* Q2 a) code ends here */

// Function prototypes
void printBoard(int, int, int *);
void mpiPrintBoard(int, int, int,int, int *);
void updateBoard(int, int, int *, int *);
int *mpiGameSetup(int, int, int *, int *, int *, FILE *);
void haloExchange(int, int, int, int, int *);
  
int globalSumStateChanges = 0;

int getHowManyRows(int rank, int size, int N) {
  // this function computes the number of rows a proces with rank "rank"
  // receives
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
  int T; //number of updates
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
  
  int *boardA = mpiGameSetup(rank, size, &N, &M, &T, fp);
  int Nlocal = getHowManyRows(rank, size, N);
  int *boardB = (int *)calloc((Nlocal + 2) * (M + 2), sizeof(int));
  // Print Initial Board
  if (rank == 0) {
    printf("Initial condition\n\n");
  }
  mpiPrintBoard(N, M, rank, size, boardA);
   
  /* Initialize t1 and t2*/
  double t1 = 0;
  double t2 = 0; 

  /* t1 is the time at this point*/
  t1 = MPI_Wtime();

  int t = 0;
  // Start Game
  //int T = 10;

  while (t < T) {
    
    ++t;
    /* Q4 f): call haloExchange on boardA data */
    if (size!=1){
      haloExchange(N, M, rank, size, boardA);
    }
    /* Q4 f) code ends here */
    updateBoard(Nlocal, M, boardA, boardB);
    
    if (globalSumStateChanges == 0) {
      if (rank == 0) {
        printf("The board has reached a steady state.\r\n");
      }
      break;
    }

    /* Q5 c): call mpiPrintBoard to print the local board  */
    if (rank == 0) {
      printf("After %d iterations\n", t);
    }
    mpiPrintBoard(N, M, rank, size, boardB);
    /* Q5 c) code ends here */
    
    if (t == T)
      break;
    
    ++t;
    /* Q4 f): call haloExchange on boardB data */
    if (size!=1){
      haloExchange(N, M, rank, size, boardB);
    }
    /* Q4 f) code ends here */
    updateBoard(Nlocal, M, boardB, boardA);
    /* Q5 c): call mpiPrintBoard to print the local board*/
    if (rank == 0) {
      printf("After %d iterations\n", t);
    }
    
    mpiPrintBoard(N, M, rank, size, boardA);
    /* Q5 c) code ends here */
  }
  
  /* t2 is the time after all board updates*/
  t2 = MPI_Wtime();

  /* Print the elapsed time*/
  //printf("Elapsed time is %f\n", t2 - t1 ); 

  // Finish
  free(boardA);
  free(boardB);
  
  /* Q2 d): use MPI_Finalize to finalize MPI */
  MPI_Finalize();
  /* Q2 d) code ends here */

  return (0);

}

// HW03 functions
int *mpiGameSetup(int rank, int size, int *sizeX, int *sizeY, int *T,  FILE *fp) {
  if (rank == 0) {
    // keep reading the file until you find $Size
    char buf[BUFSIZ];
    do {
      fgets(buf, BUFSIZ, fp);
    } while (!strstr(buf, "$Size"));
    
    // read the size
    fgets(buf, BUFSIZ, fp);
    int N, M, Ttemp;
    sscanf(buf, "%d %d", &N, &M);
    sizeX[0] = N;
    sizeY[0] = M;
    printf("just read N M %d %d \n", N, M);
    // Initialize board
    /* The line below is a part of Q3 c*/
    int *initialBoard = (int *)calloc((N + 2) * (M + 2), sizeof(int));
    
    // read number of updates
    
    
    do {
      fgets(buf, BUFSIZ, fp);
    } while (!strstr(buf, "$Updates"));
    fgets(buf, BUFSIZ, fp);
    sscanf(buf, "%d", &Ttemp);
    T[0]=Ttemp;
    
    printf("number of updates: %d\n", Ttemp);
    int numAlive;
    // next, scan for how many alive cells you have
    do {
      fgets(buf, BUFSIZ, fp);
    } while (!strstr(buf, "$Alive"));
    
    // read the number of alive cells
    fgets(buf, BUFSIZ, fp);
    sscanf(buf, "%d", &numAlive);
    printf("initial number of alive cells:  %d \n", numAlive);
    
    // allocate the alive list (one list per every dimension
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
      initialBoard[j + i * (M + 2)] = 1;
    }
    free(LiveList_i);
    free(LiveList_j);
    
    /* Q3 a): create a for-loop and use MPI_Send inside the loop
     to send the size of initialBoard to an appropriate process*/
    for (int i = 1; i < size; i++) {
      // sending to process i;
      
      int localN = getHowManyRows(i, size, N);
      int msgLength = 1;
      
      MPI_Send(&N, msgLength, MPI_INT, i, 999, MPI_COMM_WORLD);
      MPI_Send(&M, msgLength, MPI_INT, i, 998, MPI_COMM_WORLD);
      MPI_Send(&Ttemp, msgLength, MPI_INT, i, 997, MPI_COMM_WORLD);
      int mymin;
      int chunk = N / size;
      int remainder = N - chunk * size;
      
      int extra  = remainder;
      if (i<remainder){
        extra = i;
      }
      
      int start = i*chunk + extra + 1;
      
      start = start*(M+2);
      
      msgLength = localN * (M + 2);
      MPI_Send(&initialBoard[start], msgLength, MPI_INT, i, 998,
               MPI_COMM_WORLD);
    }
    /* Q3 a) code ends here */
    
    /* Q3 c): create a for-loop and use MPI_Send inside the loop
     to send rows of initialBoard to an appropriate process*/
    
    /* Q3 c) code ends here */
    
    int localN = getHowManyRows(0, size, N);
    
    int *boardA = (int *)calloc((localN + 2) * (M + 2), sizeof(int));
    for (int i = M + 2; i < (localN + 1) * (M + 2); i++) {
      boardA[i] = initialBoard[i];
    }
    return boardA;
  } else {
    MPI_Status msgStatus1, msgStatus2, msgStatus3;
    
    int N, M, Ttemp;
    /*Q3 a): use MPI_Recv to receive the size of the board (N rows, M columns)
     */
    MPI_Recv(&N, 1, MPI_INT, 0, 999, MPI_COMM_WORLD, &msgStatus1);
    MPI_Recv(&M, 1, MPI_INT, 0, 998, MPI_COMM_WORLD, &msgStatus2);
    MPI_Recv(&Ttemp, 1, MPI_INT, 0, 997, MPI_COMM_WORLD, &msgStatus3);
    /*Q3 a) ends here */
    sizeX[0] = N;
    sizeY[0] = M;
    T[0] = Ttemp;
    /*Q3 b): use appropraite algebraic operation to calculate how many rows the
     * process receives*/
    int localN = getHowManyRows(rank, size, N);
    /*Q3 b) ends here */
    /*Q3 d): allocate a localBoard of size (localN+2)*(M+2) and use MPI_Recv to
     receive the data
     and place it in localBoard*/
    int *boardA = (int *)calloc((localN + 2) * (M + 2), sizeof(int));
    
    MPI_Recv(&boardA[M + 2], localN * (M + 2), MPI_INT, 0, 998, MPI_COMM_WORLD,
             &msgStatus3);
    /*Q3 d) ends here */
    
    return boardA;
  }
}
void mpiPrintBoard(int N, int M, int rank, int size, int *board) {
  
  //  for (int RankIndex = 0; RankIndex < size; ++RankIndex) {
  //  if (rank == RankIndex) {
  /* Q5 a) call serial printBoard here */
  if (rank ==0){
    int localN;
    localN = getHowManyRows(rank, size, N);
    printBoard(localN, M, board);
    MPI_Status msgStatus1;
    
    for (int i=1; i< size; i++){
      localN = getHowManyRows(i, size, N);
      int * boardAux = (int*) calloc((localN+2)*(M+2), sizeof(int));
      MPI_Recv(&boardAux[M+2], localN*(M+2), MPI_INT, i, 999, MPI_COMM_WORLD, &msgStatus1);
      printBoard(localN, M, boardAux);
      free(boardAux);
      }
    
  }
  else{
    int localN = getHowManyRows(rank, size, N);
    MPI_Send(&board[M+2],  localN*(M+2), MPI_INT, 0, 999,MPI_COMM_WORLD);
  }
}
void haloExchange(int N, int M, int rank, int size, int *board) {
  
  int localN = getHowManyRows(rank, size, N);
  if (rank == 0) {
    
    /*Q4 b) use MPI_Send to send board data to process 1 and MPI_Recv to receive
     * board data from process 1*/
    MPI_Status msgStatus1;
    MPI_Send(&board[(localN) * (M + 2)], M + 2, MPI_INT, 1, 99, MPI_COMM_WORLD);
    MPI_Recv(&board[(localN + 1) * (M + 2)], M + 2, MPI_INT, 1, 99,
             MPI_COMM_WORLD, &msgStatus1);
    /*Q4 b) code ends here */
  } else {
    if (rank == (size - 1)) {
      
      /*Q4 c) use MPI_Recv to receive board data from second to last process and
       * use MPI_Send to send board data to second to last process*/
      
      MPI_Status msgStatus1;
      MPI_Recv(&board[0], M + 2, MPI_INT, size - 2, 99, MPI_COMM_WORLD,
               &msgStatus1);
      MPI_Send(&board[M + 2], M + 2, MPI_INT, size - 2, 99, MPI_COMM_WORLD);
      
      /*Q4 c) code ends here */
    } else { // all other processes
      
      MPI_Status msgStatus1, msgStatus2;
      /*Q4 d) use MPI_Recv to send board data from process rank+1 and use
       * MPI_Send to send board data to process rank-1 */
      MPI_Send(&board[(localN) * (M + 2)], M + 2, MPI_INT, rank + 1, 99,
               MPI_COMM_WORLD);
      MPI_Recv(&board[0], M + 2, MPI_INT, rank - 1, 99, MPI_COMM_WORLD,
               &msgStatus1);
      /*Q4 d) code ends here */
      /*Q4 e)  use MPI_Send to send board data to process rank-1 and MPI_Recv to
       * receive board data from process rank+1 */
      MPI_Send(&board[M + 2], M + 2, MPI_INT, rank - 1, 99, MPI_COMM_WORLD);
      MPI_Recv(&board[(localN + 1) * (M + 2)], M + 2, MPI_INT, rank + 1, 99,
               MPI_COMM_WORLD, &msgStatus2);
      /*Q4 e) code ends here */
    }
  }
}

// HW02 functions
void updateBoard(int N, int M, int *oldBoard, int *newBoard) {
  // Update the board

  // 2a)
  int sumStateChanges = 0;

  for (int i = 1; i < N + 1; ++i) { // starting at 1 to skip boundary layer
    for (int j = 1; j < M + 1; ++j) { // starting at 1 to skip boundary layer
      
      // Make useful indices
      int cell = j + i * (M + 2);            // Current cell
      int cellBelow = j + (i - 1) * (M + 2); // Cell below it
      int cellAbove = j + (i + 1) * (M + 2); // Cell above it
      
      // Split the sum into the 3 above, 3 below, and 2 level neighbors
      int sumBelow = oldBoard[cellBelow] + oldBoard[cellBelow - 1] +
      oldBoard[cellBelow + 1];
      int sumLevel = oldBoard[cell - 1] + oldBoard[cell + 1];
      int sumAbove = oldBoard[cellAbove] + oldBoard[cellAbove - 1] +
      oldBoard[cellAbove + 1];
      
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

      if (oldState != newState) {
        sumStateChanges = sumStateChanges + 1;
      }

    }
  }

  MPI_Allreduce(&sumStateChanges, &globalSumStateChanges, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

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