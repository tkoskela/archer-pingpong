#include <stdio.h>
#include<stdlib.h>
#include <mpi.h>

#define K 1024
#define M K*K

#define BOSS 0

#define STANDARD
#define SYNC

#undef BUFFERED
#undef STRIDED

#define BUFFSIZ M
#define START 1
#define TAG 100
#define SEND 1
#define REPEAT 1000


int buffer[BUFFSIZ];
int buffer2[BUFFSIZ];
int mpi_buffer[BUFFSIZ];
MPI_Datatype stride;
void check_error(int , char *);
void ping(int peer);
void pong(int peer);
void printlocation();

int main(int argc, char **argv)
{
  int error;
  int npe, me;
  int i;

  
  error = MPI_Init(&argc,&argv);
  check_error(error,"MPI_Init");


  for(i=0;i<BUFFSIZ;i++){
    buffer[i]=i;
    buffer2[i]=-1;
  }

  error = MPI_Comm_size(MPI_COMM_WORLD, &npe);
  check_error(error,"MPI_Comm_size");

  if (npe < 2)
    {
      printf("Error! Must run on at least 2 processes\n");

      error = MPI_Finalize();
      check_error(error,"finalise");

      exit(1);
    }

  error = MPI_Comm_rank(MPI_COMM_WORLD, &me);
  check_error(error,"MPI_Comm_rank");

  error = MPI_Buffer_attach(&mpi_buffer,sizeof(mpi_buffer));
  check_error(error,"MPI_Buffer_attach");

  error = MPI_Barrier(MPI_COMM_WORLD);
  check_error(error,"MPI_Barrier");

  if( me == BOSS ){
    printf("Running on %d processes: sending messages between ranks %d and %d\n", npe, 0, npe-1);
    printlocation();
    ping(npe-1);
  }else if(me == (npe-1)){
    printlocation();
    pong(BOSS);
  }

  error = MPI_Finalize();
  check_error(error,"finalise");

  exit(0);
}

void check_error(int code, char *string)
{
  char exp[MPI_MAX_ERROR_STRING];
  int len;

/*  printf("trace %s\n",string); */
  if( code != MPI_SUCCESS )
  {
     MPI_Error_string(code,exp,&len);
     fprintf(stderr,"MPI error [%s] %d %s\n",exp,code,string);
     exit(1);
  }
}


void ping(int peer)
{
  int error;
  MPI_Status stat;
  int decade;
  int i, size;
  double start, stop, bw, elapsed;
  FILE *out;
  FILE *bw_out;

#ifdef STANDARD
  printf("Standard Send\n");
  out=fopen("standard_time.plot","w");
  bw_out=fopen("standard_bandwidth.plot","w");
  for(decade=START;decade<=BUFFSIZ;decade *= 10)
  for(size=decade;size<=BUFFSIZ && size < (10*decade) ;size+=decade)
  {
    for(i=0;i<BUFFSIZ;i++){
       buffer[i]=i;
       buffer2[i]=-1;
    }
    start = MPI_Wtime();
    for(i=0; i<REPEAT; i++)
    {
  
      error= MPI_Send(buffer,size,MPI_INT,peer,TAG,MPI_COMM_WORLD);
      check_error(error,"send");
    
      error = MPI_Recv(buffer2,BUFFSIZ,MPI_INT,peer,TAG,MPI_COMM_WORLD,&stat);
      check_error(error,"recv");
    }
    stop = MPI_Wtime();
    for(i=0;i<size;i++){
       if(buffer[i] != buffer2[i]){
         fprintf(stderr,"Data copy Error %d %d != %d\n",i,buffer[i],buffer2[i]);
         exit(1);
       }
    }
    elapsed = stop-start;
    bw = 2.0*REPEAT*size*sizeof(int)/(1.0*elapsed*M); 
    printf("%d ping-pongs of %d integers took %g seconds,%g each, %g Mb/sec\n",
      REPEAT,size,elapsed, elapsed/REPEAT,bw);
    fprintf(out,"%d %g\n",size*sizeof(int),elapsed/REPEAT);
    fprintf(bw_out,"%d %g\n",size*sizeof(int),bw);

  }
  fclose(out);
  fclose(bw_out);
#endif
#ifdef BUFFERED
  out=fopen("buffered_time.plot","w");
  bw_out=fopen("buffered_bandwidth.plot","w");
  printf("Buffered Send\n");
  for(decade=START;decade<=BUFFSIZ;decade *= 10)
  for(size=decade;size<=BUFFSIZ && size < (10*decade) ;size+=decade)
  {
    for(i=0;i<BUFFSIZ;i++){
       buffer[i]=i;
       buffer2[i]=-1;
    }
    start = MPI_Wtime();
    for(i=0; i<REPEAT; i++)
    {
  
      error= MPI_Bsend(buffer,size,MPI_INT,peer,TAG,MPI_COMM_WORLD);
      check_error(error,"send");
    
      error = MPI_Recv(buffer2,BUFFSIZ,MPI_INT,peer,TAG,MPI_COMM_WORLD,&stat);
      check_error(error,"recv");
    }
    stop = MPI_Wtime();
    for(i=0;i<size;i++){
       if(buffer[i] != buffer2[i]){
         fprintf(stderr,"Data copy Error %d %d != %d\n",i,buffer[i],buffer2[i]);
         exit(1);
       }
    }
    elapsed = stop-start;
    bw = 2.0*REPEAT*size*sizeof(int)/(1.0*elapsed*M); 
    printf("%d ping-pongs of %d integers took %g seconds,%g each, %g Mb/sec\n",
      REPEAT,size,elapsed, elapsed/REPEAT,bw);
    fprintf(out,"%d %g\n",size*sizeof(int),elapsed/REPEAT);
    fprintf(bw_out,"%d %g\n",size*sizeof(int),bw);

  }
  fclose(out);
  fclose(bw_out);
#endif
#ifdef SYNC
  printf("Sync Send\n");
  out=fopen("sync_time.plot","w");
  bw_out=fopen("sync_bandwidth.plot","w");
  for(decade=START;decade<=BUFFSIZ;decade *= 10)
  for(size=decade;size<=BUFFSIZ && size < (10*decade) ;size+=decade)
  {
    for(i=0;i<BUFFSIZ;i++){
       buffer[i]=i;
       buffer2[i]=-1;
    }
    start = MPI_Wtime();
    for(i=0; i<REPEAT; i++)
    {
  
      error= MPI_Ssend(buffer,size,MPI_INT,peer,TAG,MPI_COMM_WORLD);
      check_error(error,"send");
    
      error = MPI_Recv(buffer2,BUFFSIZ,MPI_INT,peer,TAG,MPI_COMM_WORLD,&stat);
      check_error(error,"recv");
    }
    stop = MPI_Wtime();
    elapsed = stop-start;
    bw = 2.0*REPEAT*size*sizeof(int)/(1.0*elapsed*M); 
    printf("%d ping-pongs of %d integers took %g seconds,%g each, %g Mb/sec\n",
      REPEAT,size,elapsed, elapsed/REPEAT,bw);
    fprintf(out,"%d %g\n",size*sizeof(int),elapsed/REPEAT);
    fprintf(bw_out,"%d %g\n",size*sizeof(int),bw);
    for(i=0;i<size;i++){
       if(buffer[i] != buffer2[i]){
         fprintf(stderr,"Data copy Error %d %d != %d\n",i,buffer[i],buffer2[i]);
         exit(1);
       }
    }

  }
  fclose(out);
  fclose(bw_out);
#endif
#ifdef STRIDED
  printf("Strided Send\n");
  out=fopen("strided_time.plot","w");
  bw_out=fopen("strided_bandwidth.plot","w");
  for(decade=START;decade<=BUFFSIZ/2;decade *= 10)
  for(size=decade;size<=BUFFSIZ/2 && size < (10*decade) ;size+=decade)
  {
    for(i=0;i<BUFFSIZ;i++){
       buffer[i]=i;
       buffer2[i]=-1;
    }
    error= MPI_Type_vector(size,1,2,MPI_INT,&stride);
    check_error(error,"Make stride type");
    error = MPI_Type_commit(&stride);
    check_error(error,"commit stride type");
    start = MPI_Wtime();
    for(i=0; i<REPEAT; i++)
    {
  
      error= MPI_Send(buffer,1,stride,peer,TAG,MPI_COMM_WORLD);
      check_error(error,"send");
    
      error = MPI_Recv(buffer2,1,stride,peer,TAG,MPI_COMM_WORLD,&stat);
      check_error(error,"recv");
    }
    stop = MPI_Wtime();
    error = MPI_Type_free(&stride);
    check_error(error,"free stride type");
    for(i=0;i<size;i++){
       if(buffer[2*i] != buffer2[2*i]){
         fprintf(stderr,"Data copy Error %d %d != %d\n",2*i,buffer[2*i],buffer2[2*i]);
         exit(1);
       }
       if(-1 != buffer2[(2*i)+1]){
         fprintf(stderr,"Data copy Error %d %d != %d\n",2*i+1,-1,buffer2[2*i+1]);
         exit(1);
       }
    }
    elapsed = stop-start;
    bw = 2.0*REPEAT*size*sizeof(int)/(1.0*elapsed*M); 
    printf("%d ping-pongs of %d integers took %g seconds,%g each, %g Mb/sec\n",
      REPEAT,size,elapsed, elapsed/REPEAT,bw);
    fprintf(out,"%d %g\n",size*sizeof(int),elapsed/REPEAT);
    fprintf(bw_out,"%d %g\n",size*sizeof(int),bw);

  }
  fclose(out);
  fclose(bw_out);
#endif
}

void pong(int peer)
{
  int error;
  MPI_Status stat;
  int count;
  int running=1;
  int i;
  int size,decade;

#ifdef STANDARD
  for(decade=START;decade<=BUFFSIZ;decade *= 10)
  for(size=decade;size<=BUFFSIZ && size < (10*decade) ;size+=decade)
  {
    for(i=0; i<REPEAT; i++)
    {
      error = MPI_Recv(buffer,BUFFSIZ,MPI_INT,peer,TAG,MPI_COMM_WORLD,&stat);
      check_error(error,"recv");

      error= MPI_Send(buffer,size,MPI_INT,peer,TAG,MPI_COMM_WORLD);
      check_error(error,"send");
    }
  }
#endif
#ifdef BUFFERED
  for(decade=START;decade<=BUFFSIZ;decade *= 10)
  for(size=decade;size<=BUFFSIZ && size < (10*decade) ;size+=decade)
  {
    for(i=0; i<REPEAT; i++)
    {
      error = MPI_Recv(buffer,BUFFSIZ,MPI_INT,peer,TAG,MPI_COMM_WORLD,&stat);
      check_error(error,"recv");

      error= MPI_Bsend(buffer,size,MPI_INT,peer,TAG,MPI_COMM_WORLD);
      check_error(error,"send");
    }
  }
#endif
#ifdef SYNC
  for(decade=START;decade<=BUFFSIZ;decade *= 10)
  for(size=decade;size<=BUFFSIZ && size < (10*decade) ;size+=decade)
  {
    for(i=0; i<REPEAT; i++)
    {
      error = MPI_Recv(buffer,BUFFSIZ,MPI_INT,peer,TAG,MPI_COMM_WORLD,&stat);
      check_error(error,"recv");

      error= MPI_Ssend(buffer,size,MPI_INT,peer,TAG,MPI_COMM_WORLD);
      check_error(error,"send");
    }
  }
#endif
#ifdef STRIDED
  for(decade=START;decade<=BUFFSIZ/2;decade *= 10)
  for(size=decade;size<=BUFFSIZ/2 && size < (10*decade) ;size+=decade)
  {
    error= MPI_Type_vector(size,1,2,MPI_INT,&stride);
    check_error(error,"Make stride type");
    error = MPI_Type_commit(&stride);
    check_error(error,"commit stride type");
    for(i=0; i<REPEAT; i++)
    {
      error = MPI_Recv(buffer,1,stride,peer,TAG,MPI_COMM_WORLD,&stat);
      check_error(error,"recv");

      error= MPI_Send(buffer,1,stride,peer,TAG,MPI_COMM_WORLD);
      check_error(error,"send");
    }
    error = MPI_Type_free(&stride);
    check_error(error,"free stride type");
  }
#endif
}

