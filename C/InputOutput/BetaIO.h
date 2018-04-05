char * QToFileName(double q, bool prot)
{
  char name[50];
  char folder[50] = "beta/";
  char *filename = malloc (sizeof (char) * 50);
  sprintf(name, "%lf", q);
  name[0] = 'q';
  name[1] = '_';
  strcat(folder, name);
  if(prot)
  {
    strcat(folder, ".betaprot");
  }
  else{
    strcat(folder, ".beta");
  }
  return filename;
}

void WriteComplexArray(double q, double complex** alpha,int m,int n, bool prot){
  int i,j;
  FILE *wt;
  double complex curr;
  char * fileName;
  fileName = QToFileName(q, prot);
  //Code to do stuff with q goes here.
  printf("Write: %s\n", fileName);
  wt = fopen(fileName, "w");
  fprintf(wt, "%f\n", q);
  fprintf(wt, "%d:%d\n",m,n);
  for (i = 0; i < m; ++i)
  {
    for (j = 0; j < n; ++j)
    {
      curr = alpha[i][j];
      //curr = i + j*I;
      fprintf(wt, "%f+%f\n",crealf(curr), cimagf(curr));
    }
  }
  close(wt);
}

double complex **ReadComplexArray(double q, double complex **alpha, bool prot){
  //Outcommented code was used for testing.
  int i, j, m, n, r, ri, ii;
  double real, img;
  double complex curr;
  char * fileName;
  fileName = QToFileName(q, prot);
  printf("Read: %s\n", fileName);
  FILE *rt;

  rt = fopen(fileName, "r");
  r = fscanf(rt, "%lf", &q);
  r = fscanf(rt, "%d:%d", &m,&n);
  //printf("%f %d %d\n", q, m, n);
  for (i = 0; i < m; ++i)
  {
    for (j = 0; j < n; ++j)
    {
      r = fscanf(rt, "%lf+%lf", &real, &img);
      ri = (int) real;
      ii = (int) img;
      /*printf("i: %d j: %d Val:%lf*%lf, r: %d\n", i,j,real,img, r);
      if (ri != i || ii != j)
      {
        printf("i: %d j: %d Val:%lf*%lf, r: %d\n", i,j,real,img, r);
      }*/
      curr = real*img*I;
      alpha[i][j] = curr;
    }
  }
  close(rt);
  return alpha;
}