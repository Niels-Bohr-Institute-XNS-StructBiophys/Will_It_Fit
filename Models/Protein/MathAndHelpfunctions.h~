double Mazur(double R, double L, double b, double t, double s, double expan){
    double eps, expo, arg, eexp;
    double E2, Dxy2, Sigma2;
    eps=1.-2./t;
    expo=exp(-2.*L/b);

    // Average distance squared
    E2=b*L-0.5*b*b*(1.-expo);
    // Excluded volume expansion included
    Dxy2=E2*expan;
    Sigma2=Dxy2*tgamma((s+1)/t)/tgamma((s+3.)/t);
    arg=pow( R/sqrt(fabs(Sigma2)) , t);
    if(arg >= 30){
        eexp=0;
    }
    else{
        eexp=exp(-arg);
    }
    return t/tgamma((s+1)/t) * pow(Sigma2,-0.5*(s+1)) * pow( R, s) *eexp;
}

double SemiflexibleChainAmplitude(double q, double L, double b){
    double Result;
    double Nsegments=L/b;
    double s,t,n1,n2, iS, kS, eS, iT, kT, eT, epsilon, expan;
    double sigma2,  sinqr, mazur;
    int i,j, Lsteps=40;
    double dL, R, omega=0, omega_norm=0;
    iS=2.7; //Parameters[11];//2.3;
    eS=0.6; //Parameters[12];//1.03;
    kS= 3.0; //Parameters[13];//3.0;

    s=kS*pow(iS/Nsegments,eS);
    if(s<kS)
        s=kS;
    iT=8.7; //Parameters[14];//10.0;
    eT=0.4; //Parameters[15];//0.5;
    kT=3.4; //Parameters[16];//2.43;
    t=kT*pow(iT/Nsegments,eT);
    if(t<kT)
        t=kT;
    n1= 2.42; //1.66; //2.42;
    n2= 6.90; //3.86; //6.90;
    epsilon= 0.192; //0.192;
    expan = pow(  1.+pow(Nsegments/n1,2)+ pow(Nsegments/n2,3)  ,  epsilon/3.  );
    dL=L/Lsteps;
    for (i=0;i<Lsteps;i++){
        R=(i+0.5)*dL;
        sinqr=sin(q*R)/(q*R);
        mazur=Mazur(R, L, b, t, s, expan);
        omega+=sinqr*mazur;
        omega_norm+=mazur;
    }
    //printf("%.4lf, %4lf, %.4lf\n", s, t,  tgamma((s+1)/t)/tgamma((s+3.)/t));

    return omega/omega_norm;

}

double SemiflexibleChainAmp(int NL, double q, double L, double b){
  int i = 0;
  if(NL<q*100)
      NL=(int) q*100;
  if(NL>1500)
      NL=1500;
  double ll,dL=L/NL;
  double psi, CoilSelf=0, CoilAmp=0;
  // intensity calculation
  for(i=0;i<NL;i++){
      ll=(i+0.5)*dL;
      psi = SemiflexibleChainAmplitude(q, ll, b);
      CoilSelf += 2.* psi * (L-ll)/L /NL;
      CoilAmp  += psi /NL;
  }
  return CoilAmp;
}

double SemiflexibleChainSelf(int NL, double q, double L, double b){
  double l,dL=L/NL;
  double psi, CoilSelf=0, CoilAmp=0;
  int i;
  for(i=0;i<NL;i++){
      l=(i+0.5)*dL;
      psi = SemiflexibleChainAmplitude(q, l, b);
      CoilSelf += 2.* psi * (L-l)/L /NL;
  }
  return CoilSelf;

}

double Hammouda(double q,double Rg){
    double X=q*q*Rg*Rg;
    if(X<1e-3)
        return 1.0;
    else
        return (1-exp(-X))/X;
}

double Debye(double q,double Rg){
    double X=q*q*Rg*Rg;
    if(X<1e-3)
        return 1.0;
    else
        return 2.*(exp(-X)-1.+X)/(X*X);
}

int Sign(double x) {
    int Result;

    if (x < 0.0) {
        Result = -1;
    } else if (x > 0.0) {
        Result = 1;
    } else {
        Result = 0;
	}

    return Result;
}

double complex PolarComplexNumber(double Radius, double Phi)
{
    if (Phi == 0)
        return Radius + I * 0.0;
    else
        return Radius * (cos(Phi) + I * sin(Phi));
}

void AddScatteringFromResidue(double complex **Beta, double q, struct Residue CurrentResidue, double Contrast, double ScatteringLengthDensityOfSolvent, double DeltaB) {
    // Søren Kynde, 2012 (rewritten by Martin Cramer Pedersen, 2015)
    // This function adds the scattering of a residue to the scattering amplitude expandended by the spherical harmonic coefficients Beta_lm.
	int l;
	int m;

    double x;
    double y;
    double z;

	double Radius;
	double Theta;
	double Phi;

    double Legendre[NumberOfHarmonics + 1];
    double Bessel[NumberOfHarmonics + 1];
    double ScatteringLengthOfResidue;
	double ScatteringLengthOfDisplacedSolvent = CurrentResidue.Volume * ScatteringLengthDensityOfSolvent;
	double ExcessScatteringLength;
    if (Contrast < 0.0) {
       ScatteringLengthOfResidue = CurrentResidue.XRayScatteringLength;
       ExcessScatteringLength    = ScatteringLengthOfResidue - ScatteringLengthOfDisplacedSolvent;

       //x = CurrentResidue.xVolume ; 
       x = (CurrentResidue.xXRayScattering * ScatteringLengthOfResidue - CurrentResidue.xVolume * ScatteringLengthOfDisplacedSolvent) / ExcessScatteringLength;
       //y = CurrentResidue.yVolume ; 
       y = (CurrentResidue.yXRayScattering * ScatteringLengthOfResidue - CurrentResidue.yVolume * ScatteringLengthOfDisplacedSolvent) / ExcessScatteringLength;
       // z = CurrentResidue.zVolume ;
       z = (CurrentResidue.zXRayScattering * ScatteringLengthOfResidue - CurrentResidue.zVolume * ScatteringLengthOfDisplacedSolvent) / ExcessScatteringLength;
    } else {
       ScatteringLengthOfResidue = CurrentResidue.NeutronScatteringLength;
       ExcessScatteringLength    = ScatteringLengthOfResidue - ScatteringLengthOfDisplacedSolvent;

       //x = CurrentResidue.xVolume ; 
       x = (CurrentResidue.xNeutronScattering * ScatteringLengthOfResidue - CurrentResidue.xVolume * ScatteringLengthOfDisplacedSolvent) / ExcessScatteringLength;
       //y = CurrentResidue.yVolume ; 
       y = (CurrentResidue.yNeutronScattering * ScatteringLengthOfResidue - CurrentResidue.yVolume * ScatteringLengthOfDisplacedSolvent) / ExcessScatteringLength;
       // z = CurrentResidue.zVolume ; 
       z = (CurrentResidue.zNeutronScattering * ScatteringLengthOfResidue - CurrentResidue.zVolume * ScatteringLengthOfDisplacedSolvent) / ExcessScatteringLength;
    }
    //printf("%s, %0.3e, %f ,%f, %f\n",ScatteringLengthOfResidue,CurrentResidue.Name,x,y,z);
    Radius = sqrt(pow(x, 2) + pow(y, 2) + pow(z, 2));
    Theta  = acos(z / Radius);
    Phi    = acos(x / (Radius * sin(Theta))) * Sign(y);


	// Calculate spherical Bessel functions for l = 0 to NumberOfHarmonics

    gsl_sf_bessel_jl_array(NumberOfHarmonics, q * Radius, Bessel);

	// Calculate Legendre polynomials P_l(cos(theta)) of degree l = m to NumberOfHarmonics - store the values in Legendre[m], Legendre[m + 1], ..., Legendre[NumberOfHarmonics]
    for (m = 0; m < NumberOfHarmonics + 1; m++) {
        gsl_sf_legendre_sphPlm_array(NumberOfHarmonics, m, cos(Theta), &Legendre[m]);

        for (l = m; l < NumberOfHarmonics + 1; l++) {
            Beta[l][m] += sqrt(DeltaB)*sqrt(4.0 * M_PI) * cpow(I, l) * ExcessScatteringLength * Bessel[l] * Legendre[l] * PolarComplexNumber(1.0, -m * Phi);
        }
    }
}

void AddScatteringFromSolvent(double complex **Beta, double q, struct Residue CurrentResidue, double Contrast, double ScatteringLengthDensityOfSolvent, double DeltaB) {
    // Søren Kynde, 2012 (rewritten by Martin Cramer Pedersen, 2015)
    // This function adds the scattering of a residue to the scattering amplitude expandended by the spherical harmonic coefficients Beta_lm.
	int l;
	int m;

    double x;
    double y;
    double z;

	double Radius;
	double Theta;
	double Phi;

    double Legendre[NumberOfHarmonics + 1];
    double Bessel[NumberOfHarmonics + 1];
    double ScatteringLengthOfResidue;
    double AverageVolume = 4.*M_PI * pow(1.28,3.) / 3.  ; // From Crysol paper 

    double ScatteringLengthOfDisplacedSolvent = CurrentResidue.Volume * ScatteringLengthDensityOfSolvent;
	double ExcessScatteringLength;
    if (Contrast < 0.0) {
       ScatteringLengthOfResidue = CurrentResidue.XRayScatteringLength;
       ExcessScatteringLength    = ScatteringLengthOfDisplacedSolvent;

       x = CurrentResidue.xVolume ; //(CurrentResidue.xXRayScattering * ScatteringLengthOfResidue - CurrentResidue.xVolume * ScatteringLengthOfDisplacedSolvent) / ExcessScatteringLength;
       y = CurrentResidue.yVolume ; //(CurrentResidue.yXRayScattering * ScatteringLengthOfResidue - CurrentResidue.yVolume * ScatteringLengthOfDisplacedSolvent) / ExcessScatteringLength;
       z = CurrentResidue.zVolume ; //(CurrentResidue.zXRayScattering * ScatteringLengthOfResidue - CurrentResidue.zVolume * ScatteringLengthOfDisplacedSolvent) / ExcessScatteringLength;
    } else {
       ScatteringLengthOfResidue = CurrentResidue.NeutronScatteringLength;
       ExcessScatteringLength    = ScatteringLengthOfDisplacedSolvent;

       x = CurrentResidue.xVolume ; //(CurrentResidue.xNeutronScattering * ScatteringLengthOfResidue - CurrentResidue.xVolume * ScatteringLengthOfDisplacedSolvent) / ExcessScatteringLength;
       y = CurrentResidue.yVolume ; //(CurrentResidue.yNeutronScattering * ScatteringLengthOfResidue - CurrentResidue.yVolume * ScatteringLengthOfDisplacedSolvent) / ExcessScatteringLength;
       z = CurrentResidue.zVolume ; //(CurrentResidue.zNeutronScattering * ScatteringLengthOfResidue - CurrentResidue.zVolume * ScatteringLengthOfDisplacedSolvent) / ExcessScatteringLength;
    }

    Radius = sqrt(pow(x, 2) + pow(y, 2) + pow(z, 2));
    Theta  = acos(z / Radius);
    Phi    = acos(x / (Radius * sin(Theta))) * Sign(y);


	// Calculate spherical Bessel functions for l = 0 to NumberOfHarmonics

    gsl_sf_bessel_jl_array(NumberOfHarmonics, q * Radius, Bessel);

	// Calculate Legendre polynomials P_l(cos(theta)) of degree l = m to NumberOfHarmonics - store the values in Legendre[m], Legendre[m + 1], ..., Legendre[NumberOfHarmonics]
    for (m = 0; m < NumberOfHarmonics + 1; m++) {
        gsl_sf_legendre_sphPlm_array(NumberOfHarmonics, m, cos(Theta), &Legendre[m]);

        for (l = m; l < NumberOfHarmonics + 1; l++) {
            Beta[l][m] += sqrt(DeltaB)*sqrt(4.0 * M_PI) * cpow(I, l) * ExcessScatteringLength * exp(-1.*pow(q,2.)*pow(AverageVolume, 2./3.)/(4.*M_PI))* Bessel[l] * Legendre[l] * PolarComplexNumber(1.0, -m * Phi);
        }
    }
}

void CopyResidue(struct Residue * Original, struct Residue * Copy) {
	Copy->xVolume = Original->xVolume;
	Copy->yVolume = Original->yVolume;
	Copy->zVolume = Original->zVolume;

	Copy->xXRayScattering = Original->xXRayScattering;
	Copy->yXRayScattering = Original->yXRayScattering;
	Copy->zXRayScattering = Original->zXRayScattering;

	Copy->xNeutronScattering = Original->xNeutronScattering;
	Copy->yNeutronScattering = Original->yNeutronScattering;
	Copy->zNeutronScattering = Original->zNeutronScattering;

	Copy->XRayScatteringLength    = Original->XRayScatteringLength;
	Copy->NeutronScatteringLength = Original->NeutronScatteringLength;

	Copy->Volume    = Original->Volume;
	Copy->Name[0]   = Original->Name[0];
	Copy->Name[1]   = Original->Name[1];
	Copy->Name[2]   = Original->Name[2];
	Copy->Name[3]   = Original->Name[3];
//	printf("%s\n", Copy->Name);
	Copy->ResidueID = Original->ResidueID;
}

int sign(double x){
    int result;
    if(x<0)
        result=-1;
    else if(x>0)
        result=1;
    else
        result=0;
    return result;
}

double complex pol(double r, double phi)
{
    if(phi==0)
        return r+I*0;
    else
        return r*(cos(phi)+I*sin(phi));
}
double complex **ComplexArray(dim1,dim2)
{
    int i,ii;
    double complex ** arr;
    arr = (double complex**)malloc(dim1*sizeof(double complex*));
    for (i = 0; i< dim1; i++) {
        arr[i] = (double complex*) malloc((dim2)*sizeof(double complex));
        for (ii = 0; ii < dim2; ii++) {
            arr[i][ii] = 0;
        }
    }
    return arr;
}

void FreeComplexArray(double complex **alpha, int dim1, int dim2)
{
    int i;
    for(i=0;i<dim1;i++){
        free(alpha[i]);
    }
    free(alpha);
}

void AddScatteringFromSemiflexibleTag(double complex **beta, double Q, double x_tag, double y_tag, double z_tag,
       double DeltaB, int NL, double L, double b){
    //Søren Kynde 2012
    //This function adds the scattering of a gaussian random coil to the scattering amplitude
    //expandended by the spherical harmonic coefficients beta_lm.

    int l,m;
    double legendre[NumberOfHarmonics+1];
    double bessel[NumberOfHarmonics+1];
    double SQRT4pi=sqrt(4*M_PI);
//    double X=pow((Q*Rg_tag),2);


    double PSI = SemiflexibleChainAmp(NL, Q, L, b);//(1-exp(-X))/X;

    double rn=sqrt(pow(x_tag,2)+pow(y_tag,2)+pow(z_tag,2));
    double thetan=acos(z_tag/rn);
    double phin=acos(x_tag/(rn*sin(thetan)))*sign(y_tag);


    gsl_sf_bessel_jl_array(NumberOfHarmonics,Q*rn,bessel); // Calculate spherical bessel functions for l=0,..,Nh
    for(m=0;m<=NumberOfHarmonics;m++){
        gsl_sf_legendre_sphPlm_array(NumberOfHarmonics,m,cos(thetan),&legendre[m]); //Calculate legendre polynomials P_l(cos(theta)) of degree l=m,..., up to Nh
                                                                     //Store the values in legendre[m],legendre[m+1],...,legendre[Nh]
        for(l=m;l<NumberOfHarmonics+1;l++){
            beta[l][m]+=DeltaB*PSI*SQRT4pi*cpow(I,l)*bessel[l]*legendre[l]*pol(1,-m*phin);
        }
    }
}


double Si(double x){
    double t;
    int nsteps=100, i;
    double dt=x/(nsteps);
    double result=0;
    for(i=0;i<nsteps;i++){
        t=(i+0.5)*dt;
        result+=1./t*sin(t)*dt;
    }
    return result;
}

double RPAStructureFactor(double Itensity,double v){
  return 1./(1.+Itensity*v);

}
