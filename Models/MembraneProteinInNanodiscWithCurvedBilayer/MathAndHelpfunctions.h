void SetRotationMatrix(double RM[3][3], double rz1, double rx, double rz2){
    /*Rotation about z-axis with rz1 followed by rotation about x-axis by rx
     * followed by rotation about z-axis by rz2.
     * Positive angles give right handed rotations.*/
    RM[0][0]=cos(rz1)*cos(rz2)-sin(rz1)*cos(rx)*sin(rz2); RM[0][1]=-sin(rz1)*cos(rz2)-cos(rz1)*cos(rx)*sin(rz2); RM[0][2]= sin(rx)*sin(rz2);
    RM[1][0]=cos(rz1)*sin(rz2)+sin(rz1)*cos(rx)*cos(rz2); RM[1][1]=-sin(rz1)*sin(rz2)+cos(rz1)*cos(rx)*cos(rz2); RM[1][2]=-sin(rx)*cos(rz2);
    RM[2][0]=sin(rz1)*sin(rx)                           ; RM[2][1]= cos(rz1)*sin(rx)                           ; RM[2][2]= cos(rx)         ;
}

void SetTranslationVector(double TV[3], double s, double theta, double z){
    TV[0] = s*cos(theta);
    TV[1] = s*sin(theta);
    TV[2] = z;
}


void Orient(double* x, double* y, double* z, double rotation[3][3], double translation[3]){
    double xinit=*x;
    double yinit=*y;
    double zinit=*z;
    *x=rotation[0][0]*xinit+rotation[0][1]*yinit+rotation[0][2]*zinit+translation[0];
    *y=rotation[1][0]*xinit+rotation[1][1]*yinit+rotation[1][2]*zinit+translation[1];
    *z=rotation[2][0]*xinit+rotation[2][1]*yinit+rotation[2][2]*zinit+translation[2];
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

void FreeComplexArray(double complex **alpha, int dim1, int dim2){
    int i;
    for(i=0;i<dim1;i++){
        free(alpha[i]);
    }
    free(alpha);
}

void CopyResidue(struct Residue * Original, struct Residue * Copy){
      Copy->xVolume=Original->xVolume;
      Copy->yVolume=Original->yVolume;
      Copy->zVolume=Original->zVolume;
      
      Copy->xXRayScattering=Original->xXRayScattering;
      Copy->yXRayScattering=Original->yXRayScattering;
      Copy->zXRayScattering=Original->zXRayScattering;
      
      Copy->xNeutronScattering=Original->xNeutronScattering;
      Copy->yNeutronScattering=Original->yNeutronScattering;
      Copy->zNeutronScattering=Original->zNeutronScattering;
    
      Copy->XRayScatteringLength = Original->XRayScatteringLength;
      Copy->NeutronScatteringLength = Original->NeutronScatteringLength;
      Copy->Volume = Original->Volume;

     Copy->Name[0] = Original->Name[0];
     Copy->Name[1] = Original->Name[1];
     Copy->Name[2] = Original->Name[2];
     Copy->ResidueID = Original->ResidueID;
}

double Sinc(double x){
    double result;
    if(x==0)
        result=1.0;
    else
        result=sin(x)/x;
    return result;
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

double expand(double complex **alpha,double complex **F)
{
//Søren Kynde 2011
//This function calculates the coefficients alpha_lm of the spherical harmonics
//expansion of an analytical form factor F(theta,phi) (phi is the azimuthal angle)
    double theta,phi;
    double thetastep=M_PI/(Ntheta), phistep=2*M_PI/(Nphi);
    int ntheta,nphi;
    double complex fm[Ntheta]={0};
    double complex phase[Nh+1][Nphi];
    int l,m,i,j;
    double Int=0;
    double sinth[Ntheta];
    double w[Ntheta]={0};
    double legendre[Ntheta][Nh+1];
    size_t LegendreSize = gsl_sf_legendre_array_n(Nh);
    int LegendreIndex = 0;
    double Legendre[LegendreSize];
    double alph = 0.;
    const int LegendreMode = 0;

    for(i=0;i<Ntheta;i++){
        fm[i]=0;
    }

    for(i=0;i<Ntheta;i++){
        w[i]=0;
    }
    for(j=0;j<Ntheta;j++){ //calculate weights dtheta for theta integral
        theta=(j+.5)*thetastep;
        for(l=0;l<Ntheta/2;l++){
            w[j]+=(double) 2./(Ntheta/2)*1./(2*l+1)*sin((2*l+1)*theta);
        }
    }


    for (ntheta=0;ntheta<Ntheta;ntheta++){
        theta = (ntheta+.5)*thetastep; //Theta integration
        fm[ntheta] = 0.;
        sinth[ntheta] = sin(theta);
        gsl_sf_legendre_array(LegendreMode, Nh, cos(theta), Legendre);//Calculate all Plm(cos(th)) for l=m ... to l=N;
        for(m=0;m<Nh+1;m+=skip){
            for(nphi=0;nphi<Nphi;nphi++){//Phi integration
                phi=phistep*(nphi+.5);
                if(ntheta==0)
                    phase[m][nphi]=pol(phistep,-m*phi);
                fm[ntheta]+= F[ntheta][nphi]*phase[m][nphi];  //fm(theta)= int_0^2pi [ F(theta,phi) exp(-m*phi) dphi ]
            }
            for(l=m;l<=Nh;l+=skip){ //For disc symmetry only even harmonics contribute
                LegendreIndex = gsl_sf_legendre_array_index(l, m);
                alph = 1/sqrt(4*M_PI)*Legendre[LegendreIndex]*w[ntheta]*sinth[ntheta]*fm[ntheta];
                alpha[l][m] += alph;//1/sqrt(4*M_PI)*Legendre[legendreIndex]*w[ntheta]*sinth[ntheta]*fm[ntheta];
            }      
        }
    }
    return Int;
}

bool ResidueIsInMethylLayer(double zn, double Hmethyl){
    bool result;
    result = (fabs(zn)<Hmethyl/2.);
    return result;
}
   
bool ResidueIsInAlkylLayer(double xn, double yn, double zn, double aEndcaps, double bEndcaps, double cEndcaps,
        double ShiftOfEndcaps, double Hcore){
    bool result;
    result = (( zn>0 && pow(xn/aEndcaps,2)+pow(yn/bEndcaps,2)+pow((zn-(Hcore/2.+ShiftOfEndcaps))/cEndcaps,2)<1)
            || (zn<0 && pow(xn/aEndcaps,2)+pow(yn/bEndcaps,2)+pow((zn+(Hcore/2.+ShiftOfEndcaps))/cEndcaps,2)<1)
            || fabs(zn)<Hcore/2.);
    return result;
}

bool ResidueIsInLipidLayer(double xn, double yn, double zn, double aEndcaps, double bEndcaps, double cEndcaps,
        double ShiftOfEndcaps, double Hlipid){
    bool result;
    result = (( zn>0 && pow(xn/aEndcaps,2)+pow(yn/bEndcaps,2)+pow((zn-(Hlipid/2.+ShiftOfEndcaps))/cEndcaps,2)<1)
            || (zn <0 && pow(xn/aEndcaps,2)+pow(yn/bEndcaps,2)+pow((zn+(Hlipid/2.+ShiftOfEndcaps))/cEndcaps,2)<1)
            || fabs(zn)<Hlipid/2.);
    return result;
}

void AddScatteringFromResidue(double complex **beta, double Q, struct Residue residue, double Contrast,
       double rho_head,double rho_alkyl,double rho_methyl,double rho_belt,double rho_solvent,
       double Hlipid,double Hcore,double Hmethyl,double Hbelt, double CV_protein,
       double radius_major, double radius_minor, double cEndcaps, double ScaleFactorOfEndcaps){
    //Søren Kynde 2012
    //This function adds the scattering of a residue to the scattering amplitude
    //expandended by the spherical harmonic coefficients beta_lm.
    double xn;
    double yn;
    double zn;
    int l,m;
    double Volume=residue.Volume*CV_protein;
    //Variables used for calculating Plm functions
    int LegendreSize =  gsl_sf_legendre_array_n(Nh);
    int LegendreIndex;
    const int LegendreMode = 0;
    double Legendre[LegendreSize];

    double bessel[Nh+1];
    double rho_bg=rho_solvent;
    double SQRT4pi=sqrt(4*M_PI);
    double aEndcaps=radius_major*ScaleFactorOfEndcaps;
    double bEndcaps=radius_minor*ScaleFactorOfEndcaps;
    double ShiftOfEndcaps = - cEndcaps / aEndcaps * sqrt(pow(aEndcaps, 2) - pow(radius_major, 2));
    if(ResidueIsInMethylLayer(residue.zVolume, Hmethyl))
        rho_bg=rho_methyl;
    else if(ResidueIsInAlkylLayer(residue.xVolume, residue.yVolume, residue.zVolume, aEndcaps, bEndcaps, cEndcaps, ShiftOfEndcaps, Hcore))
        rho_bg=rho_alkyl;
    else if(ResidueIsInLipidLayer(residue.xVolume, residue.yVolume, residue.zVolume, aEndcaps, bEndcaps, cEndcaps, ShiftOfEndcaps, Hlipid))
        rho_bg=rho_head;

    double Bresidue, Bsolvent, DeltaB;
    if(Contrast<0){
       Bresidue = residue.XRayScatteringLength;
       Bsolvent = residue.Volume*rho_bg*CV_protein;
       DeltaB   = Bresidue-Bsolvent;
       xn=(residue.xXRayScattering  * Bresidue - residue.xVolume*Bsolvent)/DeltaB;
       yn=(residue.yXRayScattering  * Bresidue - residue.yVolume*Bsolvent)/DeltaB;
       zn=(residue.zXRayScattering  * Bresidue - residue.zVolume*Bsolvent)/DeltaB;
    }
    else{
       Bresidue = residue.NeutronScatteringLength;
       Bsolvent = residue.Volume*rho_bg*CV_protein;
       DeltaB   = Bresidue-Bsolvent;
       xn=(residue.xNeutronScattering*Bresidue - residue.xVolume*Bsolvent)/DeltaB;
       yn=(residue.yNeutronScattering*Bresidue - residue.yVolume*Bsolvent)/DeltaB;
       zn=(residue.zNeutronScattering*Bresidue - residue.zVolume*Bsolvent)/DeltaB;
    }

    double rn=sqrt(pow(xn,2)+pow(yn,2)+pow(zn,2));
    double thetan=acos(zn/rn);
    double phin=acos(xn/(rn*sin(thetan)))*sign(yn);
    gsl_sf_bessel_jl_array(Nh,Q*rn,bessel); // Calculate spherical bessel functions for l=0,..,Nh
    gsl_sf_legendre_array(LegendreMode, Nh,cos(thetan),Legendre); //Calculate legendre polynomials P_l(cos(theta)) of degree l=m,..., up to Nh   
    for(m=0;m<=Nh;m++){                                                               
        for(l=m;l<Nh+1;l++){
            LegendreIndex = gsl_sf_legendre_array_index(l,m);
            beta[l][m]+=SQRT4pi*cpow(I,l)*DeltaB*bessel[l]*Legendre[LegendreIndex]*pol(1,-m*phin);
        }
    }
}

/// This function is used to compute the structure of a cylinder with triaxial half ellipsoids as endcaps
double PsiEllipticCylinderWithEndcaps(double q, double Alpha, double Beta, double MajorRadius, double MinorRadius, double Height, double ScaleFactorOfEndcaps, double VerticalAxisOfEndcaps)
{
    /// Declarations
    // Dummies
    int k;

    double Dummy1;
    double Dummy2;
    double Dummy3;

    double EquivalentRadiusInCylinder;
    double EquivalentRadiusInEndcaps;

    const int TSteps = 50;
    double TMin;
    double TMax;
    double TStepSize;
    double T;

    double ReturnValue;
    double SumOverT = 0.0f;

    // Assign values to arguments
    //ScaleFactorOfEndcaps = 1.5f;
    double MajorRadiusOfCurvatureOfEndcaps = MajorRadius * ScaleFactorOfEndcaps;
    double MinorRadiusOfCurvatureOfEndcaps = MinorRadius * ScaleFactorOfEndcaps;
    double ShiftOfEndcaps;


    // Caps center inside cylidner
    ShiftOfEndcaps = - VerticalAxisOfEndcaps / MinorRadiusOfCurvatureOfEndcaps * sqrt(pow(MinorRadiusOfCurvatureOfEndcaps, 2) - pow(MinorRadius, 2));

    EquivalentRadiusInCylinder = sqrt(pow(MajorRadius * cos(Beta), 2) + pow(MinorRadius * sin(Beta), 2));
    EquivalentRadiusInEndcaps  = sqrt(pow(MajorRadiusOfCurvatureOfEndcaps * cos(Beta), 2) + pow(MinorRadiusOfCurvatureOfEndcaps * sin(Beta), 2));

    TMin = - ShiftOfEndcaps / VerticalAxisOfEndcaps;
    TMax = 1.0f;
    TStepSize = (TMax - TMin) / TSteps;

    Dummy1 = pi * MajorRadius * MinorRadius * Height *
        2.0f * sin(q * Height / 2.0f * cos(Alpha)) / (q * Height / 2.0f * cos(Alpha)) *
        j1(q * EquivalentRadiusInCylinder * sin(Alpha)) / (q * EquivalentRadiusInCylinder * sin(Alpha));

    SumOverT = 0.0f;

    for (k = 0; k < TSteps; ++k) {
        T = k * TStepSize + TMin;

            Dummy2 = cos(q * cos(Alpha) * (VerticalAxisOfEndcaps * T + ShiftOfEndcaps + Height / 2.0f)) *
            (1 - pow(T, 2)) * j1(q * EquivalentRadiusInEndcaps * sin(Alpha) * sqrt(1.0f - pow(T, 2))) /
            (q * EquivalentRadiusInEndcaps * sin(Alpha) * sqrt(1.0f - pow(T, 2)));

            SumOverT += Dummy2 * TStepSize;
    }

    Dummy3 = 4.0f * pi * MajorRadiusOfCurvatureOfEndcaps * MinorRadiusOfCurvatureOfEndcaps * VerticalAxisOfEndcaps * SumOverT;
    ReturnValue=Dummy1+Dummy3;
    if(isnan(SumOverT)){
        printf("NH2\n");
    }

    /*double Volume= pi * MajorRadiusOfCurvatureOfEndcaps * MinorRadiusOfCurvatureOfEndcaps * VerticalAxisOfEndcaps * (1.-1/3.-(-ShiftOfEndcaps)/(VerticalAxisOfEndcaps)+1/3.*pow((-ShiftOfEndcaps)/(VerticalAxisOfEndcaps),3));
    Volume*=2;
    Volume+=MajorRadius*MinorRadius*pi*Height;*/ // Calculate volume of elliptical cyllinder with two endcaps.

    return ReturnValue;
}

int DiscWithEndcaps(double complex **F, double a, double b,double L,
        double ScaleFactorOfEndcaps, double VerticalAxisOfEndcaps,
        double rho, double r0, double theta0, double phi0, double Q)
{
    //This function calculates the form factor F(Q,theta,phi) of a cylinder
    //with elliptical endcaps. The
    //cylinder may be dispaced from the origin by a (real space) vector
    //(r0, theta0, phi0).

    double cosr0q=0;
    int ntheta,nphi;
    double theta,phi;
    double thetastep=M_PI/Ntheta, phistep=2*M_PI/Nphi;

    for(ntheta=0; ntheta < Ntheta; ntheta++){
        theta=(ntheta+.5)*thetastep;
        for(nphi=0; nphi<Nphi; nphi++){
            phi=(nphi+.5)*phistep;
            if( (r0!=0.) || (theta0!=0) || (phi0!=0) )
                cosr0q=sin(theta0)*sin(theta)*cos(phi0-phi)+cos(theta0)*cos(theta);
            F[ntheta][nphi]+=rho*pol(PsiEllipticCylinderWithEndcaps(Q,theta,phi,a,b,L,ScaleFactorOfEndcaps,VerticalAxisOfEndcaps) ,-r0*Q*cosr0q);
        }
    }
    return 0;
}

int GaussCoil(double complex **F, double Rg, double rho, double volume, double r0, double theta0, double phi0, double Q)
{
    //This function calculates the form factor amplitude of a Gaussian random
    //coil based on the the Hammouda formula. The center of the coil may be
    //displaced from the origin by the vector (r0,theta0,phi0).
    double r;
    double cosr0q=0;
    int ntheta,nphi;
    double sinc,theta,phi;
    double thetastep=M_PI/Ntheta, phistep=2*M_PI/Nphi;
    double X=pow((Q*Rg),2);

    for(ntheta=0; ntheta < Ntheta; ntheta++){
        theta=(ntheta+.5)*thetastep;
        for(nphi=0; nphi<Nphi; nphi++){
            phi=(nphi+.5)*phistep;
            if( (r0!=0.) || (theta0!=0) || (phi0!=0) )
                cosr0q=sin(theta0)*sin(theta)*cos(phi0-phi)+cos(theta0)*cos(theta);
            if(Q==0)
                F[ntheta][nphi]+=volume*rho*pol(1.0,0.0);
            else
                F[ntheta][nphi]+=volume*rho*pol( (1-exp(-X))/X ,-r0*Q*cosr0q);
        }
    }
    //printf("P: %g, %g, %g, %g \n",Q,Rg,rho,volume);
    //getchar();
    return 0;

}

int Discflat(double complex **F, double a, double b,double L, double rho, double r0, double theta0, double phi0, double Q)
{
    //Søren Kynde 2012
    //This function calculates the form factor F(Q,theta,phi) of a cylinder
    //with elliptical cross-section with half axes a and b and height L. The
    //cylinder may be dispaced from the origin by a (real space) vector
    //(r0, theta0, phi0).

    double r;
    double cosr0q=0;
    int ntheta,nphi;
    double sinc,theta,phi;
    double thetastep=M_PI/Ntheta, phistep=2*M_PI/Nphi;
    double volume=a*b*M_PI*L;

    for(ntheta=0; ntheta < Ntheta; ntheta++){
        theta=(ntheta+.5)*thetastep;
        sinc=Sinc( L/2.*Q*cos(theta) );
        for(nphi=0; nphi<Nphi; nphi++){
            phi=(nphi+.5)*phistep;
            if( (r0!=0.) || (theta0!=0) || (phi0!=0) )
                cosr0q=sin(theta0)*sin(theta)*cos(phi0-phi)+cos(theta0)*cos(theta);
            r=sin(theta)*sqrt( pow(a*cos(phi),2)+pow(b*sin(phi),2) );
            if(Q*r==0)
                F[ntheta][nphi]+=volume*rho*pol(1.0*sinc,0.0);
            else
                F[ntheta][nphi]+=volume*rho*pol(2*j1(Q*r)/(Q*r)*sinc,-r0*Q*cosr0q);
        }
    }
    return 0;
}


double  NanodiscPDBModel(complex double ** alpha, double Q,
        double rho_head, double rho_alkyl, double rho_methyl, double rho_belt, double rho_solvent,
        double radius_major, double radius_minor, double Dbelt,
        double Hlipid, double Hcore, double Hmethyl, double Hbelt, 
        double VerticalAxisOfEllipsoid, double ScaleFactorOfEndcaps)
{
    double complex ** F = ComplexArray(Ntheta,Nphi);

//Calculate formfactor of Nanodisc consisting of concentric eliptical cylinders
    DiscWithEndcaps(F, radius_major, radius_minor, Hlipid,  ScaleFactorOfEndcaps, VerticalAxisOfEllipsoid,  rho_head-rho_solvent,0.,0.,0.,Q);     //Lipid bilayer
    DiscWithEndcaps(F, radius_major, radius_minor, Hcore,   ScaleFactorOfEndcaps, VerticalAxisOfEllipsoid, rho_alkyl-rho_head   ,0.,0.,0.,Q);     //Lipid core
    Discflat(F, radius_major      , radius_minor, fabs(Hmethyl), rho_methyl-rho_alkyl  , 0.,0.,0.,Q);                                             //Methyl layer
    Discflat(F, radius_major+Dbelt, radius_minor+Dbelt, Hbelt  ,  rho_belt-rho_solvent , 0.,0.,0.,Q);                                             //Belt
    Discflat(F, radius_major      , radius_minor,       Hbelt  ,-(rho_belt-rho_solvent), 0.,0.,0.,Q);                                             //Hole in belt
   
    //Expand formfactor amplitud in spherical harmonics with coefficients alpha_lm
    expand(alpha,F);

    /*
    // The Intensity can also be calculated by radial integration of F:

    double theta,phi;
    double thetastep=M_PI/(Ntheta), phistep=2*M_PI/(Nphi);
    int ntheta,nphi;
    double Int=0;
    for(ntheta=0;ntheta<Ntheta;ntheta++){ //integration over theta
        theta=thetastep*(ntheta+.5);
        for(nphi=0;nphi<Nphi;nphi++){  //integration over phi
            phi=phistep*(nphi+.5);
            Int += sin(theta)*pow(cabs(F[ntheta][nphi]),2)*phistep*thetastep;
        }
    }
    printf("%g  %g \n",Q,Int);
*/

    FreeComplexArray(F,Ntheta,Nphi);

    return 0;
}
