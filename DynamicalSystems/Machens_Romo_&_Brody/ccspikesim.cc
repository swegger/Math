#include "mex.h"
#include <math.h>
#include <stdlib.h>

/* (c) 2004 CK Machens & CD Brody */

#define retrieveScalarParam(T) \
  mxArray * ptr_##T = mxGetField(prhs[0], 0, #T); \
  if ( ptr_##T==NULL ) mexErrMsgTxt("Structure must have field '" #T "'."); \
  if ( !mxIsNumeric(ptr_##T) ) mexErrMsgTxt("'" #T "' must be numeric."); \
  if ( mxGetM(ptr_##T)!=1 || mxGetN(ptr_##T) != 1 ) \
            mexErrMsgTxt("'" #T "' must be scalar."); \
  double T = mxGetScalar(ptr_##T); 


#define retrieveScalarInteger(T) \
  mxArray * ptr_##T = mxGetField(prhs[0], 0, #T); \
  if ( ptr_##T==NULL ) mexErrMsgTxt("Structure must have field '" #T "'."); \
  if ( !mxIsNumeric(ptr_##T) ) mexErrMsgTxt("'" #T "' must be numeric."); \
  if ( mxGetM(ptr_##T)!=1 || mxGetN(ptr_##T) != 1 ) \
            mexErrMsgTxt("'" #T "' must be scalar."); \
  int T = (int)mxGetScalar(ptr_##T); 


#define retrieveDoubleVectorPtr(T,M,N)		  \
  mxArray * ptr_##T = mxGetField(prhs[0], 0, #T); \
  if ( ptr_##T==NULL ) mexErrMsgTxt("Structure must have field '" #T "'."); \
  if ( !mxIsNumeric(ptr_##T) ) mexErrMsgTxt("'" #T "' must be numeric."); \
  A = mxGetM( ptr_##T ); \
  B = mxGetN( ptr_##T ); \
  if ( !((A==M)&&(B==N)) )			   \
    if ( !((A==N)&&(B==M)&&( (M==1) || (N==1) )) ) \
      mexErrMsgTxt( "'" #T "' has wrong size." ); \
  double * T = mxGetPr(ptr_##T); \


#define retrieveCell(T) \
  mxArray * T = mxGetField(prhs[0], 0, #T); \
  if ( T == NULL ) mexErrMsgTxt("Structure must have field '" #T "'."); \
  if ( !mxIsCell( T ) ) mexErrMsgTxt(#T " must be a cell.");


#define makeVector(X, N) \
  mxArray * mx##X = mxCreateDoubleMatrix(1, N, mxREAL); \
  double * X = mxGetPr(mx##X);


#define makeIntVector(X, N) \
  int * X = (int *)mxCalloc((int) N, sizeof(int)); \
  int n##X = 0;

  
double gaussrand() {
  static int iset=0;
  static double gset;

  if ( iset==0 ) {
    double rsqu, v1, v2;
    do { 
      v1 = 2*drand48()-1;
      v2 = 2*drand48()-1;
      rsqu = v1*v1 + v2*v2;
    } while ( rsqu >= 1 || rsqu==0 );
    
    double fac = sqrt(-2*log(rsqu)/rsqu);
    gset = v1*fac; iset=1;
    return v2*fac;
  } else {
    iset = 0;
    return gset;
  }
}

void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[])
{
  // check for proper number of arguments
  if(nrhs!=1) 
    mexErrMsgTxt("One input required.");
  if(nlhs>0) 
    mexErrMsgTxt("No outputs, please : params passed by var, not value.");
  if (!mxIsStruct(prhs[0]) )
    mexErrMsgTxt("Input must be a scalar structure.");
  if ( mxGetN(prhs[0]) != 1  || mxGetM(prhs[0]) != 1 )
    mexErrMsgTxt("No structure arrays-- must be single scalar structure.");
  
  int A,B;

  // simulation parameters
  retrieveScalarParam(T);
  retrieveScalarInteger(Nt);
  retrieveScalarParam(dt);
  retrieveScalarInteger(report_every);
  retrieveScalarInteger(Nrt);
  retrieveDoubleVectorPtr( t, 1, Nt );

  // iaf parameters
  retrieveScalarParam(C);
  retrieveScalarParam(gleak);
  retrieveScalarParam(Vleak);
  retrieveScalarParam(Vthresh);
  retrieveScalarParam(Vreset);
  retrieveScalarParam(refrac);
  retrieveScalarParam(EE);
  retrieveScalarParam(EI);
  retrieveScalarParam(gaussnoise);

  // synaptic parameters
  retrieveScalarParam(Isyn_satmax);
  retrieveScalarParam(Itausyn);
  retrieveScalarParam(Igsyn);
  
  // connectivity
  retrieveScalarInteger(Nneurons);
  retrieveCell(Iweight_ids);
  retrieveCell(Iweight_vals);

  // external inputs
  retrieveDoubleVectorPtr( gE1, 1, Nt );
  retrieveDoubleVectorPtr( gE2, 1, Nt );

  // inner variables
  retrieveDoubleVectorPtr( vv,      Nneurons, 1 );
  retrieveDoubleVectorPtr( ssoutI,  Nneurons, 1 );
  retrieveDoubleVectorPtr( ssinI,   Nneurons, 1 );
  retrieveDoubleVectorPtr( last_spike_time, Nneurons, 1 );

  //reporting variables
  retrieveDoubleVectorPtr( v,       Nneurons, Nrt );
  retrieveDoubleVectorPtr( soutI,   Nneurons, Nrt );
  retrieveDoubleVectorPtr( sinI,    Nneurons, Nrt );
  retrieveDoubleVectorPtr( spikes,  Nneurons, Nrt );
  retrieveDoubleVectorPtr( nspikes, Nneurons, 1 );

  // auxiliary variables
  makeIntVector(refracted,     Nneurons);
  makeIntVector(not_refracted, Nneurons);
  makeIntVector(spikeguys,     Nneurons);

  // simulation
  double dvdt; double thisv;
  int report_k = 1;
  for ( int k=1; k<Nt; k++ ) {
    nrefracted=0; nnot_refracted=0;
    for ( int n=0; n<Nneurons; n++ ) 
      if ( t[k] - last_spike_time[n] < refrac )
	{ refracted[nrefracted++] = n; vv[n] = Vreset; }
      else
	not_refracted[nnot_refracted++] = n;      
    
    for( int n=0; n<nnot_refracted; n++ ) {
      thisv = vv[not_refracted[n]];
      
      if ( not_refracted[n]<Nneurons/2 )
	dvdt = (EE-thisv) * gE1[k];
      else
	dvdt = (EE-thisv) * gE2[k];
      dvdt += (Vleak - thisv)*gleak 
	+ (EI - thisv)*Igsyn*ssinI[not_refracted[n]]
	+ gaussnoise*gaussrand( );
      vv[not_refracted[n]] += dt*dvdt/C;
    }

    // Now, who spiked?
    nspikeguys=0;
    for ( int n=0; n<Nneurons; n++ ) 
      if ( vv[n] >= Vthresh ) spikeguys[nspikeguys++] = n;

    for ( int n=0; n<nspikeguys; n++ ) {
      mxArray *postguys;    mxArray *postweights;
      double  *mypostguys;  double  *mypostweights;
      int npost;

      vv[spikeguys[n]] = Vreset;

      // First calculate what the synaptic gating var should be
      ssoutI[spikeguys[n]] = 
	(Isyn_satmax + ssoutI[spikeguys[n]]*(1-Isyn_satmax)) *  
	exp(-(t[k] - last_spike_time[spikeguys[n]])/Itausyn);
      // Then store in the output just what should be added nex to the input
      ssoutI[spikeguys[n]] = 
	(Isyn_satmax - ssoutI[spikeguys[n]])/Isyn_satmax; 

      // Now, who do they connect to?
      postguys      = mxGetCell(Iweight_ids, spikeguys[n]);
      if ( postguys != NULL ) {
	npost = mxGetNumberOfElements(postguys);
	mypostguys    = mxGetPr(postguys);
	postweights   = mxGetCell(Iweight_vals, spikeguys[n]);
	mypostweights = mxGetPr(postweights);
	for( int p=0; p<npost; p++ )
	  ssinI[(int)mypostguys[p]] += mypostweights[p] * ssoutI[spikeguys[n]];
      }

      last_spike_time[spikeguys[n]] = t[k];
      nspikes[spikeguys[n]]++;
      spikes[((int)nspikes[spikeguys[n]]-1)*Nneurons + spikeguys[n]] = t[k];
    }

    // Now the synaptic decay
    for ( int n=0; n<Nneurons; n++ ) {
      ssinI[n] -= dt*ssinI[n]/Itausyn;
    }

    // Finally, reporting
    if ( (k-1)% ((int)report_every) == 0 ) {
      report_k++; int rn;
      for ( int n=0; n<Nneurons; n++ ) {
	rn = (report_k-1)*((int)Nneurons) + n;
	v[rn] = vv[n]; 
	sinI[rn] = ssinI[n]; soutI[rn] = ssoutI[n];
      }
    }
  }

}
