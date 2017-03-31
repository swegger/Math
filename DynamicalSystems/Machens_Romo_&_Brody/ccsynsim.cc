#include "mex.h"

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




void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[])
{
  if(nrhs!=1) 
    mexErrMsgTxt("One input required.");
  if(nlhs>0) 
    mexErrMsgTxt("No outputs, please : params passed by var, not value.");
  if (!mxIsStruct(prhs[0]) )
    mexErrMsgTxt("Input must be a scalar structure.");
  if ( mxGetN(prhs[0]) != 1  || mxGetM(prhs[0]) != 1 )
    mexErrMsgTxt("No structure arrays-- must be single scalar structure.");
  int A,B;

  retrieveScalarParam(dt);
  retrieveScalarInteger(Nt);
  retrieveDoubleVectorPtr(t, 1, Nt);

  retrieveScalarParam(satmax);
  retrieveScalarParam(tau);

  retrieveScalarInteger(nspikes);
  retrieveDoubleVectorPtr( s, 1, Nt );
  retrieveDoubleVectorPtr( spikes, 1, nspikes );

  int i=0;
  for ( int k=0; k<Nt-1; k++ ) {
    if ( spikes[i]>t[k]-dt/2 && spikes[i]<=t[k]+dt/2 ) {
      s[k+1] = s[k] + (satmax-s[k])/satmax;
      i++;
    }
    else {
      double dsdt = -s[k];
      s[k+1] = s[k] + dsdt*dt/tau;
    }
  }
}
