// Wrapping of WCSTools wcscon.c 
// This does general coord conversion and does not need a WCS

%feature("autodoc", "0");

%module (package="PyWCSTools") wcscon

%include "typemaps.i"

// %apply double *OUTPUT {double *dtheta, double *dphi, double *ptheta, double *pphi};
// void wcsconp (int sys1, int sys2, double eq1, double eq2, double ep1, double ep2, double *dtheta, double *dphi, double *ptheta, double *pphi);
// void wcsconp (int sys1, int sys2, double eq1, double eq2, double ep1, double ep2, double *dtheta, double *dphi, double *ptheta, double *pphi);


// %apply double *OUTPUT {double *dtheta, double *dphi, double *ptheta, double *pphi, double *px, double *rv};
// void wcsconv (int sys1, int sys2, double eq1, double eq2, double ep1, double ep2, double *dtheta, double *dphi, double *ptheta, double *pphi, double *px, double *rv);
// void wcsconv (int sys1, int sys2, double eq1, double eq2, double ep1, double ep2, double *dtheta, double *dphi, double *ptheta, double *pphi, double *px, double *rv);


// %apply double *OUTPUT {double *dtheta, double *dphi};
// void wcscon (int sys1, int sys2, double eq1, double eq2, double *dtheta, double *dphi, double epoch);
void wcscon (int sys1, int sys2, double eq1, double eq2, double *INOUT, double *INOUT, double epoch);

int wcscsys (char *wcstring);




