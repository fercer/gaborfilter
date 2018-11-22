/***********************************************************************************************************************************
CENTRO DE INVESTIGACION EN MATEMATICAS
DOCTORADO EN CIENCIAS DE LA COMPUTACION
FERNANDO CERVANTES SANCHEZ

FILE NAME : gabor.h

PURPOSE : Declares the functions required to filter images using the Gabor filter.

FILE REFERENCES :
Name        I / O        Description
None----       ----------

ABNORMAL TERMINATION CONDITIONS, ERROR AND WARNING MESSAGES :
None
************************************************************************************************************************************/

//#ifdef COMPILE_PYTHON
#include <Python.h>
#include <numpy/ndarraytypes.h>
#include <numpy/ufuncobject.h>
#include <numpy/npy_3kcompat.h>
//#endif


//#ifdef COMPILE_ACML
#include <acml.h>
//#endif


#include <stdio.h>
#include <math.h>
#include <limits>
#include <vector>

#define MY_PI 3.1415926535897932384626433832795028841971693993751058209749445923078164

doublecomplex ** generateGaborKernels(const unsigned int height, const unsigned int width, const unsigned int par_T, const double par_L, const unsigned int par_K);

PyObject* singleScaleGaborFilter_impl(PyArrayObject *input, const unsigned int par_T, const double par_L, const unsigned int par_K);

PyObject* multiScaleGaborFilter_impl(PyArrayObject *input, const unsigned int * par_T, const unsigned int t_scales, const double par_L, const unsigned int par_K);

static PyObject* gaborFilter(PyObject *self, PyObject *args);