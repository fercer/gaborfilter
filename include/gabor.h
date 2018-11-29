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

#ifndef GABOR_DLL_H_INCLUDED
#define GABOR_DLL_H_INCLUDED

#ifdef COMPILE_PYTHON
#include <Python.h>
#include <numpy/ndarraytypes.h>
#include <numpy/ufuncobject.h>
#include <numpy/npy_3kcompat.h>
#endif


#ifdef COMPILE_ACML
#include <acml.h>
#define ft_complex doublecomplex
#define ft_real(POINTER) (POINTER)->real
#define ft_imag(POINTER) (POINTER)->imag
#define allocate_ft_complex(ARRAY_SIZE) (doublecomplex*)malloc(ARRAY_SIZE*sizeof(doublecomplex))
#define allocate_ft_complex_pointers(ARRAY_SIZE) (doublecomplex**)malloc(ARRAY_SIZE*sizeof(doublecomplex*))
#define allocate_ft_complex_pointers_of_pointers(ARRAY_SIZE) (doublecomplex***)malloc(ARRAY_SIZE*sizeof(doublecomplex**))
#define deallocate_ft_complex(COMPLEX_ARRAY) free(COMPLEX_ARRAY)
#define ft_variables(HEIGHT, WIDTH) double* communication_work_array = (double*)malloc((4*WIDTH + 6*HEIGHT + 300) * sizeof(double)); int information_integer
#define ft_forward_setup(HEIGHT, WIDTH, INPUT, OUTPUT) dzfft2d(0, HEIGHT, WIDTH, INPUT, OUTPUT, communication_work_array, &information_integer)
#define ft_backward_setup(HEIGHT, WIDTH, INPUT, OUTPUT) zdfft2d(0, HEIGHT, WIDTH, INPUT, OUTPUT, communication_work_array, &information_integer)
#define ft_forward(HEIGHT, WIDTH, INPUT, OUTPUT) dzfft2d(1, HEIGHT, WIDTH, INPUT, OUTPUT, communication_work_array, &information_integer)
#define ft_backward(HEIGHT, WIDTH, INPUT, OUTPUT) zdfft2d(1, HEIGHT, WIDTH, INPUT, OUTPUT, communication_work_array, &information_integer)
#define ft_release_forward 
#define ft_release_backward 
#define ft_close free(communication_work_array)
#else
#include <fftw3.h>
#define ft_complex fftw_complex
#define ft_real(POINTER) *(*(POINTER))
#define ft_imag(POINTER) *(*(POINTER) + 1)
#define allocate_ft_complex(ARRAY_SIZE) (fftw_complex*)fftw_malloc(ARRAY_SIZE*sizeof(fftw_complex))
#define allocate_ft_complex_pointers(ARRAY_SIZE) (fftw_complex**)fftw_malloc(ARRAY_SIZE*sizeof(fftw_complex*))
#define allocate_ft_complex_pointers_of_pointers(ARRAY_SIZE) (fftw_complex***)fftw_malloc(ARRAY_SIZE*sizeof(fftw_complex**))
#define deallocate_ft_complex(COMPLEX_ARRAY) fftw_free(COMPLEX_ARRAY)
#define ft_variables(HEIGHT, WIDTH) fftw_plan forward_plan, backward_plan; char forward_plan_active = 0, backward_plan_active = 0
#define ft_forward_setup(HEIGHT, WIDTH, INPUT, OUTPUT) forward_plan = fftw_plan_dft_r2c_2d(HEIGHT, WIDTH, INPUT, OUTPUT,  FFTW_ESTIMATE); forward_plan_active = 1
#define ft_backward_setup(HEIGHT, WIDTH, INPUT, OUTPUT) backward_plan = fftw_plan_dft_c2r_2d(HEIGHT, WIDTH, INPUT, OUTPUT, FFTW_ESTIMATE); backward_plan_active = 1
#define ft_forward(HEIGHT, WIDTH, INPUT, OUTPUT) fftw_execute(forward_plan)
#define ft_backward(HEIGHT, WIDTH, INPUT, OUTPUT) fftw_execute(backward_plan)
#define ft_release_forward fftw_destroy_plan(forward_plan); forward_plan_active = 0
#define ft_release_backward fftw_destroy_plan(backward_plan); backward_plan_active = 0
#define ft_close if (forward_plan_active){fftw_destroy_plan(forward_plan);} if (backward_plan_active){fftw_destroy_plan(backward_plan);}
#endif

#ifndef NDEBUG
#define DEBMSG(MESSAGE) printf(MESSAGE)
#define DEBNUMMSG(MESSAGE, NUM) printf(MESSAGE, NUM);
#else
#define DEBMSG(MESSAGE) 
#define DEBNUMMSG(MESSAGE, NUM) 
#endif

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#define MY_PI 3.1415926535897932384626433832795028841971693993751058209749445923078164

#if defined(_WIN32) || defined(_WIN64)
#ifdef BUILDING_GABOR_DLL
#define GABOR_DLL __declspec(dllexport)
#else
#define GABOR_DLL __declspec(dllimport)
#endif
#else
#define GABOR_DLL 
#endif

void generateGaborKernels(ft_complex** gabor_kernels, const unsigned int height, const unsigned int width, const unsigned int par_T, const double par_L, const unsigned int par_K);
void generateHPF(ft_complex* high_pass_filter, const unsigned int height, const unsigned int width, const unsigned int par_T, const double par_L);

void gaborFilter(ft_complex* input, double* output, const unsigned int height, const unsigned int width, const unsigned int par_K, ft_complex** gabor_kernels, ft_complex* high_pass_filter);

void GABOR_DLL singleScaleGaborFilter(double * raw_input, char * mask, double * output, const unsigned int height, const unsigned width, const unsigned int par_T, const double par_L, const unsigned int par_K);
void GABOR_DLL singleScaleGaborFilter_multipleinputs(double * raw_input, const unsigned int n_inputs, char * mask, double * output, const unsigned int height, const unsigned width, const unsigned int par_T, const double par_L, const unsigned int par_K);

void GABOR_DLL multiscaleGaborFilter(double * raw_input, char * mask, double * output, const unsigned int height, const unsigned width, unsigned int * par_T, const unsigned int t_scales, const double par_L, const unsigned int par_K);
void GABOR_DLL multiscaleGaborFilter_multipleinputs(double * raw_input, const unsigned int n_inputs, char * mask, double * output, const unsigned int height, const unsigned width, unsigned int * par_T, const unsigned int t_scales, const double par_L, const unsigned int par_K);

//static PyObject* gaborFilter(PyObject *self, PyObject *args);


#endif // GABOR_DLL_H_INCLUDED
