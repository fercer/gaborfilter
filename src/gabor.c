/***********************************************************************************************************************************
CENTRO DE INVESTIGACION EN MATEMATICAS
DOCTORADO EN CIENCIAS DE LA COMPUTACION
FERNANDO CERVANTES SANCHEZ

FILE NAME : gabor.cpp

PURPOSE : Defines the functions required to filter an image using the Gabor filter.

FILE REFERENCES :
Name        I / O        Description
None----       ----------

ABNORMAL TERMINATION CONDITIONS, ERROR AND WARNING MESSAGES :
None
************************************************************************************************************************************/
#include "gabor.h"

doublecomplex ** generateGaborKernels(const unsigned int height, const unsigned int width, const unsigned int par_T, const double par_L, const unsigned int par_K)
{
	doublecomplex ** gabor_kernels = (doublecomplex**)malloc(par_K * sizeof(doublecomplex*));

	const double sx = (double)par_T / (2.0 * sqrt(2.0 * log(2.0)));
	const double sy = par_L * sx;
	const double fx = 1.0 / (double)par_T;
	const double fy = 0.0;

	for (unsigned int k = 0; k < par_K; k++)
	{
		const double ctheta = cos((double)k * (double)par_K / 180.0 * MY_PI);
		const double stheta = sin((double)k * (double)par_K / 180.0 * MY_PI);

		for (unsigned int y = 0; y < height; y++)
		{
			for (unsigned int x = 0; x < width; x++)
			{
				const double u = (2.0 * (double)x - (double)(width*width)) * MY_PI;
				const double v = (2.0 * (double)y - (double)(height*height)) * MY_PI;

				const double rotated_u = u*ctheta + v*stheta;
				const double rotated_v = -u*stheta + v*ctheta;

				const double gabor_kernel_xy = exp(-0.5 * ((sx*sx*rotated_u*rotated_u + 4.0*MY_PI*MY_PI*fx*fx) + (sy*sy*rotated_v*rotated_v + 4.0*MY_PI*MY_PI*fy*fy))) * cosh(2.0*MY_PI*())
			}
		}
	}

	return gabor_kernels;
}



PyObject* singleScaleGaborFilter_impl(PyArrayObject *input, const unsigned int par_T, const double par_L, const unsigned int par_K)
{
	
}



PyObject* multiScaleGaborFilter_impl(PyArrayObject *input, const unsigned int * par_T, const unsigned int t_scales, const double par_L, const unsigned int par_K)
{

}

/*
npy_intp extracted_patches_dimensions[] = { sample_size, input_channels, patch_size, patch_size };
PyObject * extracted_patches = PyArray_SimpleNew(4, &extracted_patches_dimensions[0], NPY_DOUBLE);
char * extracted_patches_data = ((PyArrayObject*)extracted_patches)->data;
npy_intp extracted_patches_stride = ((PyArrayObject*)extracted_patches)->strides[((PyArrayObject*)extracted_patches)->nd - 1];

PyObject *patches_tuple = PyTuple_New(2);
PyTuple_SetItem(patches_tuple, 0, extracted_patches);
PyTuple_SetItem(patches_tuple, 1, max_values);
*/

PyObject *gaborFilter(PyObject *self, PyObject *args)
{
	PyArrayObject *input;
	unsigned int patch_size;
	unsigned int stride;
	double positive_proportion;
	PyArrayObject *sampled_indices;

	/* Parse the input arguments to extract two numpy arrays: */
	if (!PyArg_ParseTuple(args, "O!IIdO!", &PyArray_Type, &input, &patch_size, &stride, &positive_proportion, &PyArray_Type, &sampled_indices))
	{
		return NULL;
	}

	return extractbalancedpatchesMax_impl(input, positive_proportion, patch_size, stride, sampled_indices);
}


static PyMethodDef gabor_methods[] = {
	{ "gaborFilter",	gaborFilter, METH_VARARGS, "applies the Gabor filter to the input image, using the parameters t, l and K passed, if the parameter t is a list, then the multiscale Gabor filter is pplied instead." },
	{ NULL, NULL, 0, NULL }
};


static struct PyModuleDef gabor_moduledef = {
	PyModuleDef_HEAD_INIT,
	"gabor",
	NULL,
	-1,
	gabor_methods,
	NULL,
	NULL,
	NULL,
	NULL
};


PyMODINIT_FUNC PyInit_gabor(void)
{
	PyObject *m;
	m = PyModule_Create(&gabor_moduledef);
	if (!m) {
		return NULL;
	}
	import_array();

	return m;
}

/*^^^^^^^^^^^^^^^^^^^^^^^^^^^^ Python interface: ^^^^^^^^^^^^^^^^^^^^^^^^^^^^*/