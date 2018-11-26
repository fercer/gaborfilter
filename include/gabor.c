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

void generateGaborKernels(ft_complex** gabor_kernels, const unsigned int height, const unsigned int width, const unsigned int par_T, const double par_L, const unsigned int par_K)
{
	const double sx_2 = ((double)par_T / (2.0 * sqrt(2.0 * log(2.0)))) * ((double)par_T / (2.0 * sqrt(2.0 * log(2.0))));
	const double sy_2 = par_L * par_L * sx_2;
	const double fx = 1.0 / (double)par_T;
	const double fx_2 = fx * fx;
	const double fy = 0.0;
	const double fy_2 = fy * fy;
	const double pi_2 = MY_PI * MY_PI;

	printf("Sx^2 = %f, Sy^2 = %f\n", sx_2, sy_2);

	for (unsigned int k = 0; k < par_K; k++)
	{		
		const double ctheta = cos((double)k * (double)par_K / 180.0 * MY_PI);
		const double stheta = sin((double)k * (double)par_K / 180.0 * MY_PI);

		for (unsigned int y = 0; y < height; y++)
		{
			const double v = 2.0*MY_PI*(double)y/(double)height - MY_PI;
			for (unsigned int x = 0; x < (width/2 + 1); x++)
			{
				const double u = 2.0*MY_PI*(double)x/(double)width - MY_PI;

				const double rotated_u = u*ctheta + v*stheta;
				const double rotated_v = -u*stheta + v*ctheta;

				// Use only the real part of the kernel:
				ft_real(*(gabor_kernels + k) + y*(width/2 + 1) + x) = exp(-0.5 * ((sx_2*rotated_u*rotated_u + 4.0*pi_2*fx_2) + (sy_2*rotated_v*rotated_v + 4.0*pi_2*fy_2))) * cosh(2.0*MY_PI*(sx_2*fx*rotated_u + sy_2*fy*rotated_v));
				ft_imag(*(gabor_kernels + k) + y*(width/2 + 1) + x) = 0.0;
			}
		}
	}
	
	DEBMSG("Gabor filters generated successfully...\n");
}



void generateHPF(ft_complex* high_pass_filter, const unsigned int height, const unsigned int width, const unsigned int par_T, const double par_L)
{	
	const double sx_2 = ((double)par_T / (2.0 * sqrt(2.0 * log(2.0)))) * ((double)par_T / (2.0 * sqrt(2.0 * log(2.0))));
	const double sy_2 = par_L * par_L * sx_2;
	
	for (unsigned int y = 0; y < height; y++)
	{
		const double v = 2.0*MY_PI*(double)y/(double)height - MY_PI;
		for (unsigned int x = 0; x < (width/2 + 1); x++)
		{
			const double u = 2.0*MY_PI*(double)x/(double)width - MY_PI;
			ft_real(high_pass_filter + y*(width/2 + 1) + x) = 1.0 - exp(-(1/2)*(sy_2*(u*u + v*v)));
			ft_imag(high_pass_filter + y*(width/2 + 1) + x) = 0.0;
		}
	}
	
	DEBMSG("High pass filter generated successfully...\n");
}



void gaborFilter(ft_complex* input, double* output, const unsigned int height, const unsigned int width, const unsigned int par_K, doublecomplex** gabor_kernels, doublecomplex* high_pass_filter)
{
	ft_variables(height, width);
	
	double ** gabor_responses_to_kernels = (double**) malloc(par_K * sizeof(double*));
	ft_complex* fourier_response_to_kernel = allocate_ft_complex(height*(width/2+1));
	
	ft_complex* fr_ptr;
	for (unsigned int k = 0; k < par_K; k++)
	{	
		ft_complex* hpf_ptr = high_pass_filter;
		ft_complex* i_ptr = input;
		fr_ptr = fourier_response_to_kernel;
		ft_complex* gk_ptr = *(gabor_kernels + k);
		
		
		for (unsigned int xy = 0; xy < height*(width/2+1); xy++, i_ptr++, hpf_ptr++, fr_ptr++, gk_ptr++)
		{
			const double hpf_i_r = ft_real(hpf_ptr) * ft_real(i_ptr) - ft_imag(hpf_ptr) * ft_imag(i_ptr);
			const double hpf_i_i = ft_real(hpf_ptr) * ft_imag(i_ptr) + ft_imag(hpf_ptr) * ft_real(i_ptr);
			
			ft_real(fr_ptr) = hpf_i_r*ft_real(gk_ptr) - hpf_i_i*ft_imag(gk_ptr);
			ft_imag(fr_ptr) = hpf_i_r*ft_imag(gk_ptr) + hpf_i_i*ft_real(gk_ptr);
		}
		
		*(gabor_responses_to_kernels + k) = (double*) malloc(height * width * sizeof(double));

		ft_backward_setup(height, width, fourier_response_to_kernel, *(gabor_responses_to_kernels + k));
		ft_backward(height, width, fourier_response_to_kernel, *(gabor_responses_to_kernels + k));
		
		ft_release_backward;
	}
	ft_close;
	deallocate_ft_complex(fourier_response_to_kernel);

	DEBMSG("Gabor kernels applied successfully...\n");
	
	memcpy(output, *gabor_responses_to_kernels, height*width*sizeof(double));
	free(*gabor_responses_to_kernels);
	
	for (unsigned int k = 1; k < par_K; k++)
	{
		double* o_ptr = output;
		double* gr_ptr = *(gabor_responses_to_kernels + k);
		
		for (unsigned int xy = 0; xy < height*width; xy++, o_ptr++, gr_ptr++)
		{
			if (*o_ptr > *gr_ptr)
			{
				*o_ptr = *gr_ptr;
			}
		}
		free(*(gabor_responses_to_kernels + k));
	}	
	free(gabor_responses_to_kernels);
	DEBMSG("Gabor filter applied successfully...\n");
}



void singleScaleGaborFilter(double * raw_input, char * mask, double * output, const unsigned int height, const unsigned width, const unsigned int par_T, const double par_L, const unsigned int par_K)
{
	ft_variables(height, width);

	ft_complex* input = allocate_ft_complex(height * (width/2 + 1));
	/* Perform the Fourier Transform: */
	ft_forward_setup(height, width, raw_input, input);
	DEBMSG("Fourier transform setup ...\n");
	ft_forward(height, width, raw_input, input);
	DEBMSG("Raw input transformed with Fourier successfully...\n");
	ft_release_forward;
	
	
	// Compute the filter kernels:
	ft_complex** gabor_kernels = allocate_ft_complex_pointers(par_K);
	for (unsigned int k = 0; k < par_K; k++)
	{
		*(gabor_kernels+k) = allocate_ft_complex(height * (width/2+1));
	}	
	ft_complex* high_pass_filter = allocate_ft_complex(height*(width/2+1));
	
	generateGaborKernels(gabor_kernels, height, width, par_T, par_L, par_K);
	generateHPF(high_pass_filter, height, width, par_T, par_L);

	// Apply the single-scale filter:
	gaborFilter(input, output, height, width, par_K, gabor_kernels, high_pass_filter);

	// Apply the mask, the mask(x, y) must be {0, 1}:
	/*
	char * m_ptr = mask;
	double *o_ptr = output;
	for (unsigned int xy = 0; xy < height*width; xy++, m_ptr++, o_ptr++)
	{
		*o_ptr = *o_ptr * (double)*m_ptr;
	}
	*/
	
	// Free all memory used for image the filtering:
	for (unsigned int k = 0; k < par_K; k++)
	{
		deallocate_ft_complex(*(gabor_kernels + k));
	}
	deallocate_ft_complex(gabor_kernels);
	deallocate_ft_complex(high_pass_filter);
	deallocate_ft_complex(input);
	ft_close;
}



void singleScaleGaborFilter_multipleinputs(double * raw_input, const unsigned int n_inputs, char * mask, double * output, const unsigned int height, const unsigned width, const unsigned int par_T, const double par_L, const unsigned int par_K)
{
	ft_variables(height, width);
	
	ft_complex* input = allocate_ft_complex(height * (width/2 + 1));
	
	// Compute the filter kernels:
	ft_complex** gabor_kernels = allocate_ft_complex_pointers(par_K);
	for (unsigned int k = 0; k < par_K; k++)
	{
		*(gabor_kernels+k) = allocate_ft_complex(height * (width/2+1));
	}	
	ft_complex* high_pass_filter = allocate_ft_complex(height*(width/2+1));
	
	generateGaborKernels(gabor_kernels, height, width, par_T, par_L, par_K);
	generateHPF(high_pass_filter, height, width, par_T, par_L);

	for (unsigned int i = 0; i < n_inputs; i++)
	{			
		/* Perform the Fourier Transform: */
		ft_forward_setup(height, width, raw_input+i*height*width, input);
		ft_forward(height, width, raw_input+i*height*width, input);
		ft_release_forward;
		
		// Apply the single-scale filter:
		gaborFilter(input, output + height*width*i, height, width, par_K, gabor_kernels, high_pass_filter);

		// Apply the mask, the mask(x, y) must be {0, 1}:
		char * m_ptr = mask;
		double *o_ptr = output + height*width*i;
		for (unsigned int xy = 0; xy < height*width; xy++, m_ptr++, o_ptr++)
		{
			*o_ptr = *o_ptr * (double)*m_ptr;
		}
	}

	// Free all memory used for image the filtering:
	for (unsigned int k = 0; k < par_K; k++)
	{
		deallocate_ft_complex(*(gabor_kernels + k));
	}
	deallocate_ft_complex(gabor_kernels);
	deallocate_ft_complex(high_pass_filter);
	deallocate_ft_complex(input);
	ft_close;
}



void multiscaleGaborFilter(double * raw_input, char * mask, double * output, const unsigned int height, const unsigned width, unsigned int * par_T, const unsigned int t_scales, const double par_L, const unsigned int par_K)
{
	ft_variables(height, width);
	
	ft_complex* input = allocate_ft_complex(height * (width/2 + 1));
	
	/* Perform the Fourier Transform: */
	ft_forward_setup(height, width, raw_input, input);
	ft_forward(height, width, raw_input, input);
	ft_release_forward;
	
	// Compute the filter kernels:
	ft_complex** gabor_kernels = allocate_ft_complex_pointers(par_K);
	for (unsigned int k = 0; k < par_K; k++)
	{
		*(gabor_kernels+k) = allocate_ft_complex(height * (width/2+1));
	}	
	ft_complex* high_pass_filter = allocate_ft_complex(height*(width/2+1));
	
	for (unsigned int t = 0; t < t_scales; t++)
	{
		generateGaborKernels(gabor_kernels, height, width, *(par_T+t), par_L, par_K);
		generateHPF(high_pass_filter, height, width, *(par_T+t), par_L);

		// Apply the single-scale filter:
		gaborFilter(input, output + height*width*t, height, width, par_K, gabor_kernels, high_pass_filter);

		// Apply the mask, the mask(x, y) must be {0, 1}:
		char * m_ptr = mask;
		double *o_ptr = output + height*width*t;
		for (unsigned int xy = 0; xy < height*width; xy++, m_ptr++, o_ptr++)
		{
			*o_ptr = *o_ptr * (double)*m_ptr;
		}
	}
	
	// Free all memory used for image the filtering:
	for (unsigned int k = 0; k < par_K; k++)
	{
		deallocate_ft_complex(*(gabor_kernels + k));
	}
	deallocate_ft_complex(gabor_kernels);
	deallocate_ft_complex(high_pass_filter);
	deallocate_ft_complex(input);
	ft_close;
}




void multiscaleGaborFilter_multipleinputs(double * raw_input, const unsigned int n_inputs, char * mask, double * output, const unsigned int height, const unsigned width, unsigned int * par_T, const unsigned int t_scales, const double par_L, const unsigned int par_K)
{
	ft_variables(height, width);
	
	ft_complex* input = allocate_ft_complex(height * (width/2 + 1));
	
	// Compute the filter kernels:
	ft_complex*** gabor_kernels = allocate_ft_complex_pointers_of_pointers(t_scales);
	ft_complex** high_pass_filter = allocate_ft_complex_pointers(t_scales);
	for (unsigned int t = 0; t < t_scales; t++)
	{
		*(gabor_kernels+t) = allocate_ft_complex_pointers(par_K);		
		for (unsigned int k = 0; k < par_K; k++)
		{
			*(*(gabor_kernels+t)+k) = allocate_ft_complex(height * (width/2+1));
		}
		*(high_pass_filter+t) = allocate_ft_complex(height*(width/2+1));
		
		generateGaborKernels(*(gabor_kernels+t), height, width, *(par_T+t), par_L, par_K);
		generateHPF(*(high_pass_filter+t), height, width, *(par_T+t), par_L);
	}

	for (unsigned int i = 0; i < n_inputs; i++)
	{
		/* Perform the Fourier Transform: */
		ft_forward_setup(height, width, raw_input+height*width*i, input);
		ft_forward(height, width, raw_input+height*width*i, input);
		ft_release_forward;
		
		for (unsigned int t = 0; t < t_scales; t++)
		{
			// Apply the single-scale filter:
			gaborFilter(input, output + height*width*t, height, width, par_K, *(gabor_kernels+t), *(high_pass_filter+t));

			// Apply the mask, the mask(x, y) must be {0, 1}:
			char * m_ptr = mask + height*width*i;
			double *o_ptr = output + height*width*t + height*width*t_scales*i;
			for (unsigned int xy = 0; xy < height*width; xy++, m_ptr++, o_ptr++)
			{
				*o_ptr = *o_ptr * (double)*m_ptr;
			}
		}
	}
	
	// Free all memory used for image the filtering:
	for (unsigned int t = 0; t < t_scales; t++)
	{
		for (unsigned int k = 0; k < par_K; k++)
		{
			deallocate_ft_complex(*(*(gabor_kernels + t) + k));
		}
		deallocate_ft_complex(*(gabor_kernels+t));
		deallocate_ft_complex(*(high_pass_filter+t));
	}
	deallocate_ft_complex(gabor_kernels);
	deallocate_ft_complex(high_pass_filter);
	deallocate_ft_complex(input);
	ft_close;
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
/*
PyObject *gaborFilter(PyObject *self, PyObject *args)
{
	PyArrayObject *input;
	unsigned int patch_size;
	unsigned int stride;
	double positive_proportion;
	PyArrayObject *sampled_indices;

	// Parse the input arguments to extract two numpy arrays:
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
*/
