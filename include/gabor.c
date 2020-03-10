/***********************************************************************************************************************************
CENTRO DE INVESTIGACION EN MATEMATICAS
DOCTORADO EN CIENCIAS DE LA COMPUTACION
FERNANDO CERVANTES SANCHEZ

FILE NAME : gabor.c

PURPOSE : Defines the functions required to filter an image using the Gabor filter.

FILE REFERENCES :
Name        I / O        Description
None----       ----------

ABNORMAL TERMINATION CONDITIONS, ERROR AND WARNING MESSAGES :
None
************************************************************************************************************************************/
#include "gabor.h"

ft_complex* generateGaborKernel(const unsigned int nearest_2p_dim, const unsigned int par_T, const double par_L, const double theta)
{
	const double sx = (double)par_T / (2.0*sqrt(2.0*log(2.0)));
	const double sy = par_L * sx;
	const double sx_2 = sx*sx;
	const double sy_2 = sy*sy;
	const double fx = 1.0 / (double)par_T;

    /* Compute the Gabor filter in the spatial domain, then transform with Fourier to frequency domain */
    const double ctheta = cos(theta);
    const double stheta = sin(theta);

    double * gabor_spatial = (double*)calloc(nearest_2p_dim*nearest_2p_dim, sizeof(double));
    const double center_y = (double)nearest_2p_dim / 2.0;
    const double center_x = (double)nearest_2p_dim / 2.0;
    const double scale = 1.0/(2.0*MY_PI*sx*sy);

    double centered_x, centered_y, u, v;

    for (unsigned int y = 0; y < nearest_2p_dim; y++)
    {
        centered_y = (double)y - center_y;
        for (unsigned int x = 0; x < nearest_2p_dim; x++)
        {
            centered_x = (double)x - center_x;

            u = ctheta * centered_x - stheta * centered_y;
            v = ctheta * centered_y + stheta * centered_x;

            *(gabor_spatial + y*nearest_2p_dim + x) = scale * exp(-0.5*((u*u)/sx_2 + (v*v)/sy_2)) * cos(2.0*MY_PI*fx*u);
        }
    }

    /* Perform the Fourier Transform: */
    ft_variables(nearest_2p_dim, nearest_2p_dim);
    ft_complex * gabor_frequency = allocate_ft_complex(nearest_2p_dim * (nearest_2p_dim/2 + 1));

    ft_forward_setup(nearest_2p_dim, nearest_2p_dim, gabor_spatial, gabor_frequency);
    ft_forward(nearest_2p_dim, nearest_2p_dim, gabor_spatial, gabor_frequency);
    ft_release_forward;
    free(gabor_spatial);

    /* Remove the complex component by passing the information to the real part of the filter */
    for (unsigned int y = 0; y < nearest_2p_dim; y++)
    {
        centered_y = (double)y - center_y;
        for (unsigned int x = 0; x < nearest_2p_dim/2 + 1; x++)
        {
            ft_real(gabor_frequency + y*(nearest_2p_dim/2 + 1) + x) = sqrt(ft_real(gabor_frequency + y*(nearest_2p_dim/2 + 1) + x)*ft_real(gabor_frequency + y*(nearest_2p_dim/2 + 1) + x) + ft_imag(gabor_frequency + y*(nearest_2p_dim/2 + 1) + x)*ft_imag(gabor_frequency + y*(nearest_2p_dim/2 + 1) + x));
            ft_imag(gabor_frequency + y*(nearest_2p_dim/2 + 1) + x) = 0.0;
        }
    }

    return gabor_frequency;
}



ft_complex**  GABOR_DLL generateGaborBank(const unsigned int nearest_2p_dim, const unsigned int par_T, const double par_L, const unsigned int K)
{
    ft_complex ** gabor_bank = allocate_ft_complex_pointers(K);

    for (unsigned int k = 0; k < K; k++)
    {
        *(gabor_bank + k) = generateGaborKernel(nearest_2p_dim, par_T, par_L, (double)k/(double)K * MY_PI);
    }

    return gabor_bank;
}



ft_complex * generateHPF(const unsigned int nearest_2p_dim, const unsigned int par_T, const double par_L)
{	
	const double sx_2 = ((double)par_T / (2.0 * sqrt(2.0 * log(2.0)))) * ((double)par_T / (2.0 * sqrt(2.0 * log(2.0))));
	const double sy_2 = par_L * par_L * sx_2;

    ft_complex * high_pass_filter = allocate_ft_complex(nearest_2p_dim * (nearest_2p_dim/2 + 1));

    const double scale = 2.0*MY_PI/(double)nearest_2p_dim;
    const double center_x = (double)nearest_2p_dim/2.0;
    const double center_y = (double)nearest_2p_dim/2.0;
    double centered_x, centered_y;

	for (unsigned int y = 0; y < nearest_2p_dim/2; y++)
	{
        centered_y = ((double)y - center_y)*scale;
	    for (unsigned int x = 0; x < nearest_2p_dim/2+1; x++)
	    {
            centered_x = ((double)x - center_x) * scale;

	        ft_real(high_pass_filter + y*(nearest_2p_dim/2+1) + x) = 1.0 - exp(-0.5*sy_2*(centered_x**2 + centered_y**2));
	        ft_imag(high_pass_filter + y*(nearest_2p_dim/2+1) + x) = 0.0;

	        ft_real(high_pass_filter + (nearest_2p_dim - y - 1)*(nearest_2p_dim/2+1) + x) = 1.0 - exp(-0.5*sy_2*(centered_x**2 + centered_y**2));
	        ft_imag(high_pass_filter + (nearest_2p_dim - y - 1)*(nearest_2p_dim/2+1) + x) = 0.0;
	    }
	}

	return high_pass_filter;
}



ft_complex ** zeroPaddingImages(double * raw_input_data, const unsigned int n_imgs, const unsigned int height, const unsigned int width, const unsigned int nearest_2p_dim)
{
    ft_complex ** zero_padded_img_bank = allocate_ft_complex_pointers(n_imgs);
    double * zero_padded_img = (double*) calloc(nearest_2p_dim * nearest_2p_dim, sizeof(double));
    const unsigned int offset_x = (unsigned int)((double)(nearest_2p_dim - width)/2.0);
    const unsigned int offset_y = (unsigned int)((double)(nearest_2p_dim - height)/2.0);

    for (unsigned int i; i < n_imgs; i++)
    {
        for (unsigned int y = 0; y < height; y++)
        {
            memcpy(zero_padded_img + (y + offset_y)*nearest_2p_dim + offset_x, raw_input + i*height*width + y*width, width*sizeof(double));
        }

        /* Perform the Fourier Transform: */
        ft_variables(nearest_2p_dim, nearest_2p_dim);
        *(zero_padded_img_bank + i) =  allocate_ft_complex(nearest_2p_dim * (nearest_2p_dim/2 + 1));

        ft_forward_setup(nearest_2p_dim, nearest_2p_dim, zero_padded_img, *(zero_padded_img_bank + i));
        ft_forward(nearest_2p_dim, nearest_2p_dim, zero_padded_img, *(zero_padded_img_bank + i));
        ft_release_forward;
    }

    free(zero_padded_img);

	return zero_padded_img_bank;
}



void gaborFilter_impl(ft_complex* input, double* output, const unsigned int height, const unsigned int width, const unsigned int par_K, ft_complex** gabor_kernels, ft_complex* high_pass_filter)
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

	memcpy(output, *gabor_responses_to_kernels, height*width*sizeof(double));
	free(*gabor_responses_to_kernels);
	
	for (unsigned int k = 1; k < par_K; k++)
	{
		double* o_ptr = output;
		double* gr_ptr = *(gabor_responses_to_kernels + k);
		
		for (unsigned int xy = 0; xy < height*width; xy++, o_ptr++, gr_ptr++)
		{
			if (*o_ptr < *gr_ptr)
			{
				*o_ptr = *gr_ptr;
			}
		}
		free(*(gabor_responses_to_kernels + k));
	}	
	free(gabor_responses_to_kernels);
}



void multiscaleFilter(double * raw_input_data, double * gabor_response_data, const unsigned int n_imgs, const unsigned int height, const unsigned int width, unsigned int * par_tau_data, const unsigned int tau_scales, double * par_l_data, const unsigned int l_scales, par_K, unsigned char compute_max)
{

}


#ifdef BUILDING_PYTHON_MODULE
PyObject * gaborFilter(PyObject *self, PyObject *args)
{
    PyArrayObject *raw_input;
    char *raw_input_data = NULL;

    npy_intp n_imgs = 1;
    npy_intp height, width;

    PyArrayObject *multiscale_par_l = NULL;
    PyArrayObject *multiscale_par_tau = NULL;

    unsigned int par_K;
	unsigned char compute_max = 0;

	/** untrimmed kernel indicator and template src is an optional argument, which is specified after '|' */
    if (!PyArg_ParseTuple(args, "O!O!O!I|b", &PyArray_Type, &raw_input, &PyArray_Type, &multiscale_par_l, &PyArray_Type, &multiscale_par_tau, &par_K, &compute_max))
    {
        return NULL;
    }

    double * par_l_data = (double*)((PyArrayObject*)multiscale_par_l)->data;
    unsigned int * par_tau_data = (unsigned int*)((PyArrayObject*)multiscale_par_tau)->data;

    unsigned int l_scales = (PyArrayObject*)multiscale_par_l->size;
    unsigned int tau_scales = (PyArrayObject*)multiscale_par_tau->size;

    if (((PyArrayObject*)raw_input)->nd > 3)
    {
        n_imgs = ((PyArrayObject*)raw_input)->dimensions[0];
		/* Do not consider channels */
        height = ((PyArrayObject*)raw_input)->dimensions[2];
        width = ((PyArrayObject*)raw_input)->dimensions[3];
    }
    else if (((PyArrayObject*)raw_input)->nd > 2)
    {
        n_imgs = ((PyArrayObject*)raw_input)->dimensions[0];
        height = ((PyArrayObject*)raw_input)->dimensions[1];
        width = ((PyArrayObject*)raw_input)->dimensions[2];
    }
    else
    {
        n_imgs = 1;
        height = ((PyArrayObject*)raw_input)->dimensions[0];
        width = ((PyArrayObject*)raw_input)->dimensions[1];
    }

    raw_input_data  = (double*)((PyArrayObject*)raw_input)->data;

    npy_intp gabor_response_shape[] = { n_imgs, tau_scales * (compute_max ? 1 : par_K) * l_scales * (compute_max ? 1 : par_K), height, width };
    PyObject * gabor_response = PyArray_SimpleNew(4, &gabor_response_shape[0], NPY_DOUBLE);

    double * gabor_response_data = (double*)((PyArrayObject*)gabor_response)->data;
	multiscaleFilter(raw_input_data, gabor_response_data, n_imgs, height, width, par_tau_data, tau_scales, par_l_data, l_scales, par_K, compute_max);

    return gabor_response;
}


static PyMethodDef gabor_methods[] = {
	{ "gaborFilter", gaborFilter, METH_VARARGS, "applies the Gabor filter to the input image, using the parameters t, l and K passed, if the parameter t is a list, then the multiscale Gabor filter is pplied instead." },
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
#endif // BUILDING_PYTHON_MODULE