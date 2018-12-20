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
	const double sx_2 = (double)(par_T * par_T) / (8.0 * log(2.0));
	const double sy_2 = par_L * par_L * sx_2;
	const double fx = 1.0 / (double)par_T;
	const double fx_2 = fx * fx;
	const double fy = 0.0;
	const double fy_2 = fy * fy;
	const double pi_2 = MY_PI * MY_PI;

	for (unsigned int k = 0; k < par_K; k++)
	{	
		const double ctheta = cos((double)k / (double)par_K * MY_PI);
		const double stheta = sin((double)k / (double)par_K * MY_PI);
		
		for (unsigned int y = height/2; y < height; y++)
		{
			const double v = MY_PI*(2.0*(double)y/(double)height - 1.0);
			{
				const double u = -MY_PI;

				const double rotated_u = u*ctheta + v*stheta;
				const double rotated_v = -u*stheta + v*ctheta;

				// Use only the real part of the kernel:
				ft_real(*(gabor_kernels + k) + (y - height/2)*(width/2 + 1) + width/2) = exp(-(sx_2*(rotated_u*rotated_u + 4.0*pi_2*fx_2) + sy_2*(rotated_v*rotated_v + 4.0*pi_2*fy_2))/2.0) * cosh(2.0*MY_PI*(sx_2*fx*rotated_u + sy_2*fy*rotated_v));
				ft_imag(*(gabor_kernels + k) + (y - height/2)*(width/2 + 1) + width/2) = 0.0;
			}
			for (unsigned int x = width/2; x < width; x++)
			{
				const double u = MY_PI*(2.0*(double)x/(double)width - 1.0);

				const double rotated_u = u*ctheta + v*stheta;
				const double rotated_v = -u*stheta + v*ctheta;

				// Use only the real part of the kernel:
				ft_real(*(gabor_kernels + k) + (y - height/2)*(width/2 + 1) + x - width/2) = exp(-(sx_2*(rotated_u*rotated_u + 4.0*pi_2*fx_2) + sy_2*(rotated_v*rotated_v + 4.0*pi_2*fy_2))/2.0) * cosh(2.0*MY_PI*(sx_2*fx*rotated_u + sy_2*fy*rotated_v));
				ft_imag(*(gabor_kernels + k) + (y - height/2)*(width/2 + 1) + x - width/2) = 0.0;
			}
		}
		
		for (unsigned int y = 0; y < height/2; y++)
		{
			const double v = MY_PI*(2.0*(double)y/(double)height - 1.0);
			{
				const double u = -MY_PI;

				const double rotated_u = u*ctheta + v*stheta;
				const double rotated_v = -u*stheta + v*ctheta;

				// Use only the real part of the kernel:
				ft_real(*(gabor_kernels + k) + (y + height/2)*(width/2 + 1) + width/2) = exp(-(sx_2*(rotated_u*rotated_u + 4.0*pi_2*fx_2) + sy_2*(rotated_v*rotated_v + 4.0*pi_2*fy_2))/2.0) * cosh(2.0*MY_PI*(sx_2*fx*rotated_u + sy_2*fy*rotated_v));
				ft_imag(*(gabor_kernels + k) + (y + height/2)*(width/2 + 1) + width/2) = 0.0;
			}
			for (unsigned int x = width/2; x < width; x++)
			{
				const double u = MY_PI*(2.0*(double)x/(double)width - 1.0);

				const double rotated_u = u*ctheta + v*stheta;
				const double rotated_v = -u*stheta + v*ctheta;

				// Use only the real part of the kernel:
				ft_real(*(gabor_kernels + k) + (y + height/2)*(width/2 + 1) + x - width/2) = exp(-(sx_2*(rotated_u*rotated_u + 4.0*pi_2*fx_2) + sy_2*(rotated_v*rotated_v + 4.0*pi_2*fy_2))/2.0) * cosh(2.0*MY_PI*(sx_2*fx*rotated_u + sy_2*fy*rotated_v));
				ft_imag(*(gabor_kernels + k) + (y + height/2)*(width/2 + 1) + x - width/2) = 0.0;
			}
		}
	}	
}



void generateHPF(ft_complex* high_pass_filter, const unsigned int height, const unsigned int width, const unsigned int par_T, const double par_L)
{	
	const double sx_2 = ((double)par_T / (2.0 * sqrt(2.0 * log(2.0)))) * ((double)par_T / (2.0 * sqrt(2.0 * log(2.0))));
	const double sy_2 = par_L * par_L * sx_2;

	for (unsigned int y = height/2; y < height; y++)
	{
		const double v = 2.0*MY_PI*(double)y/(double)height - MY_PI;
		{
			const double u = -MY_PI;
			ft_real(high_pass_filter + (y - height/2)*(width/2 + 1) + width/2) = 1.0 - exp(-(sy_2*(u*u + v*v))/2.0);
			ft_imag(high_pass_filter + (y - height/2)*(width/2 + 1) + width/2) = 0.0;
		}
		for (unsigned int x = width/2; x < width; x++)
		{
			const double u = 2.0*MY_PI*(double)x/(double)width - MY_PI;
			ft_real(high_pass_filter + (y - height/2)*(width/2 + 1) + x - width/2) = 1.0 - exp(-(sy_2*(u*u + v*v))/2.0);
			ft_imag(high_pass_filter + (y - height/2)*(width/2 + 1) + x - width/2) = 0.0;
		}
	}
	
	for (unsigned int y = 0; y < height/2; y++)
	{
		const double v = 2.0*MY_PI*(double)y/(double)height - MY_PI;
		{
			const double u = -MY_PI;
			ft_real(high_pass_filter + (y + height/2)*(width/2 + 1) + width/2) = 1.0 - exp(-(sy_2*(u*u + v*v))/2.0);
			ft_imag(high_pass_filter + (y + height/2)*(width/2 + 1) + width/2) = 0.0;
		}
		for (unsigned int x = width/2; x < width; x++)
		{
			const double u = 2.0*MY_PI*(double)x/(double)width - MY_PI;
			ft_real(high_pass_filter + (y + height/2)*(width/2 + 1) + x - width/2) = 1.0 - exp(-(sy_2*(u*u + v*v))/2.0);
			ft_imag(high_pass_filter + (y + height/2)*(width/2 + 1) + x - width/2) = 0.0;
		}
	}
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
			const double hpf_i_r = ft_real(hpf_ptr) * ft_real(i_ptr) / (double)(height*width) - ft_imag(hpf_ptr) * ft_imag(i_ptr) / (double)(height*width);
			const double hpf_i_i = ft_real(hpf_ptr) * ft_imag(i_ptr) / (double)(height*width) + ft_imag(hpf_ptr) * ft_real(i_ptr) / (double)(height*width);
			
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


void gaborFilterWithAngles_impl(ft_complex* input, double* output, double * angles_output, const unsigned int height, const unsigned int width, const unsigned int par_K, ft_complex** gabor_kernels, ft_complex* high_pass_filter)
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
			const double hpf_i_r = ft_real(hpf_ptr) * ft_real(i_ptr) / (double)(height*width) - ft_imag(hpf_ptr) * ft_imag(i_ptr) / (double)(height*width);
			const double hpf_i_i = ft_real(hpf_ptr) * ft_imag(i_ptr) / (double)(height*width) + ft_imag(hpf_ptr) * ft_real(i_ptr) / (double)(height*width);
			
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
		double* a_ptr = angles_output;
		double* gr_ptr = *(gabor_responses_to_kernels + k);
		
		for (unsigned int xy = 0; xy < height*width; xy++, o_ptr++, a_ptr++, gr_ptr++)
		{
			if (*o_ptr < *gr_ptr)
			{
				*o_ptr = *gr_ptr;
				*a_ptr = (double)k / (double)par_K * 2.0 * MY_PI;
			}
		}
		free(*(gabor_responses_to_kernels + k));
	}	
	free(gabor_responses_to_kernels);
}


void singleScaleGaborFilter(double * raw_input, char * mask, double * output, const unsigned int height, const unsigned width, const unsigned int par_T, const double par_L, const unsigned int par_K)
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
	
	generateGaborKernels(gabor_kernels, height, width, par_T, par_L, par_K);
	generateHPF(high_pass_filter, height, width, par_T, par_L);

	// Apply the single-scale filter:
	gaborFilter_impl(input, output, height, width, par_K, gabor_kernels, high_pass_filter);

	// Apply the mask, the mask(x, y) must be {0, 1}:
	char * m_ptr = mask;
	double *o_ptr = output;
	for (unsigned int xy = 0; xy < height*width; xy++, m_ptr++, o_ptr++)
	{
		*o_ptr = *o_ptr * (double)*m_ptr;
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
		gaborFilter_impl(input, output + height*width*i, height, width, par_K, gabor_kernels, high_pass_filter);

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
		gaborFilter_impl(input, output + height*width*t, height, width, par_K, gabor_kernels, high_pass_filter);

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
			gaborFilter_impl(input, output + height*width*t + height*width*t_scales*i, height, width, par_K, *(gabor_kernels+t), *(high_pass_filter+t));

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


void singleScaleGaborFilterWithAngles(double * raw_input, char * mask, double * output, double * angles_output, const unsigned int height, const unsigned width, const unsigned int par_T, const double par_L, const unsigned int par_K)
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
	
	generateGaborKernels(gabor_kernels, height, width, par_T, par_L, par_K);
	generateHPF(high_pass_filter, height, width, par_T, par_L);

	// Apply the single-scale filter:
	gaborFilterWithAngles_impl(input, output, angles_output, height, width, par_K, gabor_kernels, high_pass_filter);

	// Apply the mask, the mask(x, y) must be {0, 1}:
	char * m_ptr = mask;
	double *o_ptr = output;
	double *a_ptr = angles_output;
	
	for (unsigned int xy = 0; xy < height*width; xy++, m_ptr++, o_ptr++, a_ptr++)
	{
		*o_ptr = *o_ptr * (double)*m_ptr;
		*a_ptr = *a_ptr * (double)*m_ptr;
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


void singleScaleGaborFilterWithAngles_multipleinputs(double * raw_input, const unsigned int n_inputs, char * mask, double * output, double * angles_output, const unsigned int height, const unsigned width, const unsigned int par_T, const double par_L, const unsigned int par_K)
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
		gaborFilterWithAngles_impl(input, output + height*width*i, angles_output + height*width*i, height, width, par_K, gabor_kernels, high_pass_filter);

		// Apply the mask, the mask(x, y) must be {0, 1}:
		char * m_ptr = mask;
		double *o_ptr = output + height*width*i;
		double *a_ptr = angles_output + height*width*i;
		for (unsigned int xy = 0; xy < height*width; xy++, m_ptr++, o_ptr++, a_ptr++)
		{
			*o_ptr = *o_ptr * (double)*m_ptr;
			*a_ptr = *a_ptr * (double)*m_ptr;
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


void multiscaleGaborFilterWithAngles(double * raw_input, char * mask, double * output, double * angles_output, const unsigned int height, const unsigned width, unsigned int * par_T, const unsigned int t_scales, const double par_L, const unsigned int par_K)
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
		gaborFilterWithAngles_impl(input, output + height*width*t, angles_output + height*width*t, height, width, par_K, gabor_kernels, high_pass_filter);

		// Apply the mask, the mask(x, y) must be {0, 1}:
		char * m_ptr = mask;
		double *o_ptr = output + height*width*t;
		double *a_ptr = angles_output + height*width*t;
		
		for (unsigned int xy = 0; xy < height*width; xy++, m_ptr++, o_ptr++, a_ptr++)
		{
			*o_ptr = *o_ptr * (double)*m_ptr;
			*a_ptr = *a_ptr * (double)*m_ptr;
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


void multiscaleGaborFilterWithAngles_multipleinputs(double * raw_input, const unsigned int n_inputs, char * mask, double * output, double * angles_output, const unsigned int height, const unsigned width, unsigned int * par_T, const unsigned int t_scales, const double par_L, const unsigned int par_K)
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
			gaborFilterWithAngles_impl(input, output + height*width*t + height*width*t_scales*i, angles_output + height*width*t + height*width*t_scales*i, height, width, par_K, *(gabor_kernels+t), *(high_pass_filter+t));

			// Apply the mask, the mask(x, y) must be {0, 1}:
			char * m_ptr = mask + height*width*i;
			double *o_ptr = output + height*width*t + height*width*t_scales*i;
			double *a_ptr = angles_output + height*width*t + height*width*t_scales*i;
			for (unsigned int xy = 0; xy < height*width; xy++, m_ptr++, o_ptr++, a_ptr++)
			{
				*o_ptr = *o_ptr * (double)*m_ptr;
				*a_ptr = *a_ptr * (double)*m_ptr;
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



#ifdef BUILDING_PYTHON_MODULE
PyObject *gaborFilter(PyObject *self, PyObject *args)
{
	PyArrayObject *raw_input;
	char *raw_input_data = NULL;
	npy_intp raw_input_stride, n_imgs = 1;
	npy_intp height, width;
	
	PyArrayObject *multiscale_par_T = NULL;
	char * par_T_data = NULL;
	npy_intp par_T_stride, par_T_scales = 1;
	
	PyArrayObject *mask = NULL;
	char * mask_data = NULL;
	npy_intp mask_stride;
	
	double par_L;
	unsigned int par_K;
	
	if (!PyArg_ParseTuple(args, "O!O!dIO!", &PyArray_Type, &raw_input, &PyArray_Type, &multiscale_par_T, &par_L, &par_K, &PyArray_Type, &mask))
	{
		return NULL;
	}
	
	par_T_data = ((PyArrayObject*)multiscale_par_T)->data;
	par_T_scales = ((PyArrayObject*)multiscale_par_T)->dimensions[0];
	par_T_stride = ((PyArrayObject*)multiscale_par_T)->strides[((PyArrayObject*)multiscale_par_T)->nd-1];

	if (((PyArrayObject*)raw_input)->nd > 2)
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
	
	DEBNUMMSG("N imgs: %i\n", (int)n_imgs);	
	DEBNUMMSG("T scales: %i\n", (int)par_T_scales);
	
	raw_input_data  = ((PyArrayObject*)raw_input)->data;
	raw_input_stride = ((PyArrayObject*)raw_input)->strides[((PyArrayObject*)raw_input)->nd - 1];
	
	mask_data = ((PyArrayObject*)mask)->data;
	mask_stride = ((PyArrayObject*)mask)->strides[((PyArrayObject*)mask)->nd - 1];

	PyObject * gabor_response = NULL;
	if (n_imgs > 1 && par_T_scales > 1)	
	{
		npy_intp gabor_response_shape[] = { n_imgs, par_T_scales, height, width };		
		gabor_response = PyArray_SimpleNew(4, &gabor_response_shape[0], NPY_DOUBLE);
	}
	else if (n_imgs > 1 && par_T_scales == 1)
	{
		npy_intp gabor_response_shape[] = { n_imgs, height, width };		
		gabor_response = PyArray_SimpleNew(3, &gabor_response_shape[0], NPY_DOUBLE);
	}
	else if (n_imgs == 1 && par_T_scales > 1)
	{
		npy_intp gabor_response_shape[] = { par_T_scales, height, width };		
		gabor_response = PyArray_SimpleNew(3, &gabor_response_shape[0], NPY_DOUBLE);
	}
	else
	{
		npy_intp gabor_response_shape[] = { height, width };		
		gabor_response = PyArray_SimpleNew(2, &gabor_response_shape[0], NPY_DOUBLE);
	}
	
	char * gabor_response_data = ((PyArrayObject*)gabor_response)->data;
	npy_intp gabor_response_stride = ((PyArrayObject*)gabor_response)->strides[((PyArrayObject*)gabor_response)->nd-1];
	
	if (n_imgs > 1)
	{
		DEBMSG("Multiscale gabor filtering over multiple images\n");
		multiscaleGaborFilter_multipleinputs((double*)raw_input_data, (unsigned int)n_imgs, mask_data, (double*)gabor_response_data, (unsigned int)height, (unsigned int)width, (int*)par_T_data, (unsigned int)par_T_scales, (double)par_L, (unsigned int)par_K);
	}
	else
	{
		DEBMSG("Multiscale gabor filtering over a single image\n");
		multiscaleGaborFilter((double*)raw_input_data, mask_data, (double*)gabor_response_data, (unsigned int)height, (unsigned int)width, (int*)par_T_data, (unsigned int)par_T_scales, (double)par_L, (unsigned int)par_K);
	}
	
	return gabor_response;
}
#endif


#ifdef BUILDING_PYTHON_MODULE
PyObject *gaborFilterWithAngles(PyObject *self, PyObject *args)
{
	PyArrayObject *raw_input;
	char *raw_input_data = NULL;
	npy_intp raw_input_stride, n_imgs = 1;
	npy_intp height, width;
	
	PyArrayObject *multiscale_par_T = NULL;
	char * par_T_data = NULL;
	npy_intp par_T_stride, par_T_scales = 1;
	
	PyArrayObject *mask = NULL;
	char * mask_data = NULL;
	npy_intp mask_stride;
	
	double par_L;
	unsigned int par_K;
	
	if (!PyArg_ParseTuple(args, "O!O!dIO!", &PyArray_Type, &raw_input, &PyArray_Type, &multiscale_par_T, &par_L, &par_K, &PyArray_Type, &mask))
	{
		return NULL;
	}
	
	par_T_data = ((PyArrayObject*)multiscale_par_T)->data;
	par_T_scales = ((PyArrayObject*)multiscale_par_T)->dimensions[0];
	par_T_stride = ((PyArrayObject*)multiscale_par_T)->strides[((PyArrayObject*)multiscale_par_T)->nd-1];

	if (((PyArrayObject*)raw_input)->nd > 2)
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
	
	raw_input_data  = ((PyArrayObject*)raw_input)->data;
	raw_input_stride = ((PyArrayObject*)raw_input)->strides[((PyArrayObject*)raw_input)->nd - 1];
	
	mask_data = ((PyArrayObject*)mask)->data;
	mask_stride = ((PyArrayObject*)mask)->strides[((PyArrayObject*)mask)->nd - 1];

	PyObject * gabor_response = NULL;
	PyObject * gabor_angles_response = NULL;
	if (n_imgs > 1 && par_T_scales > 1)	
	{
		npy_intp gabor_response_shape[] = { n_imgs, par_T_scales, height, width };
		gabor_response = PyArray_SimpleNew(4, &gabor_response_shape[0], NPY_DOUBLE);
		gabor_angles_response = PyArray_SimpleNew(4, &gabor_response_shape[0], NPY_DOUBLE);
	}
	else if (n_imgs > 1 && par_T_scales == 1)
	{
		npy_intp gabor_response_shape[] = { n_imgs, height, width };		
		gabor_response = PyArray_SimpleNew(3, &gabor_response_shape[0], NPY_DOUBLE);
		gabor_angles_response = PyArray_SimpleNew(3, &gabor_response_shape[0], NPY_DOUBLE);
	}
	else if (n_imgs == 1 && par_T_scales > 1)
	{
		npy_intp gabor_response_shape[] = { par_T_scales, height, width };		
		gabor_response = PyArray_SimpleNew(3, &gabor_response_shape[0], NPY_DOUBLE);
		gabor_angles_response = PyArray_SimpleNew(3, &gabor_response_shape[0], NPY_DOUBLE);
	}
	else
	{
		npy_intp gabor_response_shape[] = { height, width };		
		gabor_response = PyArray_SimpleNew(2, &gabor_response_shape[0], NPY_DOUBLE);
		gabor_angles_response = PyArray_SimpleNew(2, &gabor_response_shape[0], NPY_DOUBLE);
	}
	
	char * gabor_response_data = ((PyArrayObject*)gabor_response)->data;
	npy_intp gabor_response_stride = ((PyArrayObject*)gabor_response)->strides[((PyArrayObject*)gabor_response)->nd-1];
	
	char * gabor_angles_response_data = ((PyArrayObject*)gabor_angles_response)->data;
	npy_intp gabor_angles_response_stride = ((PyArrayObject*)gabor_angles_response)->strides[((PyArrayObject*)gabor_angles_response)->nd-1];
	
	if (n_imgs > 1)
	{
		DEBMSG("Multiscale gabor filtering over multiple images\n");
		multiscaleGaborFilterWithAngles_multipleinputs((double*)raw_input_data, (unsigned int)n_imgs, mask_data, (double*)gabor_response_data, (double*)gabor_angles_response_data, (unsigned int)height, (unsigned int)width, (int*)par_T_data, (unsigned int)par_T_scales, (double)par_L, (unsigned int)par_K);
	}
	else
	{
		DEBMSG("Multiscale gabor filtering over a single image\n");
		multiscaleGaborFilterWithAngles((double*)raw_input_data, mask_data, (double*)gabor_response_data, (double*)gabor_angles_response_data, (unsigned int)height, (unsigned int)width, (int*)par_T_data, (unsigned int)par_T_scales, (double)par_L, (unsigned int)par_K);
	}
	
	PyObject *gabor_responses_tuple = PyTuple_New(2);
	PyTuple_SetItem(gabor_responses_tuple, 0, gabor_response);
	PyTuple_SetItem(gabor_responses_tuple, 1, gabor_angles_response);
	
	return gabor_responses_tuple;
}
#endif


#ifdef BUILDING_PYTHON_MODULE
static PyMethodDef gabor_methods[] = {
	{ "gaborFilter", gaborFilter, METH_VARARGS, "applies the Gabor filter to the input image, using the parameters t, l and K passed, if the parameter t is a list, then the multiscale Gabor filter is applied instead." },
	{ "gaborFilterWithAngles", gaborFilterWithAngles, METH_VARARGS, "applies the Gabor filter to the input image, using the parameters t, l and K passed, if the parameter t is a list, then the multiscale Gabor filter is applied instead, the angles of the max response are saved too." },
	{ NULL, NULL, 0, NULL }
};
#endif


#ifdef BUILDING_PYTHON_MODULE
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
#endif


#ifdef BUILDING_PYTHON_MODULE
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
#endif
