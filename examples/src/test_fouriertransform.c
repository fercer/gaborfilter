#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <fftw3.h>

int main(int argc, char * argv[])
{    
    if (argc < 2)
    {
        printf("\nExample: Filter an input image using the Gabor filter\n");
        printf("First argument is the path to a .PGM image in raw format\n");
        return -1;
    }
  
    FILE * fp_img = fopen(argv[1], "r");
    if (!fp_img)
    {
        printf("<<Error: The file \"%s\" could not be loaded>>\n", argv[1]);
        return -1;
    }
    
    char magic_number[3];
    fgets(magic_number, 3, fp_img);
    if (strcmp(magic_number, "P5"))
    {
        printf("<<Error: The input file is not of the required format, its format is: \"%s\" instead of \"P5\">>\n", magic_number);
        fclose(fp_img);
        return -1;
    }
    fgetc(fp_img);
    
    char commentary[64];
    fgets(commentary, 64, fp_img);
    printf("Commentary: %s\n", commentary);
    int height, width;
    if (*commentary == '#')
    {
        fscanf(fp_img, "%i", &height);
        fscanf(fp_img, "%i", &width);
    }
    else
    {
        height = atoi(commentary);
        char * width_ptr = strchr(commentary, ' ') + 1;
        width = atoi(width_ptr);
    }
    
    int max_value;
    fscanf(fp_img,"%i", &max_value);
    
    printf("The image has dimension: %ix%i pixels, with maximum value of: %i\n", height, width, max_value);
    
    int read_value;
    long double * img = (long double*) malloc(height * width * sizeof(long double));
    fftwl_complex* fimg = (fftwl_complex*) fftwl_malloc(height * (width/2 + 1) * sizeof(fftwl_complex));
    long double * ifimg = (long double*) malloc(height * width * sizeof(long double));

    for (unsigned int x = 0; x < width; x++)
    {
        for (unsigned int y = 0; y < height; y++)
        {
            read_value = fgetc(fp_img);
            *(img + y*width + x) = 1.0 - (long double)read_value / 255.0;
        }
    }
    
    fclose(fp_img);

    printf("Image loaded ...\n");
    
    fftwl_plan forward_plan;
	
    /* Initialize the function: */
    printf("Preparing to perform the Fourier transform ...\n");
    forward_plan = fftwl_plan_dft_r2c_2d(height, width, img, fimg, FFTW_ESTIMATE);
    printf("Ready to perform the Fourier transform ...\n");
    fftwl_execute(forward_plan);
    printf("Fourier transform applied successfully ...\n");
    fftwl_destroy_plan(forward_plan);
    
    free(img);
    
    fp_img = fopen("fouriertransform.bin", "wb");
    
    fwrite(&height, sizeof(int), 1, fp_img);
    fwrite(&width, sizeof(int), 1, fp_img);
    
    double real_part, imag_part;
    for (unsigned int y=0; y<height; y++)
    {
        for (unsigned int x=0; x<(width/2+1); x++)
        {        
            real_part = (double)*(*(fimg + y*(width/2 + 1) + x));
            fwrite(&real_part, sizeof(double), 1, fp_img);
        }
    }
    for (unsigned int y=0; y<height; y++)
    {
        for (unsigned int x=0; x<(width/2+1); x++)    
        {        
            imag_part = (double)*(*(fimg + y*(width/2 + 1) + x)+1);
            fwrite(&imag_part, sizeof(double), 1, fp_img);
        }
    }
    
    fclose(fp_img);

    fftwl_plan backward_plan;
    /* Initialize the function: */
    printf("Preparing to perform the inverse Fourier transform ...\n");
    backward_plan = fftwl_plan_dft_c2r_2d(height, width, fimg, ifimg, FFTW_ESTIMATE);
    printf("Ready to perform the inverse Fourier transform ...\n");
    fftwl_execute(backward_plan);
    printf("Inverse Fourier transform applied successfully ...\n");

    fftwl_free(fimg);
    
    fp_img = fopen("ift.bin", "wb");
    
    fwrite(&height, sizeof(int), 1, fp_img);
    fwrite(&width, sizeof(int), 1, fp_img);
    for (unsigned int xy; xy < height*width; xy++)
    {
        real_part = (double) *(ifimg+xy);
        fwrite(&real_part, sizeof(double), 1, fp_img);
    }
    
    fclose(fp_img);
    
    free(ifimg);
    
    return 0;
}
