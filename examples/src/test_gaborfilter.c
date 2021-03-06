#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <gabor.h>

int main(int argc, char * argv[])
{    
    if (argc < 2)
    {
        printf("\nExample: Filter an input image using the Gabor filter\n");
        printf("First argument is the path to a .PGM image in raw format\n");
        printf("Second to fourth arguments are the parameters for the Gabor filter\n");
        return -1;
    }
 
    int par_K = 45;
    double par_L = 2.65;
    int par_T = 13;
 
    FILE * fp_img = fopen(argv[1], "r");
    if (!fp_img)
    {
        printf("<<Error: The file \"%s\" could not be loaded>>\n", argv[1]);
        return -1;
    }
    
    char magic_number[3];
    fgets(magic_number, 3, fp_img);
    if (strcmp(magic_number, "P5") && strcmp(magic_number, "P2"))
    {
        printf("<<Error: The input file is not of the required format, its format is: \"%s\" instead of \"P5/P2\">>\n", magic_number);
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
    char * mask = (char*) malloc(height * width * sizeof(char));
    double * img = (double*) malloc(height * width * sizeof(double));
    double * resp = (double*) malloc(height * width * sizeof(double));

    if (magic_number[1] == '2')
    {
        printf("Reading the image from an ascii file\n");
        for (unsigned int xy = 0; xy < height*width; xy++)
        {
            fscanf(fp_img, "%i", &read_value);
            *(img + xy) = 1.0 - (double)read_value / 255.0;
            *(mask+xy) = (char)1;
        }
    }
    else
    {
        printf("Reading the image from a raw file\n");
        for (unsigned int xy = 0; xy < height*width; xy++)
        {
            read_value = fgetc(fp_img);
            *(img + xy) = 1.0 - (double)read_value / 255.0;
            *(mask+xy) = (char)1;
        }
    }
    
    fclose(fp_img);

    printf("Image loaded (first, middle and last values are: %.16f, %.16f, %.16f)...\n", *img, *(img+(150-1)*300+150-1), *(img + height*width-1));

    singleScaleGaborFilter(img, mask, resp, height, width, par_T, par_L, par_K);

    free(img);
    free(mask);
    
    fp_img = fopen("gabor_resp.bin", "wb");
    
    fwrite(&height, sizeof(int), 1, fp_img);
    fwrite(&width, sizeof(int), 1, fp_img);
    fwrite(resp, sizeof(double), height*width, fp_img);
    
    fclose(fp_img);

    free(resp);   
    
    printf("Example run successfully ...\n");
    
    return 0;
}
