//
// File:       hw05.c
//
// Abstract:   covert grey scale of image using OpenCL which
//             replace r, g, b with a new L lunimance 
//             
//

////////////////////////////////////////////////////////////////////////////////

#include <fcntl.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include <math.h>
#include <unistd.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <CL/opencl.h>
#include "im1.h"

////////////////////////////////////////////////////////////////////////////////

// Use a static data size for simplicity
//
////////////////////////////////////////////////////////////////////////////////

// Simple compute kernel which computes the square of an input array 
//
const char *KernelSource = "\n" \
"__kernel void rgb2grey(								 \n" \
"  __global float *input,							   	 \n" \
"  __global float *output,							   	 \n" \
"  const int w,										 \n" \
"  const int h)										 \n" \
"{										   	 \n" \
"    size_t id0 = get_global_id( 0 ) ;						   	 \n" \
"    size_t id1 = get_global_id( 1 );						   	 \n" \
"    if ( id0 < 0 || id1 < 0 || id0 >= w || id1 >= h ) return;				 \n" \
"    unsigned int idx = ((id1 * w) + id0) * 3;					   	 \n" \
"    float L = 0.2126f * input[idx] + 0.7152f * input[idx+1] + 0.0722f * input[idx+2];	 \n" \
"    output[idx] = L;									 \n" \
"    output[idx+1] = L;		     							 \n" \
"    output[idx+2] = L;		     							 \n" \
"}											 \n" \
"\n";

////////////////////////////////////////////////////////////////////////////////

int main(int argc, char** argv)
{
    int err;                            // error code returned from api calls
      
    float *data;			// original data set given to device
    float *results;		        // results returned from device
    
    cl_device_id device_id[2];          // compute device id
    cl_context context;                 // compute context
    cl_command_queue commands;          // compute command queue
    cl_program program;                 // compute program
    cl_kernel kernel;                   // compute kernel
    
    cl_mem input;                       // device memory used for the input array
    cl_mem output;                      // device memory used for the output array
    
    
    // load data
    printf("reading openEXR file %s\n", argv[1]);
    int w, h;			       // the width & height of the image, used frequently!
    readOpenEXRFile (argv[1], &data, w, h);

    // set 
    size_t global[2] ;	       // global domain size for our calculation
    size_t local[2]  = {32, 32};                  // local work group size for our calculation

   
    //
    // Get platform info
    //
    cl_platform_id platforms[32];
    cl_uint num_platforms;
    clGetPlatformIDs (32, platforms, &num_platforms );
    
    if (num_platforms == 0) {
        printf("Error: cant find any platforms!\n");
        return EXIT_FAILURE;
    }
    else {
        printf ("Got %d platforms to choose from\n", num_platforms);
    }

     //
    // how many devices are there?
    //
    cl_uint num_devs = 0;
    err = clGetDeviceIDs(platforms[0], CL_DEVICE_TYPE_ALL, 1, NULL, &num_devs);
    if (err != CL_SUCCESS)
    {
        printf("Error: Failed to create a device group!\n");
        return EXIT_FAILURE;
    }
    else {
        printf ("Got %d devices to choose from\n", num_devs);
    }
    

    //
    // What kind of device are we interested in?
    //
    bool use_gpu = true; // specify whether to use CPU or GPU
    
    if (use_gpu)
        printf("Using GPU\n");
    else
        printf("Using CPU\n");
    
    
    cl_uint num_devs_of_class = 1;
    if (use_gpu)
        num_devs_of_class = 2;

    
    //
    // Get a compute device:
    //
    err = clGetDeviceIDs(platforms[0], use_gpu ? CL_DEVICE_TYPE_GPU : CL_DEVICE_TYPE_CPU, num_devs_of_class, device_id, NULL);
    if (err != CL_SUCCESS)
    {
        printf("Error: Failed to create a device group!\n");
        return EXIT_FAILURE;
    }
 
    cl_uint whichdev = 0;
    
    if (use_gpu) {
        whichdev = 0;  // 0 or 1 when choosing which GPU to use
    }
    
    printf("Choosing  id %d \n", whichdev);
    
    //
    // Create a compute context 
    //
    context = clCreateContext(0, 1, &device_id[whichdev], NULL, NULL, &err);
    if (!context)
    {
        printf("Error: Failed to create a compute context!\n");
        return EXIT_FAILURE;
    }

    
    //
    // Create a command queue for the device
    //
    commands = clCreateCommandQueue(context, device_id[whichdev], 0, &err);
    if (!commands)
    {
        printf("Error: Failed to create a command commands!\n");
        return EXIT_FAILURE;
    }
    
    //
    // Create the compute program from the source buffer
    //
    program = clCreateProgramWithSource(context, 1, (const char **) & KernelSource, NULL, &err);
    if (!program)
    {
        printf("Error: Failed to create compute program!\n");
        return EXIT_FAILURE;
    }

    
    //
    // Build the program executable
    //
    err = clBuildProgram(program, 0, NULL, NULL, NULL, NULL);
    if (err != CL_SUCCESS)
    {
        size_t len;
        char buffer[2048];

        printf("Error: Failed to build program executable!\n");
        clGetProgramBuildInfo(program, device_id[whichdev], CL_PROGRAM_BUILD_LOG, sizeof(buffer), buffer, &len);
        printf("%s\n", buffer);
        exit(1);
    }

    
    //
    // Create the compute kernel in the program we wish to run
    //
    kernel = clCreateKernel(program, "rgb2grey", &err);
    if (!kernel || err != CL_SUCCESS)
    {
        printf("Error: Failed to create compute kernel!\n");
        exit(1);
    }

    
    //
    // Create the input and output arrays in device memory for our calculation
    //
    input = clCreateBuffer(context,  CL_MEM_READ_ONLY,  sizeof(float) * 3 * w * h, NULL, NULL);
    output = clCreateBuffer(context, CL_MEM_WRITE_ONLY, sizeof(float) * 3 * w * h, NULL, NULL);
    if (!input || !output)
    {
        printf("Error: Failed to allocate device memory!\n");
        exit(1);
    }    
    
    
    //
    // Write our data set into the input array in device memory 
    //
    err = clEnqueueWriteBuffer(commands, input, CL_TRUE, 0, sizeof(float) * 3 * w * h, data, 0, NULL, NULL);
    if (err != CL_SUCCESS)
    {
        printf("Error: Failed to write to source array!\n");
        exit(1);
    }

    
    //
    // Set the arguments to our compute kernel
    //
    err = 0;
    err  = clSetKernelArg(kernel, 0, sizeof(cl_mem), &input);
    err |= clSetKernelArg(kernel, 1, sizeof(cl_mem), &output);
    err |= clSetKernelArg(kernel, 2, sizeof(int), &w);
    err |= clSetKernelArg(kernel, 3, sizeof(int), &h);
    if (err != CL_SUCCESS)
    {
        printf("Error: Failed to set kernel arguments! %d\n", err);
        exit(1);
    }

    
    //
    // Get the maximum work group size for executing the kernel on the device
    //
   // err = clGetKernelWorkGroupInfo(kernel, device_id[whichdev], CL_KERNEL_WORK_GROUP_SIZE, sizeof(local), local, NULL);
   // if (err != CL_SUCCESS)
   // {
   //     printf("Error: Failed to retrieve kernel work group info! %d\n", err);
   //     exit(1);
   // }
   // printf("Got %ld, %ld  for the maximum work group size\n", local[0], local[1]);

    
    //
    // Execute the kernel over the entire range of our 1d input data set
    // using the maximum number of work group items for this device
    //
    global[0] = ceil(w/ 32.0) * 32; 
    global[1] = ceil(h /32.0) * 32;
    printf ("global dim %ld, %ld \n", global[0], global[1]);
     
    err = clEnqueueNDRangeKernel(commands, kernel, 2, NULL, global, local, 0, NULL, NULL);
    if (err)
    {
        printf("Error: Failed to execute kernel!\n");
        return EXIT_FAILURE;
    }

    
    //
    // Wait for the command commands to get serviced before reading back results
    //
    clFinish(commands);
    
    
    //
    // Read back the results from the device to verify the output
    //
    err = clEnqueueReadBuffer( commands, output, CL_TRUE, 0, sizeof(float) * 3 * w * h, results, 0, NULL, NULL );  
    if (err != CL_SUCCESS)
    {
        printf("Error: Failed to read output array! %d\n", err);
        exit(1);
    }
    
    
    //
    // Check the results
    //
    printf("writing output image hw05_openCL.exr\n");
    writeOpenEXRFile ("hw05_openCL.exr", results, w, h);

    // Shutdown and cleanup
    //
    clReleaseMemObject(input);
    clReleaseMemObject(output);
    clReleaseProgram(program);
    clReleaseKernel(kernel);
    clReleaseCommandQueue(commands);
    clReleaseContext(context);

    return 0;
}

