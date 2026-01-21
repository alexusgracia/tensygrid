/**
 * @file main.c
 * @author Enrico Uhlenberg (enrico.uhlenberg@gmail.com)
 * @brief Simulation routine for explicit time discrete multilinear models in CPN1 decomposition. (MTI models for short)
 * @date 2025-05-19
 * @details
 *  Very brief theory: <br> <br>
 *      The makeup of an MTI model is comparable to a linear state-space model (LTI). 
 *      A simple  LTI model may consist of two matrices A and B, which describe the linear dynamics of the system. (A describing state-state realtions and B state-input relations)
 *      MTIs extend this concept by allowing multiple states and inputs to have combined influence on the system. In the linear case a difference equation (with states x1,x2 and input u) might look like this:
 *         <br>
 *              \f$x1(t+1) = a*x1(k) + b*x2(k) + c*u(k)\f$
 *          <br>
 *      where a in MTI  with the additional multilinear terms could be:
 *         <br>
 *              \f$x1(t+1) = a*x1(k) + b*x2(k) + c*u(k) + d*x1(k)*x2(k) + e*x2(k)*u(k) + f*u(k)*x1(k)\f$
 *        <br>
 * 
 * Brief description of the operation:
 * Assuming we have system with 
 * - n states (x1,x2,...,xn)
 * - m inputs (u1,u2,...,um)
 * - r rank (basically the Number of different multilinear combinations of states and inputs)
 * our Model is descibed completetly by the matrices 
 * - U ((n+m) x r) containing the information about multilinear relations
 * - Phi (n x r) containing scaling factors, relating the multilinear terms to the individual states
 * 
 * To obtain a next state vector, the following steps are performed:
 * 1. The current state vector is appended to the input vector, creating a combined state-input vector xu.
 * 2. This vector is then multilplied (elementwise) with every column of U
 * 3. For every column j of U, every element if multiplied to create the current multilinear terms.
 * 4. The resulting vector is then right-multiplied with Phi resulting a new state vector.
 * 
 * Note this description is a simplification, detailed descriptions of the algorithm can be found in the implementations.
 * 
 * 
 * The matrices U and Phi can be assumed to be sparse and are stored in CSC format. The main performance gains in this software are achieved by 
 * optimizing memory layout and utilizing SIMD operations with this specific operation in mind.  
 * 
 * Copyright 2025 Enrico Uhlenberg
 * Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:
 * 1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
 * 2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.
 * 3. Neither the name of the copyright holder nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS “AS IS” AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */


#include <stdio.h>
#include <stdlib.h> 
#include <stdint.h>
#include <stdbool.h> 
#include <stdarg.h>
#include <string.h> 
#include <inttypes.h>
#include <time.h>
#include <math.h>
#include <omp.h>
#include <immintrin.h>
#include "mex.h"


// =============================================================================
// Config
// =============================================================================

#define USING_SINGLE_PRECISION 0   ///< Precompiler directive. If set to "1" indicating float, the input data will be converted to float. 
#define DEBUG 0                ///< Precompiler directive. If set to "1" debug messages will be printed to the console. 
#define MEX 1                  ///< Precompiler directive. If set to "1" the code will be compiled as a MEX function for MATLAB.
#define oneMKL 1            ///< Precompiler directive. If set to "1" the code will use oneMKL for SIMD operations.

// =============================================================================
// Constants
// =============================================================================

/**
 * @defgroup PRECISION_SIMD_MACROS
 *  @brief Macros for SIMD operations
 *  @details These macros are used to define the SIMD types and operations based on the selected precision.
 * @{
 */
#if USING_SINGLE_PRECISION
    #define PRECISION float
    #define SIMD_WORD __m256
    #define SIMD_INDICES __m256i
    #define SIMD_LOAD(p) _mm256_loadu_ps((p))
    #define SIMD_STORE(p, v) _mm256_storeu_ps((p), (v))
    #define SIMD_STREAM(p, v) _mm256_stream_ps((p), (v))
    #define SIMD_SET(src,i,ia) _mm256_set_ps((src[(size_t)((i) + 7) * (ia)]), (src[(size_t)((i) + 6) * (ia)]), (src[(size_t)((i) + 5) * (ia)]), (src[(size_t)((i) + 4) * (ia)]),(src[(size_t)((i) + 3) * (ia)]), (src[(size_t)((i) + 2) * (ia)]), (src[(size_t)((i) + 1) * (ia)]), (src[(size_t)(i) * (ia)]))
    #define SIMD_FMADD(a,b,c) _mm256_fmadd_ps((a),(b),(c))
    #define SIMD_GATHER(a,b,c) _mm256_i32gather_ps((a),(b),(c))
    #define SIMD_MUL(a,b) _mm256_mul_ps((a),(b))
    #define SIMD_ADD(a,b) _mm256_add_ps((a),(b))
    #define SIMD_SET1(a) _mm256_set1_ps((a))
    #define SIMD_LOAD_INDICES(p) _mm256_loadu_si256((p))
#else
    #define PRECISION double
    #define SIMD_WORD __m256d
    #define SIMD_INDICES __m128i
    #define SIMD_LOAD(p) _mm256_loadu_pd((p))
    #define SIMD_STORE(p, v) _mm256_store_pd((p), (v))
    #define SIMD_STREAM(p, v) _mm256_stream_pd((p), (v))
    #define SIMD_SET(src,i,ia) _mm256_set_pd((src[(size_t)((i) + 3) * (ia)]), (src[(size_t)((i) + 2) * (ia)]), (src[(size_t)((i) + 1) * (ia)]), (src[(size_t)(i) * (ia)]))
    #define SIMD_FMADD(a,b,c) _mm256_fmadd_pd((a),(b),(c))
    #define SIMD_GATHER(a,b,c) _mm256_i32gather_pd((a),(b),(c))
    #define SIMD_MUL(a,b) _mm256_mul_pd((a),(b))
    #define SIMD_ADD(a,b) _mm256_add_pd((a),(b))
    #define SIMD_SET1(a) _mm256_set1_pd((a))
    #define SIMD_LOAD_INDICES(p) _mm_loadu_si128((p))
#endif
#define ELEMENTS_PER_SIMD_LANE (sizeof(__m256d) / sizeof(PRECISION))

/** @} */



// =========================
// Typedefs
// =========================
/**
 * @brief Struct for storing the sparse matrix in Compressed Sparse Column (CSC) format.
 * @details 
 *      The CSC format stores the non-zeros contiguously in a 1D array, and uses two additional arrays to store the row indices and column pointers.
 *      This is  preferred to the otherwise more commen CSR, since data for this concrete operation is mostly accessed column-wise.
 *      For reference, see https://www.intel.com/content/www/us/en/docs/onemkl/developer-reference-fortran/2023-0/sparse-blas-csc-matrix-storage-format.html
 */
typedef struct __attribute__((aligned(32))) {
    int32_t m;          ///< Number of rows
    int32_t n;          ///< Number of columns
    int32_t *col_ptr;   ///< Column pointers
    int32_t *rows;      ///< Row indices
    int32_t nnz;        ///< Number of non-zero values
    PRECISION *vals;    ///< Contiguous array of non-zero values. Can be double or float
} SparseMatrix_t;




// =========================
// Function Prototypes
// =========================

void run_step();
/**
 * @brief Validate the input arguments for the MEX function.
 * @details Will report errors to matlab if the input arguments are not valid.
 */
bool mexInputValidation(void);

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]);
/**
 * @brief Multiplicatively accumulate the values of a column in a sparse matrix using OpenMP directives.
 * 
 * @param nnz_val Number of non-zero values in the column.
 * @param u_vals_ptr Pointer to the values of the column.
 * @param xu_unpacked_ptr Pointer to the unpacked state-input vector.
 * @return The result of the multiplicative accumulation.
 */
#pragma omp declare simd
static inline PRECISION accumulate_columns_omp(
    int32_t nnz_val,
    const PRECISION* u_vals_ptr,
    const PRECISION* xu_unpacked_ptr
);

/**
 * @brief Gather scattered vector and accumulate the prodcut of that vector and a column in a sparse matrix using AVX intrinsics.
 * 
 * @param nnz Number of non-zero values in the column.
 * @param u Pointer to the values of the column.
 * @param xu Pointer to the unpacked state-input vector.
 * @param index_vector Pointer to the index vector for gathering.
 * @return The result of the multiplicative accumulation.
 */
static inline PRECISION gather_and_accumulate_columns_avx2( 
    size_t nnz,            // Number of non-zero values in current column
    const PRECISION* u,     // Pointer to the values of current column
    const PRECISION* xu,    // Pointer to the unpacked state-input vector
    const int32_t* index_vector // Pointer to the index vector for gathering 
);

/**
 * @brief Copy values from a strided source array to a destination array using AVX2 intrinsics.
 * 
 * @param n Number of elements to copy.
 * @param src Source array.
 * @param ia Stride of the source array.
 * @param dst Destination array.
 */
void strided_gather_avx2(
    size_t n, 
    const PRECISION* src, 
    size_t ia, 
    PRECISION* dst
);

/**
 * @brief Add a scalar value to each element of an input vector and store the result in an output vector using AVX2 intrinsics.
 * 
 * @param n Number of elements in the input vector.
 * @param input_vector Input vector.
 * @param scalar Scalar value to add.
 * @param output_vector Output vector.
 */
void add_scalar_avx2(
    size_t n, 
    const PRECISION* input_vector, 
    PRECISION scalar, 
    PRECISION* output_vector
);


/**
 * @brief Helper function to calculate the padded size of an array for aligned alloc.
 */
static inline size_t calculate_padded_size(
    size_t num_elements, 
    size_t element_size, 
    size_t alignment
);

/**
 * @brief Helper function to allocate aligned memory for an array.
 */
static inline void safe_aligned_alloc(
    void ** ptr,
    size_t num_elements, 
    size_t element_size, 
    size_t alignment, 
    const char* var_name
);


/**
 * @brief Load sparse matrix data from a binary file in CSC format and double precision.
 * 
 * @param filename Name of the file to load.
 * @param mat Pointer to the SparseMatrix_t struct to store the loaded data.
 * @param vals Pointer to the array to store the values.
 * @return true if loading was successful, false otherwise.
 */
void load_sparse_matrix_data(
    const char *filename, 
    SparseMatrix_t *mat
);

/**
 * @brief Allocate memory for the simulation data structures.
 * 
 * @return true if memory allocation was successful, false otherwise.
 */
void allocate_memory(void);

/**
 * @brief Deallocate memory for the simulation data structures.
 * 
 */
void deallocate_memory(void);

/**
 * @brief Run the simulation loop.
 * 
 * This function performs the main simulation loop
 * @return true if the simulation ran successfully, false otherwise.
 */
 void run_sim(void);

/**
 * @brief Compare and display simulation results vs reference for first 2 states and 10 timesteps.
 * 
 * This function compares the simulated results in x_result with the reference results in XRef
 * and displays them in a formatted table for analysis.
 */
void compare_and_display_results(void);

/**
 * @brief Output complete simulation results as CSV to stdout.
 * 
 * Optimized for performance with minimal formatting overhead.
 * Each row represents one timestep, columns represent states.
 */
void output_results_csv_stdout(void);


/**
 * @brief Print an error message and exit the program. 
 * 
 * @param fmt Format string for the error message.
 * @param ... Additional arguments for the format string.
 */
void die(const char * fmt, ...);

/**
 * @brief Convert MATLAB sparse matrix to SparseMatrix_t struct.
 * 
 * @details Memory ownership model for MEX mode:
 * - col_ptr and rows: Copied to mxMalloc'd memory, must be freed with mxFree()
 * - vals: Direct alias to MATLAB-managed memory, never free manually
 * 
 * @param mx_array MATLAB sparse matrix array.
 * @param sparse_mat Pointer to SparseMatrix_t struct to initialize.
 */
void matlab_sparse_to_struct(const mxArray *mx_array, SparseMatrix_t *sparse_mat);

/**
 * @brief Entry point for the MEX function.
 * 
 * @param nlhs Number of left-hand side arguments (outputs).
 * @param plhs Array of left-hand side arguments (outputs).
 * @param nrhs Number of right-hand side arguments (inputs).
 * @param prhs Array of right-hand side arguments (inputs). 
 */

// =============================================================================
// Globals
// =============================================================================

// Structs for sparse matrices. XRef and URef are admittedly dense, but they are stored in a sparse format for consistency
static SparseMatrix_t U; 
static SparseMatrix_t Phi; 
static SparseMatrix_t XRef; 
static SparseMatrix_t URef;

// Number of timesteps to simulate (later on)
static size_t t_num = 60;

// Buffers for storing data without CSC indexing data 
// (either because it's a flat array, or because the structure si identical to an already existing matrix)
static PRECISION *x0                   = NULL; 
static PRECISION *Uxucum               = NULL;
static PRECISION *xu                   = NULL;
static PRECISION *xu_minus_1           = NULL; 
static PRECISION **x_result            = NULL;

int main(int argc, char *argv[])
{
    const char *model_name = "cstr_concat50n400";
    bool output_results = false;
    
    if (argc > 1) {
        int temp_t_num = atoi(argv[1]);
        if (temp_t_num > 0) {
            t_num = (size_t)temp_t_num;
        } else {
            printf("Warning: Invalid t_num value '%s'. Using default value %zu.\n", argv[1], t_num);
        }
    }
    
    if (argc > 2) {
        size_t model_len = strlen(argv[2]);
        if (model_len > 0 && model_len < 100) {
            bool valid_name = true;
            for (size_t i = 0; i < model_len; i++) {
                char c = argv[2][i];
                if (!((c >= 'a' && c <= 'z') || (c >= 'A' && c <= 'Z') || 
                      (c >= '0' && c <= '9') || c == '_')) {
                    valid_name = false;
                    break;
                }
            }
            
            if (valid_name) {
                model_name = argv[2];
            } else {
                printf("Warning: Invalid model name '%s'. Using default '%s'.\n", argv[2], model_name);
            }
        } else {
            printf("Warning: Model name too long or empty. Using default '%s'.\n", model_name);
        }
    }
    
    if (argc > 3 && strcmp(argv[3], "--output") == 0) {
        output_results = true;
    }

    if (!output_results) {
        #if USING_SINGLE_PRECISION
            printf("INFO: Using SINGLE precision computation path.\n");
        #else
            printf("INFO: Using DOUBLE precision computation path.\n");
        #endif
        
        printf("INFO: Using model '%s'\n", model_name);
    }
    
    char file_path[256];
    
    snprintf(file_path, sizeof(file_path), "data/%s_U.bin", model_name);
    load_sparse_matrix_data(file_path, &U);
    
    snprintf(file_path, sizeof(file_path), "data/%s_Phi.bin", model_name);
    load_sparse_matrix_data(file_path, &Phi);
    
    snprintf(file_path, sizeof(file_path), "data/%s_uVec.bin", model_name);
    load_sparse_matrix_data(file_path, &URef);
    
    snprintf(file_path, sizeof(file_path), "data/%s_xVec.bin", model_name);
    load_sparse_matrix_data(file_path, &XRef);

    allocate_memory();

    double begin = omp_get_wtime();

    run_sim(); 
    
    double end = omp_get_wtime();
    
    if (output_results) {
        // Performance mode: only output CSV data
        output_results_csv_stdout();
    } else {
        // Normal mode: existing output
        printf("Main loop finished.\n");
        double time_spent = (double)(end - begin);
        printf("0: %f \n", x_result[0][0]);
        printf("1: %f \n", x_result[1][0]);
        printf("2: %f \n", x_result[2][0]);
        printf("Simulation time: %f seconds.\n", time_spent);
        
        compare_and_display_results();
    }

    return EXIT_SUCCESS;
}

// =============================================================================
// Function Implementations 
// =============================================================================

#pragma omp declare simd
static inline PRECISION accumulate_columns_omp(
    int32_t nnz_val, 
    const PRECISION* u_vals_ptr, 
    const PRECISION* xu_unpacked_ptr 
) {
    PRECISION column_product = 1.0;
    #pragma omp simd reduction(*:column_product)
    for (int32_t k = 0; k < nnz_val; ++k) {
        column_product *= (u_vals_ptr[k] * xu_unpacked_ptr[k] + 1.0);
    }
    return column_product;
}

static inline PRECISION gather_and_accumulate_columns_avx2( 
    size_t nnz,             // Number of non-zero values in current column
    const PRECISION* u,     // Pointer to the values of current column
    const PRECISION* xu,    // Pointer to the unpacked state-input vector
    const int* index_vector // Pointer to the index vector for gathering 
) {
    PRECISION       column_product_scalar = 1.0;                    // Scalar product accumulator 
    SIMD_WORD       product_acc_vec = SIMD_SET1((PRECISION) 1.0);   // SIMD product accumulator
    const SIMD_WORD ones = SIMD_SET1((PRECISION) 1.0);              // Vector of ones for SIMD operations
    size_t         k = 0;                                           // Loop index       
    
    // Buffer for storing the result of the SIMD operations
    PRECISION temp_prod_array[ELEMENTS_PER_SIMD_LANE] __attribute__((aligned(32)));
    
    // Limit for k, because we can only process full SIMD lanes
    size_t limit_avx2 = nnz - (nnz % ELEMENTS_PER_SIMD_LANE); 

    for (; k < limit_avx2; k += ELEMENTS_PER_SIMD_LANE) 
    {
        // Load values from U
        //printf("Loading SIMD vector from u at index %zu\n", k);
        //fflush(stdout); // Ensure the message is printed before any potential crash
        const SIMD_WORD u_vec = SIMD_LOAD(u + k);  
        // Load index vector 
        const SIMD_INDICES inds = SIMD_LOAD_INDICES((const SIMD_INDICES*)(index_vector + k));   
        // Gather scattered elements  
        const SIMD_WORD xu_packed = SIMD_GATHER(xu, inds, sizeof(PRECISION));                   
        // Fused multiply-add operation (U*(xu_minus_one)+1)
        const SIMD_WORD term_result = SIMD_FMADD(u_vec, xu_packed, ones);                       
        // Multiply the result with the accumulator
        product_acc_vec = SIMD_MUL(product_acc_vec, term_result);                               
    }

    // Copy packed vector to memory for scalar multiplication
    SIMD_STORE(temp_prod_array, product_acc_vec);
    
    // Multiply the SIMD lanes to a scalar, 4 or 8 elements depending on precision mode
    for(size_t i = 0; i < ELEMENTS_PER_SIMD_LANE; i++) 
    {
        column_product_scalar *= temp_prod_array[i];
    }

    // Process the remaining elements that don't fit into a full SIMD lane
    for (; k < nnz; k++) 
    {
        column_product_scalar *= (u[k] * xu[index_vector[k]] + 1.0);
    }

    return column_product_scalar;
}

void strided_gather_avx2(
    size_t n, 
    const PRECISION* src, 
    size_t ia, 
    PRECISION* dst
) { 
    size_t i = 0;
    size_t limit_avx2 = n - (n % ELEMENTS_PER_SIMD_LANE);
    for (; i < limit_avx2; i += ELEMENTS_PER_SIMD_LANE) 
    {
        SIMD_WORD a_vec = SIMD_SET(src,i,ia);
        SIMD_STORE(dst + i, a_vec);
    }
    for (; i < n; i++)
    {
        dst[i] = src[(size_t)i * ia];
    } 
}

void add_scalar_avx2(
    size_t n, 
    const PRECISION* input_vector, 
    PRECISION scalar, 
    PRECISION* output_vector
) { 
    size_t i = 0;
    SIMD_WORD scalar_vec = SIMD_SET1(scalar);
    size_t limit_avx2 = n - (n % ELEMENTS_PER_SIMD_LANE);
    
    for (; i < limit_avx2; i += ELEMENTS_PER_SIMD_LANE) 
    {
        SIMD_WORD input_vec = SIMD_LOAD(input_vector + i);
        SIMD_WORD result_vec = SIMD_ADD(input_vec, scalar_vec);
        SIMD_STORE(output_vector + i, result_vec);
    }
    for (; i < n; i++) 
    {
        output_vector[i] = input_vector[i] + scalar;
    }
}




void load_sparse_matrix_data(const char *filename, SparseMatrix_t *mat)
{
    FILE *fp = fopen(filename, "rb");
    if (!fp) 
    { 
        die("Failed to load matrix file %s", filename); 
    }

    mat->m =0; 
    mat->n =0;
    mat->col_ptr = NULL;
    mat->rows = NULL;
    mat->nnz = 0;
    double *vals = NULL; 
    // Read dimensions
    if (
        (fread((void*)&mat->m, sizeof(int32_t), 1, fp) != 1) 
        || 
        (fread((void*)&mat->n, sizeof(int32_t), 1, fp) != 1)
    ) {
        die("Dimension Error in %s", filename); 
    }

    safe_aligned_alloc((void**)&mat->col_ptr ,(size_t)(mat->n) + 1, sizeof(int32_t), 64, "col_ptr");
    if (!mat->col_ptr && ((mat->n) + 1) > 0) 
    { 
        die("Dimension Error in %s", filename);
    }
    // Read column pointers
    if (
        fread((void*)mat->col_ptr, sizeof(int32_t), (size_t)(mat->n) + 1, fp) != (size_t)((mat->n) + 1)
    ) 
    {
        die("Dimension Error in %s", filename);
    }

    mat->nnz = (mat->n >= 0 && mat->col_ptr) ? ((mat->col_ptr)[mat->n] - 1) : -1;
    if (mat->nnz < 0) 
    {
        die("Invalid nnz (%d) for %s\n", mat->nnz, filename);
    }

    safe_aligned_alloc((void**)&mat->rows ,mat->nnz, sizeof(int32_t), 64, "rows");
    safe_aligned_alloc((void**)&vals ,mat->nnz, sizeof(double), 64, "valsdouble"); 

    if (mat->nnz > 0) {
        if (!(mat->rows) || !(vals)) 
        {
            die("Data alloc failed for %s\n", filename);
        }

        // Read row indices and values
        if (
            (fread((void*)mat->rows, sizeof(int32_t), mat->nnz, fp) != (size_t)mat->nnz) 
            ||
            (fread((void*)vals,     sizeof(double),  mat->nnz, fp) != (size_t)mat->nnz)
        ) 
        { 
            die("Data read failed for %s\n", filename);
        }
    }
    fclose(fp);
    
    for(int32_t i = 0; i < mat->n + 1; ++i) 
    {
        mat->col_ptr[i] -= 1; // Convert to 0-based indexing
    }
    for(int32_t i = 0; i < mat->nnz; ++i) 
    {
        mat->rows[i] -= 1; // Convert to 0-based indexing
    }

    #if USING_SINGLE_PRECISION
        safe_aligned_alloc((void**)&mat->vals,mat->nnz, sizeof(float), 64, "U_vals");
        for (int32_t k = 0; k < mat->nnz; ++k) 
        {
            mat->vals[k] = (float)vals[k];
        }
        printf("Data loaded (as double) and converted to float from %s (m=%" PRId32 ", n=%" PRId32 ", nnz=%d).\n", filename, mat->m, mat->n, mat->nnz);
    #else
        mat->vals = vals;
        printf( "Data loaded (as double) from %s (m=%" PRId32 ", n=%" PRId32 ", nnz=%d).\n", filename, mat->m, mat->n, mat->nnz);
    #endif
}

static inline size_t calculate_padded_size(size_t num_elements, size_t element_size, size_t alignment) {
    if (num_elements == 0 || element_size == 0) 
    {
        return 0;
    }

    if (alignment == 0)
    {
        alignment = 1;
    } 
    size_t total_size = num_elements * element_size;
    size_t remainder = total_size % alignment;
    return (remainder == 0) ? total_size : (total_size + alignment - remainder);
}

static inline void safe_aligned_alloc(
    void ** ptr,
    size_t num_elements, 
    size_t element_size, 
    size_t alignment, 
    const char* var_name
) {
    if (num_elements == 0) 
    {
        *ptr = NULL;
    }
    
    if (alignment == 0 || (alignment & (alignment - 1)) != 0 || alignment % sizeof(void*) != 0) 
    {
        die(" Error: Invalid alignment %zu for %s.\n", alignment, var_name);
    }
    size_t padded_size = calculate_padded_size(num_elements, element_size, alignment);
    if (padded_size == 0 && (num_elements * element_size > 0)) 
    {
        die("Error: padded_size is 0 for non-zero request for %s.\n", var_name);
    }
    
    *ptr = NULL;
    #ifdef _WIN32
        *ptr = _aligned_malloc(padded_size, alignment);
    #else
        #if MEX 
        *ptr = mxMalloc(padded_size);
        #else
        if (posix_memalign(ptr, alignment, padded_size) != 0) 
        {
            *ptr = NULL; 
        }
        #endif
    #endif


}

void die(const char * fmt, ...)                           
{                  
    #if MEX
    char errBuf[100];
    va_list ap;                                                                      
    va_start(ap,fmt);   
    sprintf(errBuf,fmt,va_arg(ap,const char *));
    mexErrMsgIdAndTxt("MTISIM:die",errBuf);
    va_end(ap);  
    #else
    va_list ap;                                                                      
    va_start(ap,fmt);         
    vprintf(fmt, ap);                                                        
    va_end(ap);                                                                      
    exit(EXIT_FAILURE);                                                              
    #endif
    
}



void allocate_memory(void)
{

    // Initial state vector
    safe_aligned_alloc((void**)&x0,Phi.m, sizeof(PRECISION), 64, "x0");
   

    // Allocate memory for the accumulation vector and the result table
    safe_aligned_alloc((void**)&Uxucum,U.n, sizeof(PRECISION), 64, "Uxucum");
    safe_aligned_alloc((void**)&x_result,t_num, sizeof(PRECISION*), 64, "x_result table");
    for (size_t i = 0; i < t_num; i++) 
    {
        safe_aligned_alloc((void**)&x_result[i],Phi.m, sizeof(PRECISION), 64, "x_result_row");
        if (x_result[i]) 
        {
            memset(x_result[i], 0, Phi.m * sizeof(PRECISION));
        }
        else
        {
            die("Error allocating memory for x_result[%i]", i);
        }
    }
    safe_aligned_alloc((void**)&xu,U.m, sizeof(PRECISION), 64, "xu");
    
    safe_aligned_alloc((void**)&xu_minus_1,U.m, sizeof(PRECISION), 64, "xu_minus_1");
    
    

    if ((U.nnz > 0 && !U.vals)) {
        die("Memory allocation failed: U.vals is NULL but U.nnz > 0. Exiting.\n");
    }
    if ((Phi.nnz > 0 && !Phi.vals)) {
        die("Memory allocation failed: Phi.vals is NULL but Phi.nnz > 0. Exiting.\n");
    }
    if (!URef.vals) {
        die("Memory allocation failed: URef.vals is NULL. Exiting.\n");
    }
    if (!XRef.vals) {
        die("Memory allocation failed: XRef.vals is NULL. Exiting.\n");
    }
    if (!Uxucum) {
        die("Memory allocation failed: Uxucum is NULL. Exiting.\n");
    }
    if (!x_result) {
        die("Memory allocation failed: x_result is NULL. Exiting.\n");
    }
    if (!xu) {
        die("Memory allocation failed: xu is NULL. Exiting.\n");
    }
    if (!xu_minus_1) {
        die("Memory allocation failed: xu_minus_1 is NULL. Exiting.\n");
    }

}

void run_sim(void)
{
    memcpy(
        x_result[0], 
        XRef.vals, 
        Phi.m * sizeof(PRECISION)
    );
    // Append the current input to the state vector  
    memcpy(
        (xu + (size_t)Phi.m), 
        URef.vals, 
        URef.n * sizeof(PRECISION)
    );
    for (size_t i = 1; i < t_num; i++) 
    {
        // Copy the previous result into the current state vector
        memcpy(xu, x_result[i-1], Phi.m * sizeof(PRECISION));

        // Append the current input to the state vector  
        memcpy(
        (xu + (size_t)Phi.m), 
        URef.vals + ((i-1) * URef.m) , 
        URef.m * sizeof(PRECISION)
        );
        //strided_gather_avx2(
        //    URef.n,                         // Number of elements to gather
        //    URef.vals + (i % URef.m) - 1,   // Pointer to the source array (URef.vals) Modulo URef.m to cycle through inputs 
        //    URef.m,                         // Stride of the source array (URef.m)
        //    (xu + (size_t)XRef.n)           // Pointer to area for inputs in the the destination array (xu)
        //);
        
        // Compute the offset input-state vector
        add_scalar_avx2(U.m, xu, -1.0, xu_minus_1);
        
        // Clear the result vector for good measure
        memset(x_result[i], 0, Phi.m * sizeof(PRECISION));
        
        //Loop through columns
        for (int32_t j = 0; j < U.n; ++j) 
        {
            // col_ptr[j] points to the start of the column j
            // col_ptr[j+1] points to the end of the column j
            const int32_t U_curInd = U.col_ptr[j];                       // Index to start of values of column j
            const int32_t U_nnz = U.col_ptr[j + 1] - U_curInd;           // Number of non-zero values in column j
            const int32_t phi_curInd = Phi.col_ptr[j];                   // Index to start of values of column j in Phi
            const int32_t phi_nnz = Phi.col_ptr[j + 1] - phi_curInd;     // Number of non-zero values in column j in Phi

            Uxucum[j] = gather_and_accumulate_columns_avx2(
                U_nnz,                              // Number of non-zero values in column j
                &U.vals[U_curInd],                  // Pointer to the values of column j 
                xu_minus_1,                         // Pointer to the state-input vector 
                (const int*)&U.rows[U_curInd]       // Pointer to the index vector for gathering 
            );

            // Scale current column j by non-zeros in Phi and additively accumulate the result in x_result[i]
            for(int32_t l = 0; l < phi_nnz; l++) 
            {
                int32_t target_row = Phi.rows[phi_curInd + l];
                if (target_row >= 0 && target_row < Phi.m) 
                {
                   x_result[i][target_row] += Uxucum[j] * Phi.vals[phi_curInd + l];
                }
            }
        }
    }
}

void run_step(void)
{

    memcpy(
        xu, 
        XRef.vals, 
        Phi.m * sizeof(PRECISION)
    );
    // Append the current input to the state vector  
    memcpy(
        (xu + (size_t)Phi.m), 
        URef.vals, 
        URef.m * sizeof(PRECISION)
    );
    
    // Compute the offset input-state vector
    add_scalar_avx2(U.m, xu, -1.0, xu_minus_1);
    
    // Clear the result vector for good measure
    memset(x_result[0], 0, Phi.m * sizeof(PRECISION));
    //Loop through columns
    for (int32_t j = 0; j < U.n; ++j) 
    {
        // col_ptr[j] points to the start of the column j
        // col_ptr[j+1] points to the end of the column j
        const int32_t U_curInd = U.col_ptr[j];                       // Index to start of values of column j
        const int32_t U_nnz = U.col_ptr[j + 1] - U_curInd;           // Number of non-zero values in column j
        const int32_t phi_curInd = Phi.col_ptr[j];                   // Index to start of values of column j in Phi
        const int32_t phi_nnz = Phi.col_ptr[j + 1] - phi_curInd;     // Number of non-zero values in column j in Phi
        Uxucum[j] = gather_and_accumulate_columns_avx2(
            U_nnz,                              // Number of non-zero values in column j
            &U.vals[U_curInd],                  // Pointer to the values of column j 
            xu_minus_1,                         // Pointer to the state-input vector 
            (const int*)&U.rows[U_curInd]       // Pointer to the index vector for gathering 
        );
        // Scale current column j by non-zeros in Phi and additively accumulate the result in x_result[i]
        for(int32_t l = 0; l < phi_nnz; l++) 
        {
            int32_t target_row = Phi.rows[phi_curInd + l];
            if (target_row >= 0 && target_row < Phi.m) 
            {
               x_result[0][target_row] += Uxucum[j] * Phi.vals[phi_curInd + l];
            }
        }
    }
}


void deallocate_memory(void)
{
    // Free the allocated memory for the buffers
    if (x0 != NULL) 
    {
        #if MEX
        mxFree(x0);
        #else
        free(x0); 
        #endif
        x0 = NULL;
    }
    /*
    if (Uxucum != NULL)
    {
        #if MEX
        mxFree(Uxucum);
        #else
        free(Uxucum);
        #endif
        Uxucum = NULL;
        
    } 
    if (xu != NULL)
    {
        #if MEX
        mxFree(xu);
        #else
        free(xu);
        #endif
        xu = NULL;
    } 
    if (xu_minus_1 != NULL) 
    {
        #if MEX
        mxFree(xu_minus_1);
        #else
        free(xu_minus_1);
        #endif
        xu_minus_1 = NULL;
    }
    // Free the result table
    if (x_result != NULL) 
    {
        for (size_t i = 0; i < t_num; i++) 
        {
            if (x_result[i] != NULL) 
            {
                #if MEX
                mxFree(x_result[i]);
                #else
                free(x_result[i]);
                #endif
                x_result[i] = NULL;
            }
        }
        #if MEX
        mxFree(x_result);
        #else
        free(x_result);
        #endif
        x_result = NULL;
    }*/
}

void compare_and_display_results(void)
{
    const size_t num_timesteps = 1;
    const size_t num_states = (size_t)Phi.m;  // Use all states from the system
    
    PRECISION total_squared_error = 0.0;
    size_t total_comparisons = 0;
    
    // Calculate mean error across all states and timesteps
    for (size_t state = 0; state < num_states; state++) {
        for (size_t t = 0; t < num_timesteps && t < t_num; t++) {
            // Get simulated value
            PRECISION sim_val = x_result[t][state];
            
            // Get reference value (stored column-wise in XRef)
            PRECISION ref_val = XRef.vals[state * XRef.m + t];
            
            // Calculate squared error
            PRECISION error = sim_val - ref_val;
            total_squared_error += error * error;
            total_comparisons++;
        }
    }
    
    // Calculate root mean square error
    PRECISION mean_error = sqrt(total_squared_error / total_comparisons);
    
    printf("Mean error: %.6e\n", mean_error);
}

void output_results_csv_stdout(void) 
{
    size_t num_states = (size_t)XRef.n;
    
    setvbuf(stdout, NULL, _IOFBF, 65536);
    
    for (size_t t = 0; t < t_num; t++) {
        for (size_t s = 0; s < num_states; s++) {
            printf("%.15g", x_result[t][s]); 
            if (s < num_states - 1) {
                putchar(',');
            }
        }
        putchar('\n');
    }
    
    fflush(stdout);
}

void matlab_sparse_to_struct(const mxArray *mx_array, SparseMatrix_t *sparse_mat)
{
    // Get matrix dimensions
    mwSize m = mxGetM(mx_array);
    mwSize n = mxGetN(mx_array);
    mwIndex nnz = mxGetNzmax(mx_array);
    
    // Get sparse matrix structure pointers
    mwIndex *ir = mxGetIr(mx_array);    // Row indices
    mwIndex *jc = mxGetJc(mx_array);    // Column pointers
    double *vals = mxGetPr(mx_array);   // Values
    
    // Initialize struct fields
    sparse_mat->m = (int32_t)m;
    sparse_mat->n = (int32_t)n;
    sparse_mat->nnz = (int32_t)nnz;
    
    // Allocate and copy column pointers
    sparse_mat->col_ptr = (int32_t*)mxMalloc((n + 1) * sizeof(int32_t));
    for (mwIndex i = 0; i <= n; i++) {
        sparse_mat->col_ptr[i] = (int32_t)jc[i];
    }
    
    // Allocate and copy row indices
    sparse_mat->rows = (int32_t*)mxMalloc(nnz * sizeof(int32_t));
    for (mwIndex i = 0; i < nnz; i++) {
        sparse_mat->rows[i] = (int32_t)ir[i];
    }
    
    // Values pointer (no copy needed since PRECISION == double)
    sparse_mat->vals = vals;
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{


    if(prhs == NULL || nrhs < 1) {
        mexErrMsgIdAndTxt("MTISIM:mexFunction", "No input provided. Please provide a sparse matrix as input.");
    }
    if(nrhs > 5) {
        mexErrMsgIdAndTxt("MTISIM:mexFunction", "Too many inputs,expected max 5 inputs");
    }
    const mxArray *U_mx = prhs[0];
    const mxArray *phi_mx = prhs[1]; 
    const mxArray *uRef_mx = prhs[2];
    const mxArray *xRef_mx = prhs[3];

    t_num = (nrhs > 4) ? (size_t)mxGetScalar(prhs[4]) : 1; // Default to 1 if not provided 

    if(!mxIsSparse(U_mx) || !mxIsSparse(phi_mx)) {
        mexErrMsgIdAndTxt("MTISIM:mexFunction", "First two inputs must be sparse matrices.");
    }

    matlab_sparse_to_struct(U_mx, &U);
    matlab_sparse_to_struct(phi_mx, &Phi);
    if(mxGetN(xRef_mx) > 1 && mxGetM(xRef_mx) > 1)
    {
        mexErrMsgIdAndTxt("MTISIM:mexFunction", "xRef must be a vector.");
    }
    if(mxGetM(uRef_mx) != U.m - Phi.m) {
        mexErrMsgIdAndTxt("MTISIM:mexFunction", "uRef must have the same number of rows %d as U.m %d - Phi.m %d. Every column should contain one timestep of the input vector.",mxGetM(uRef_mx), U.m, Phi.m);
    } 
    URef.vals = mxGetDoubles(uRef_mx);
    URef.m = (int32_t)mxGetM(uRef_mx);
    URef.n = (int32_t)mxGetN(uRef_mx);
    URef.nnz = (int32_t)mxGetNumberOfElements(uRef_mx);
    
    XRef.vals = mxGetDoubles(xRef_mx);
    XRef.m = (int32_t)mxGetM(xRef_mx);
    XRef.n = (int32_t)mxGetN(xRef_mx);
    XRef.nnz = (int32_t)mxGetNumberOfElements(xRef_mx);
    
    //matlab_sparse_to_struct(uRef_mx, &URef);
    //matlab_sparse_to_struct(xRef_mx, &XRef);

    // Validate all input data before proceeding with simulation
    if (!mexInputValidation()) {
        // TODO
        return;
    }

#if DEBUG
     double begin = omp_get_wtime();
#endif
    
     allocate_memory();
    
    if(t_num == 1)
    {
        run_step();
    }
    else 
    {
        run_sim();
    }
#if DEBUG
    double end = omp_get_wtime();
    double time_spent = (double)(end - begin);
    mexPrintf("Simulation completed in %f seconds.\n", time_spent);
#endif
    // Create output matrix for x_result (t_num x XRef.n)
    if (nlhs > 0) {
        plhs[0] = mxCreateDoubleMatrix((mwSize)Phi.m, (mwSize)t_num, mxREAL);
        double *output_ptr = mxGetPr(plhs[0]);

        // Copy x_result to MATLAB output matrix (column-major order)
        for (size_t t = 0; t < t_num; t++) {
            for (size_t s = 0; s < (size_t)Phi.m; s++) {
                output_ptr[s + t * XRef.n] = x_result[t][s];
                
            }
        }    
    }
    // Clean up allocated memory
    return;
}

bool mexInputValidation(void)
{
    
    if (U.n != Phi.n) {
        mexErrMsgIdAndTxt("MTISIM:validation", "Dimension mismatch: U and Phi must have same number of columns (U.n=%d, Phi.n=%d)", U.n, Phi.n);
        return false;
    }
    
   
    
    #if DEBUG   
    size_t estimated_memory = t_num * XRef.n * sizeof(PRECISION); // x_result memory
    estimated_memory += (U.nnz + Phi.nnz + URef.nnz + XRef.nnz) * sizeof(PRECISION); // matrix values
    if (estimated_memory > 2LL * 1024 * 1024 * 1024) { // 2GB limit
        mexErrMsgIdAndTxt("MTISIM:validation", "Estimated memory usage (%.2f GB) exceeds safety limit of 2GB", 
                         estimated_memory / (1024.0 * 1024.0 * 1024.0));
        return false;
    }
    
    mexPrintf("Input validation passed:\n");
    mexPrintf("  U: %dx%d, nnz=%d\n", U.m, U.n, U.nnz);
    mexPrintf("  Phi: %dx%d, nnz=%d\n", Phi.m, Phi.n, Phi.nnz);
    mexPrintf("  uRef: %dx%d, nnz=%d\n", URef.m, URef.n, URef.nnz);
    mexPrintf("  xRef: %dx%d, nnz=%d\n", XRef.m, XRef.n, XRef.nnz);
    mexPrintf("  Timesteps: %zu\n", t_num);
    mexPrintf("  Estimated memory: %.2f MB\n", estimated_memory / (1024.0 * 1024.0));
    #endif
    
    return true;
}
/** @file */
