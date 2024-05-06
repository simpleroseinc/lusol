#pragma once

extern "C" void* create_ctx(int pivot_model);
extern "C" void assign_matrix(void* ctx, int* row_indices, int* col_indices, double* values, int nnz);
extern "C" void factorize_matrix(void* ctx);
extern "C" void col_solve(void* ctx, double* col);
extern "C" void row_solve(void* ctx, double* row);
extern "C" void destroy_ctx(void* ctx);
