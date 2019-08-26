/**\file
 * Read text data file into memory as an array of numbers in double format.
 */
#ifndef __MREADARRAY_H__
#define __MREADARRAY_H__

/** Handle structure for mrdary.
 * We organize arrays in the C row major fashion.
 * If the handle is not properly initialized or is freed, marray==NULL
 */
typedef struct mrdary_handle
{
    FILE *fp;
    char *linebuf;
    size_t column;
    size_t row;
    size_t rowmax;
    double *marray;
    double *min;
    double *max;
} mrdary_hdl;

/** initialize with a file.  Number of columns is automatically determined */
mrdary_hdl *mrdary_init_f(const char *fname, size_t rowmax);
/** Read the entire file at once.
 * @return The number of rows read.
 */
size_t mrdary_read_all(mrdary_hdl *hdl);
/** Free the handle */
int mrdary_free(mrdary_hdl *hdl);
/** For get and set a value in the array.
 * @return A pointer to the memory at row m and column n.
 */
double *mrdary_value_mn(mrdary_hdl *hdl, size_t m, size_t n);
/** Get min value of a column.
 * @return A pointer to the min value of column i.
 */
double *mrdary_min(mrdary_hdl *hdl, size_t i);
/** Get max value of a column.
 * @return A pointer to the max value of column i.
 */
double *mrdary_max(mrdary_hdl *hdl, size_t i);

#endif /* __MREADARRAY_H__ */
