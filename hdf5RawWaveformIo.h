#ifndef __HDF5IO_H__
#define __HDF5IO_H__

#include <hdf5.h>

#define HDF5IO_NAME_BUF_SIZE 256

struct hdf5io_waveform_file 
{
    hid_t waveFid;
    size_t nPt;
    size_t nCh;
    size_t nWfmPerChunk;
    size_t nEvents;
};

struct hdf5io_waveform_event
{
    size_t eventId;
    /* wavBuf should point to a contiguous 2D array, mapped as
     * ch1..ch2..ch3..ch4 (row-major).  Omitting one or more ch? is
     * allowed in accordance with chMask.*/
    char *wavBuf;
};

/* nWfmPerChunk: waveforms are stored in 2D arrays.  To optimize
 * performance, n waveforms are grouped together to be put in the same
 * array, then the (n+1)th waveform is put into the next grouped
 * array, and so forth. */
struct hdf5io_waveform_file *hdf5io_open_file(const char *fname, size_t nWfmPerChunk,
                                              size_t nCh);
struct hdf5io_waveform_file *hdf5io_open_file_for_read(const char *fname);
int hdf5io_close_file(struct hdf5io_waveform_file *wavFile);
/* flush also writes nEvents to the file */
int hdf5io_flush_file(struct hdf5io_waveform_file *wavFile);

int hdf5io_write_waveform_attribute_in_file_header(struct hdf5io_waveform_file *wavFile,
                                                   struct waveform_attribute *wavAttr);
int hdf5io_read_waveform_attribute_in_file_header(struct hdf5io_waveform_file *wavFile,
                                                   struct waveform_attribute *wavAttr);
int hdf5io_write_event(struct hdf5io_waveform_file *wavFile,
                       struct hdf5io_waveform_event *wavEvent);
int hdf5io_read_event(struct hdf5io_waveform_file *wavFile,
                      struct hdf5io_waveform_event *wavEvent);
size_t hdf5io_get_number_of_events(struct hdf5io_waveform_file *wavFile);

#endif /* __HDF5IO_H__ */
