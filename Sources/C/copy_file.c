///  @file copy_file.c
///  @brief Functions to fast copy files.
///
///  MPI and calls to fork do not mix well. In place of slow copies this makes
///  Kernel function calls to the OS fast copy systems by passing a fork to the
///  cp command.
//******************************************************************************

#include <sys/errno.h>

#if defined(__APPLE__)
#include <copyfile.h>
#else

#include <fcntl.h>
#include <sys/sendfile.h>
#include <sys/stat.h>

///  Flag to copy the enture file.
#define COPYFILE_ALL 0

//------------------------------------------------------------------------------
///  @brief Copy a file from the source to the destination.
///
///  Linux does not define the copyfile function. Define one using the
///  <a href="http://stackoverflow.com/questions/2180079/how-can-i-copy-a-file-on-unix-using-c/2180157#2180157">sendfile(2)</a>
///  function.
///
///  @param[in] src      Path to the source file to copy.
///  @param[in] dest     Path to the destination file to copy.
///  @param[in] unused_1 Unused argument.
///  @param[in] unused_2 Unused argument.
///  @return Error status.
//------------------------------------------------------------------------------
size_t copyfile(const char *src, const char *dest, const int unused_1, const int unused_2) {
    
    // Open source file as read only.
    int read_file = open(src, O_RDONLY);

    // Need the permissions and file size values of the source file.
    struct stat stat_buffer;
    fstat(read_file, &stat_buffer);

    // Open the destination file as write only with the same permissions as the
    // source file.
    int write_file = open(dest, O_WRONLY | O_CREAT, stat_buffer.st_mode);
    
    // Copy the file.
    off_t bytesCopied = 0;
    size_t error = sendfile(write_file, read_file, &bytesCopied, stat_buffer.st_size);
    
    // Close file handles.
    close(read_file);
    close(write_file);
}

#endif

//------------------------------------------------------------------------------
///  @brief Copy a file from the source to the destination.
///
///  Copy file returns the error status. Check errno to see the error that
///  occurred.
///
///  @param[in] src  Path to the source file to copy.
///  @param[in] dest Path to the destination file to copy.
///  @return Error code.
//------------------------------------------------------------------------------
int copy_file_c(const char *src, const char *dest) {
    if (copyfile(src, dest, 0, COPYFILE_ALL) < 0) {
        return errno;
    } else {
        return 0;
    }
}
