#include <sys/stat.h>



void c_stat(const char *filename, int *mode, int *exist, int *isdir, int *time) {
    int k;
    struct stat buf;

    k = stat(filename, &buf);
    if (k != 0) {
        *mode = 0;
        *exist = 0;
        *isdir = 0;
        *time = 0;
    } else {
        *mode = buf.st_mode;
        if (*mode == 0) *exist = 0; else *exist = 1;
        if (S_ISDIR(buf.st_mode)) *isdir = 1; else *isdir = 0;
        *time = buf.st_mtime;
    }
}



int c_mkdir(const char *dirname) {
    return mkdir(dirname, 0777);
}
