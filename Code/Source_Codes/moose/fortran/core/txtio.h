#define WRITETXT(fmt, arg) \
    write (unit, fmt, iostat=iostat, iomsg=iomsg) arg; \
    if (iostat /= 0) return

#define NEWRECORD() \
    write (unit, '(/)', iostat=iostat, iomsg=iomsg); \
    if (iostat /= 0) return
