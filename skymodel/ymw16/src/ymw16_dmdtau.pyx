cdef extern from "cn.h":
     void dmdtau(double gl, double gb ,double dordm, double DM_Host, int ndir, int np, int vbs, char *dirname, char *text, double *dmord)



# np=1 for Galactic pulsars
# vbs=0 for no verbosity
# dirname is location of data files
# text="" for no extra display text

def dmdtau_c(gl, gb, dordm, ndir, dirname):
    cdef double dmord

    cdef char* dirname_c = dirname
    cdef char* text_c = ''
    np=1
    vbs=0
    DM_Host=0

    dmdtau(gl, gb, dordm, DM_Host, ndir, np, vbs, dirname_c, text_c, &dmord)

    return dmord
