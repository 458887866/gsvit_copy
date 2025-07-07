bin_PROGRAMS=gsvit2d
if USE_CUDA
gsvit2d_SOURCES=\
   1dpool.c \
   1dpool.h \
   boundary.c \
   boundary.h \
   constants.h \
   dcube.c \
   dcube.h \
   dline.c \
   dline.h \
   farfield.c \
   farfield.h \
   fcube.c \
   fcube.h \
   gpu.c \
   gpu.h \
   gpu_kernels.cu \
   gpu_kernels.h \
   icube.c \
   icube.h \
   output.c \
   output.h \
   plan.c \
   plan.h \
   pool.c \
   pool.h \
   settings.c \
   settings.h \
   source.c \
   source.h \
   svit2d.c 
else
gsvit2d_SOURCES=\
   1dpool.c \
   1dpool.h \
   boundary.c \
   boundary.h \
   constants.h \
   dcube.c \
   dcube.h \
   dline.c \
   dline.h \
   farfield.c \
   farfield.h \
   fcube.c \
   fcube.h \
   gpu.h \
   icube.c \
   icube.h \
   output.c \
   output.h \
   plan.c \
   plan.h \
   pool.c \
   pool.h \
   settings.c \
   settings.h \
   source.c \
   source.h \
   svit2d.c 
endif


#if USE_CUDA
#xsvit_SOURCES=\
#   1dpool.c \
#   1dpool.h \
#   boundary.c \
#   boundary.h \
#   constants.h \
#   dcube.c \
#   dcube.h \
#   dline.c \
#   dline.h \
#   farfield.c \
#   farfield.h \
#   fcube.c \
#   fcube.h \
#   gpu.c \
#   gpu.h \
#   gpu_kernels.cu \
#   gpu_kernels.h \
#   icube.c \
#   icube.h \
#   output.c \
#   output.h \
#   pfarfield.c \
#   pfarfield.h \
#   plan.c \
#   plan.h \
#   pool.c \
#   pool.h \
#   settings.c \
#   settings.h \
#   source.c \
#   source.h \
#   xsvit.h \
#   xsvit.c \
#   xsvitdialogs.h \
#   xsvitdialogs.c \
#   tests.c \
#   zpool.c \
#   zpool.h 
#else
#xsvit_SOURCES=\
#   1dpool.c \
#   1dpool.h \
#   boundary.c \
#   boundary.h \
#   constants.h \
#   dcube.c \
#   dcube.h \
#   dline.c \
#   dline.h \
#   farfield.c \
#   farfield.h \
#   fcube.c \
#   fcube.h \
#   gpu.h \
#   icube.c \
#   icube.h \
#   output.c \
#   output.h \
#   pfarfield.c \
#   pfarfield.h \
#   plan.c \
#   plan.h \
#   pool.c \
#   pool.h \
#   settings.c \
#   settings.h \
#   source.c \
#   source.h \
#   xsvit.h \
#   xsvit.c \
#   xsvitdialogs.h \
#   xsvitdialogs.c \
#   tests.h \
#   tests.c \
#   zpool.c \
#   zpool.h
#endif
  
if USE_CUDA
gsvit2d_LDADD=$(GLIB_LIBS) $(GWYDDION_LIBS) $(CUDA_LIBS)
gsvit2d_CFLAGS=$(GLIB_CFLAGS) $(CUDA_CFLAGS) $(GWYDDION_CFLAGS) $(UCUDA_CFLAGS) $(GTK_CFLAGS)
else
gsvit2d_LDADD=$(GLIB_LIBS) $(GWYDDION_LIBS)
gsvit2d_CFLAGS=$(GLIB_CFLAGS) $(GWYDDION_CFLAGS) $(UCUDA_CFLAGS) $(GTK_CFLAGS) 
endif

#xsvit2d_LDADD=$(gsvit2d_LDADD)
#xsvit2d_CFLAGS=$(gsvit2d_CFLAGS)

.cu.o:
	$(NVCC) -o $@ -c $< $(NVCCFLAGS)
.cu.lo:
	$(top_srcdir)/cudalt.py $@ $(NVCC) -c $(NVCCFLAGS) $<

