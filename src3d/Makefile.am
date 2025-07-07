bin_PROGRAMS=gsvit xsvit
if USE_CUDA
gsvit_SOURCES=\
   1dpool.c \
   1dpool.h \
   boundary.c \
   boundary.h \
   constants.h \
   dcube.c \
   dcube.h \
   defaults.h \
   dline.c \
   dline.h \
   farfield.c \
   farfield.h \
   fcube.c \
   fcube.h \
   global.c \
   global.h \
   gpu.c \
   gpu.h \
   gpu_kernels.cu \
   gpu_kernels.h \
   icube.c \
   icube.h \
   messages.h \
   modifiers.c \
   modifiers.h \
   output.c \
   output.h \
   pfarfield.c \
   pfarfield.h \
   plan.c \
   plan.h \
   pool.c \
   pool.h \
   settings.c \
   settings.h \
   source.c \
   source.h \
   svit.c \
   subgrid.c \
   subgrid.h \
   tests.c \
   yee.c \
   yee.h \
   zpool.c \
   zpool.h 
else
gsvit_SOURCES=\
   1dpool.c \
   1dpool.h \
   boundary.c \
   boundary.h \
   constants.h \
   dcube.c \
   dcube.h \
   defaults.h \
   defaults.h \
   dline.c \
   dline.h \
   farfield.c \
   farfield.h \
   fcube.c \
   fcube.h \
   global.c \
   global.h \
   gpu.h \
   icube.c \
   icube.h \
   messages.h \
   modifiers.c \
   modifiers.h \
   output.c \
   output.h \
   pfarfield.c \
   pfarfield.h \
   plan.c \
   plan.h \
   pool.c \
   pool.h \
   settings.c \
   settings.h \
   source.c \
   source.h \
   svit.c \
   subgrid.c \
   subgrid.h \
   tests.h \
   tests.c \
   yee.c \
   yee.h \
   zpool.c \
   zpool.h
endif

if USE_CUDA
xsvit_SOURCES=\
   1dpool.c \
   1dpool.h \
   boundary.c \
   boundary.h \
   constants.h \
   dcube.c \
   dcube.h \
   defaults.h \
   dline.c \
   dline.h \
   farfield.c \
   farfield.h \
   fcube.c \
   fcube.h \
   global.c \
   global.h \
   gpu.c \
   gpu.h \
   gpu_kernels.cu \
   gpu_kernels.h \
   icube.c \
   icube.h \
   messages.h \
   modifiers.c \
   modifiers.h \
   output.c \
   output.h \
   pfarfield.c \
   pfarfield.h \
   plan.c \
   plan.h \
   pool.c \
   pool.h \
   settings.c \
   settings.h \
   source.c \
   source.h \
   subgrid.c \
   subgrid.h \
   write.c\
   write.h\
   xs_dialogs.h \
   xs_dialogs.c \
   xs_layer_treeview.c \
   xs_layer_treeview.h \
   xs_main_menu.c \
   xs_main_menu.h \
   xs_mat_controls.c\
   xs_mat_controls.h\
   xs_mat_tree.c\
   xs_mat_tree.h\
   xs_par_controls.c\
   xs_par_controls.h\
   xs_par_tree.c\
   xs_par_tree.h\
   xs_scene.c\
   xs_scene.h\
   xsvit.c \
   xsvit.h \
   tests.c \
   yee.c \
   yee.h \
   zpool.c \
   zpool.h 
else
xsvit_SOURCES=\
   1dpool.c \
   1dpool.h \
   boundary.c \
   boundary.h \
   constants.h \
   dcube.c \
   dcube.h \
   defaults.h \
   dline.c \
   dline.h \
   farfield.c \
   farfield.h \
   fcube.c \
   fcube.h \
   global.c \
   global.h \
   gpu.h \
   icube.c \
   icube.h \
   messages.h \
   modifiers.c \
   modifiers.h \
   output.c \
   output.h \
   pfarfield.c \
   pfarfield.h \
   plan.c \
   plan.h \
   pool.c \
   pool.h \
   settings.c \
   settings.h \
   source.c \
   source.h \
   subgrid.c \
   subgrid.h \
   write.c\
   write.h\
   xs_dialogs.h \
   xs_dialogs.c \
   xs_layer_treeview.c \
   xs_layer_treeview.h \
   xs_main_menu.c \
   xs_main_menu.h \
   xs_mat_controls.c\
   xs_mat_controls.h\
   xs_mat_tree.c\
   xs_mat_tree.h\
   xs_par_controls.c\
   xs_par_controls.h\
   xs_par_tree.c\
   xs_par_tree.h\
   xs_scene.c\
   xs_scene.h\
   xsvit.h \
   xsvit.c \
   tests.h \
   tests.c \
   yee.c \
   yee.h \
   zpool.c \
   zpool.h
endif
  
if USE_CUDA
gsvit_LDADD=$(GLIB_LIBS) $(GWYDDION_LIBS) $(CUDA_LIBS)
gsvit_CFLAGS=$(GLIB_CFLAGS) $(CUDA_CFLAGS) $(GWYDDION_CFLAGS) $(UCUDA_CFLAGS) $(GTK_CFLAGS)
else
gsvit_LDADD=$(GLIB_LIBS) $(GWYDDION_LIBS)
gsvit_CFLAGS=$(GLIB_CFLAGS) $(GWYDDION_CFLAGS) $(UCUDA_CFLAGS) $(GTK_CFLAGS) 
endif

xsvit_LDADD=$(gsvit_LDADD) 
xsvit_CFLAGS=$(gsvit_CFLAGS) 

.cu.o:
	$(NVCC) -o $@ -c $< $(NVCCFLAGS)
.cu.lo:
	$(top_srcdir)/cudalt.py $@ $(NVCC) -c $(NVCCFLAGS) $<

