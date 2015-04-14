#
# Makefile for producing Assignment #0 output
#

EXE=src/imgpro

IMGS=output/globos_brightness_0.5.jpg \
     output/globos_brightness_1.0.jpg \
     output/globos_brightness_1.5.jpg \
     output/islas_del_rosario_0.5.jpg \
     output/islas_del_rosario_1.0.jpg \
     output/islas_del_rosario_1.5.jpg \
     output/praia_da_sereia_0.5.jpg \
     output/praia_da_sereia_1.0.jpg \
     output/praia_da_sereia_1.5.jpg \
     output/praia_do_itapua_0.5.jpg \
     output/praia_do_itapua_1.0.jpg \
     output/praia_do_itapua_1.5.jpg 

all: $(IMGS)

output/globos_brightness_0.5.jpg: input/globos_de_colores.jpg
	$(EXE) $^ $@ -brightness 0.5

output/globos_brightness_1.0.jpg: input/globos_de_colores.jpg
	$(EXE) $^ $@ -brightness 1.0

output/globos_brightness_1.5.jpg: input/globos_de_colores.jpg
	$(EXE) $^ $@ -brightness 1.5

output/islas_del_rosario_0.5.jpg: input/islas_del_rosario.jpg
	$(EXE) $^ $@ -brightness 0.5

output/islas_del_rosario_1.0.jpg: input/islas_del_rosario.jpg
	$(EXE) $^ $@ -brightness 1.0

output/islas_del_rosario_1.5.jpg: input/islas_del_rosario.jpg
	$(EXE) $^ $@ -brightness 1.5

output/praia_da_sereia_0.5.jpg: input/praia_da_sereia.jpg
	$(EXE) $^ $@ -brightness 0.5

output/praia_da_sereia_1.0.jpg: input/praia_da_sereia.jpg
	$(EXE) $^ $@ -brightness 1.0

output/praia_da_sereia_1.5.jpg: input/praia_da_sereia.jpg
	$(EXE) $^ $@ -brightness 1.5

output/praia_do_itapua_0.5.jpg: input/praia_do_itapua.jpg
	$(EXE) $^ $@ -brightness 0.5

output/praia_do_itapua_1.0.jpg: input/praia_do_itapua.jpg
	$(EXE) $^ $@ -brightness 1.0

output/praia_do_itapua_1.5.jpg: input/praia_do_itapua.jpg
	$(EXE) $^ $@ -brightness 1.5

clean:
	rm -f $(IMGS)
