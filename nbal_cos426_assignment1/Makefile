#
# Makefile for producing Assignment #0 output
#

EXE=src/imgpro

IMGS=	clean	\
output/Gamma/gamma1.jpg\
output/Gamma/gamma2.jpg\
\
output/Compositing/globos_composite.jpg\
\
output/Extract/globos_extract1.jpg\
output/Extract/globos_extract2.jpg\
output/Extract/globos_extract3.jpg\
\
output/Bilateral/globos_bilateral.jpg\
output/Bilateral/globos_bilateral2.jpg\
output/Bilateral/globos_bilateral3.jpg\
\
output/Scale/globos_scale.jpg\
output/Scale/globos_scaleblur.jpg\
output/Scale/globos_scalebilerp1.jpg\
output/Scale/globos_scale2.jpg\
output/Scale/globos_scaleblur2.jpg\
output/Scale/globos_scalebilerp2.jpg\
\
output/Saturation/globos_saturation1.jpg\
output/Saturation/globos_saturation2.jpg\
output/Saturation/globos_saturation3.jpg\
output/Saturation/globos_saturation4.jpg\
\
output/Sharpen/globos_sharpen.jpg\
\
output/EdgeDetect/globos_edgedetector.jpg\
\
output/Rotate/globos_rotate1.jpg\
output/Rotate/globos_rotate2.jpg\
output/Rotate/globos_rotate3.jpg\
output/Rotate/globos_rotate1-1.jpg\
output/Rotate/globos_rotate2-1.jpg\
output/Rotate/globos_rotate3-1.jpg\
\
output/Median/globos_median1.jpg\
output/Median/globos_median2.jpg\
output/Median/globos_median3.jpg\
\
output/Fun/globos_fun1-1.jpg\
output/Fun/globos_fun2-1.jpg\
output/Fun/globos_fun3-1.jpg\
\
output/Motion_Blur/fast1.jpg\
output/Motion_Blur/fast2.jpg\
output/Motion_Blur/fast3.jpg\
output/Motion_Blur/fast4.jpg\
\
output/Dithering/globos_floyd1.jpg\
output/Dithering/globos_floyd2.jpg\
output/Dithering/globos_floyd3.jpg\
output/Dithering/globos_floyd4.jpg\
output/Dithering/globos_floyd5.jpg\
\
output/Dithering/globos_ordereddither1.jpg\
output/Dithering/globos_ordereddither2.jpg\
output/Dithering/globos_ordereddither3.jpg\
output/Dithering/globos_ordereddither4.jpg\
output/Dithering/globos_ordereddither5.jpg\
\
output/Dithering/globos_randomdither1.jpg\
output/Dithering/globos_randomdither2.jpg\
output/Dithering/globos_randomdither3.jpg\
output/Dithering/globos_randomdither4.jpg\
output/Dithering/globos_randomdither5.jpg\
\
output/Dithering/globos_quantize1.jpg\
output/Dithering/globos_quantize2.jpg\
output/Dithering/globos_quantize3.jpg\
output/Dithering/globos_quantize4.jpg\
output/Dithering/globos_quantize5.jpg\
\
output/Black_and_White/globos_blackandwhite.jpg\
\
output/Contrast/globos_contrast1.jpg\
output/Contrast/globos_contrast2.jpg\
output/Contrast/globos_contrast3.jpg\
output/Contrast/globos_contrast4.jpg\
\
output/Morph/trans/1.jpg\
output/Morph/trans/2.jpg\
output/Morph/trans/3.jpg\
output/Morph/trans/4.jpg\
output/Morph/trans/5.jpg\
output/Morph/trans/6.jpg\
output/Morph/trans/7.jpg\
output/Morph/trans/8.jpg\
output/Morph/trans/9.jpg\
output/Morph/trans/99.jpg\
output/Morph/trans/999.jpg\
output/Morph/trans/9999.jpg\
output/Morph/trans/99999.jpg\
output/Morph/trans/999999.jpg\
output/Morph/trans/9999999.jpg\
\
output/Blur/globos_blur1.jpg\
output/Blur/globos_blur2.jpg\
output/Blur/globos_blur3.jpg\

all: $(IMGS)
output/Gamma/gamma1.jpg: input/globos_de_colores.jpg
	$(EXE) $^ $@ -gamma 2
output/Gamma/gamma2.jpg: input/globos_de_colores.jpg
	$(EXE) $^ $@ -gamma .1

output/Compositing/globos_composite.jpg: input/GUYS.jpg
	$(EXE) $^ $@ -composite input/GUYSMATTE.jpg input/GUY.jpg input/GUYMATTE.jpg 0

output/Extract/globos_extract1.jpg: input/globos_de_colores.jpg
	$(EXE) $^ $@ -extract 0
output/Extract/globos_extract2.jpg: input/globos_de_colores.jpg
	$(EXE) $^ $@ -extract 1
output/Extract/globos_extract3.jpg: input/globos_de_colores.jpg
	$(EXE) $^ $@ -extract 2

output/Bilateral/globos_bilateral.jpg: input/luca.jpg
	$(EXE) $^ $@ -bilateral 0 0
output/Bilateral/globos_bilateral2.jpg: input/luca.jpg
	$(EXE) $^ $@ -bilateral .7 .1
output/Bilateral/globos_bilateral3.jpg: input/luca.jpg
	$(EXE) $^ $@ -bilateral 3 .1

output/Scale/globos_scale.jpg: input/globos_de_colores.jpg
	$(EXE) $^ $@ -sampling 1 -scale .5 .5
output/Scale/globos_scaleblur.jpg: input/globos_de_colores.jpg
	$(EXE) $^ $@ -sampling 2 -scale .5 .5 
output/Scale/globos_scalebilerp1.jpg: input/globos_de_colores.jpg
	$(EXE) $^ $@ -sampling 3 -scale .5 .5
output/Scale/globos_scale2.jpg: input/globos_de_colores.jpg
	$(EXE) $^ $@ -sampling 1 -scale 1.2 1.2
output/Scale/globos_scaleblur2.jpg: input/globos_de_colores.jpg
	$(EXE) $^ $@ -sampling 2 -scale 1.2 1.2 
output/Scale/globos_scalebilerp2.jpg: input/globos_de_colores.jpg
	$(EXE) $^ $@ -sampling 3 -scale 1.2 1.2

output/Saturation/globos_saturation1.jpg: input/globos_de_colores.jpg
	$(EXE) $^ $@ -saturation -1
output/Saturation/globos_saturation2.jpg: input/globos_de_colores.jpg
	$(EXE) $^ $@ -saturation 0
output/Saturation/globos_saturation3.jpg: input/globos_de_colores.jpg
	$(EXE) $^ $@ -saturation 1
output/Saturation/globos_saturation4.jpg: input/globos_de_colores.jpg
	$(EXE) $^ $@ -saturation 2

output/Sharpen/globos_sharpen.jpg: input/globos_de_colores.jpg
	$(EXE) $^ $@ -sharpen
	
output/EdgeDetect/globos_edgedetector.jpg: input/globos_de_colores.jpg
	$(EXE) $^ $@ -edge

output/Rotate/globos_rotate1.jpg: input/globos_de_colores.jpg
	$(EXE) $^ $@ -sampling 1 -rotate 1
output/Rotate/globos_rotate2.jpg: input/globos_de_colores.jpg
	$(EXE) $^ $@ -sampling 2 -rotate 1
output/Rotate/globos_rotate3.jpg: input/globos_de_colores.jpg
	$(EXE) $^ $@ -sampling 3 -rotate 1
output/Rotate/globos_rotate1-1.jpg: input/globos_de_colores.jpg
	$(EXE) $^ $@ -sampling 1 -rotate 2
output/Rotate/globos_rotate2-1.jpg: input/globos_de_colores.jpg
	$(EXE) $^ $@ -sampling 2 -rotate 2
output/Rotate/globos_rotate3-1.jpg: input/globos_de_colores.jpg
	$(EXE) $^ $@ -sampling 3 -rotate 2

output/Blur/globos_blur1.jpg: input/luca.jpg
	$(EXE) $^ $@ -blur 0
output/Blur/globos_blur2.jpg: input/luca.jpg
	$(EXE) $^ $@ -blur 2
output/Blur/globos_blur3.jpg: input/luca.jpg
	$(EXE) $^ $@ -blur 8

output/Median/globos_median1.jpg: input/image7.jpg
	$(EXE) $^ $@ -median 0
output/Median/globos_median2.jpg: input/image7.jpg
	$(EXE) $^ $@ -median .2
output/Median/globos_median3.jpg: input/image7.jpg
	$(EXE) $^ $@ -median .9

output/Fun/globos_fun1-1.jpg: input/globos_de_colores.jpg
	$(EXE) $^ $@ -sampling 1 -fun
output/Fun/globos_fun2-1.jpg: input/globos_de_colores.jpg
	$(EXE) $^ $@ -sampling 2 -fun
output/Fun/globos_fun3-1.jpg: input/globos_de_colores.jpg
	$(EXE) $^ $@ -sampling 3 -fun
	

output/Motion_Blur/fast1.jpg: input/fast.jpg
	$(EXE) $^ $@ -motionblur 0
output/Motion_Blur/fast2.jpg: input/fast.jpg
	$(EXE) $^ $@ -motionblur 10
output/Motion_Blur/fast3.jpg: input/fast.jpg
	$(EXE) $^ $@ -motionblur 20
output/Motion_Blur/fast4.jpg: input/fast.jpg
	$(EXE) $^ $@ -motionblur 30

output/Dithering/globos_floyd1.jpg: input/globos_de_colores.jpg
	$(EXE) $^ $@ -dither 2 1
output/Dithering/globos_floyd2.jpg: input/globos_de_colores.jpg
	$(EXE) $^ $@ -dither 2 2
output/Dithering/globos_floyd3.jpg: input/globos_de_colores.jpg
	$(EXE) $^ $@ -dither 2 3
output/Dithering/globos_floyd4.jpg: input/globos_de_colores.jpg
	$(EXE) $^ $@ -dither 2 4
output/Dithering/globos_floyd5.jpg: input/globos_de_colores.jpg
	$(EXE) $^ $@ -dither 2 5

output/Dithering/globos_ordereddither1.jpg: input/globos_de_colores.jpg
	$(EXE) $^ $@ -dither 1 1
output/Dithering/globos_ordereddither2.jpg: input/globos_de_colores.jpg
	$(EXE) $^ $@ -dither 1 2
output/Dithering/globos_ordereddither3.jpg: input/globos_de_colores.jpg
	$(EXE) $^ $@ -dither 1 3
output/Dithering/globos_ordereddither4.jpg: input/globos_de_colores.jpg
	$(EXE) $^ $@ -dither 1 4
output/Dithering/globos_ordereddither5.jpg: input/globos_de_colores.jpg
	$(EXE) $^ $@ -dither 1 5

output/Dithering/globos_randomdither1.jpg: input/globos_de_colores.jpg
	$(EXE) $^ $@ -dither 0 1
output/Dithering/globos_randomdither2.jpg: input/globos_de_colores.jpg
	$(EXE) $^ $@ -dither 0 2
output/Dithering/globos_randomdither3.jpg: input/globos_de_colores.jpg
	$(EXE) $^ $@ -dither 0 3
output/Dithering/globos_randomdither4.jpg: input/globos_de_colores.jpg
	$(EXE) $^ $@ -dither 0 4
output/Dithering/globos_randomdither5.jpg: input/globos_de_colores.jpg
	$(EXE) $^ $@ -dither 0 5

output/Dithering/globos_quantize1.jpg: input/globos_de_colores.jpg
	$(EXE) $^ $@ -quantize 1
output/Dithering/globos_quantize2.jpg: input/globos_de_colores.jpg
	$(EXE) $^ $@ -quantize 2
output/Dithering/globos_quantize3.jpg: input/globos_de_colores.jpg
	$(EXE) $^ $@ -quantize 3
output/Dithering/globos_quantize4.jpg: input/globos_de_colores.jpg
	$(EXE) $^ $@ -quantize 4
output/Dithering/globos_quantize5.jpg: input/globos_de_colores.jpg
	$(EXE) $^ $@ -quantize 5

output/Black_and_White/globos_blackandwhite.jpg: input/globos_de_colores.jpg
	$(EXE) $^ $@ -blackandwhite

output/Contrast/globos_contrast1.jpg: input/globos_de_colores.jpg
	$(EXE) $^ $@ -contrast -.5
output/Contrast/globos_contrast2.jpg: input/globos_de_colores.jpg
	$(EXE) $^ $@ -contrast 0
output/Contrast/globos_contrast3.jpg: input/globos_de_colores.jpg
	$(EXE) $^ $@ -contrast 1
output/Contrast/globos_contrast4.jpg: input/globos_de_colores.jpg
	$(EXE) $^ $@ -contrast 3

output/Morph/trans/1.jpg: input/fishA.jpg
	$(EXE) $^ $@ -morph input/fishB.jpg input/input.txt 0
output/Morph/trans/2.jpg: input/fishA.jpg
	$(EXE) $^ $@ -morph input/fishB.jpg input/input.txt .1
output/Morph/trans/3.jpg: input/fishA.jpg
	$(EXE) $^ $@ -morph input/fishB.jpg input/input.txt .15
output/Morph/trans/4.jpg: input/fishA.jpg
	$(EXE) $^ $@ -morph input/fishB.jpg input/input.txt .2
output/Morph/trans/5.jpg: input/fishA.jpg
	$(EXE) $^ $@ -morph input/fishB.jpg input/input.txt .25
output/Morph/trans/6.jpg: input/fishA.jpg
	$(EXE) $^ $@ -morph input/fishB.jpg input/input.txt .3
output/Morph/trans/7.jpg: input/fishA.jpg
	$(EXE) $^ $@ -morph input/fishB.jpg input/input.txt .35
output/Morph/trans/8.jpg: input/fishA.jpg
	$(EXE) $^ $@ -morph input/fishB.jpg input/input.txt .5
output/Morph/trans/9.jpg: input/fishA.jpg
	$(EXE) $^ $@ -morph input/fishB.jpg input/input.txt .55
output/Morph/trans/99.jpg: input/fishA.jpg
	$(EXE) $^ $@ -morph input/fishB.jpg input/input.txt .6
output/Morph/trans/999.jpg: input/fishA.jpg
	$(EXE) $^ $@ -morph input/fishB.jpg input/input.txt .65
output/Morph/trans/9999.jpg: input/fishA.jpg
	$(EXE) $^ $@ -morph input/fishB.jpg input/input.txt .75
output/Morph/trans/99999.jpg: input/fishA.jpg
	$(EXE) $^ $@ -morph input/fishB.jpg input/input.txt .85
output/Morph/trans/999999.jpg: input/fishA.jpg
	$(EXE) $^ $@ -morph input/fishB.jpg input/input.txt .95
output/Morph/trans/9999999.jpg: input/fishA.jpg
	$(EXE) $^ $@ -morph input/fishB.jpg input/input.txt 1
clean:
	rm -f $(IMGS)
