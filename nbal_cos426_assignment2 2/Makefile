#
# Makefile for producing Assignment #2 output
#

EXE=src/meshpro

IMGS=	art/TwilightZone.jpg		\
		output/split4.jpg			\
		output/inflate4.jpg		\
		output/rev1.jpg			\
		output/split1.jpg			\
		output/split2.jpg			\
		output/split3.jpg			\
		output/truncate1.jpg		\
		output/truncate2.jpg		\
		output/truncate3.jpg		\
		output/truncate4.jpg		\
		output/sharpen1.jpg		\
		output/sharpen2.jpg		\
		output/sharpen3.jpg		\
		output/smooth1.jpg		\
		output/smooth2.jpg		\
		output/smooth3.jpg		\
		output/inflate1.jpg		  \
		output/inflate2.jpg		  \
		output/inflate3.jpg 		  \
	   output/noise1.jpg			  \
	   output/noise2.jpg			  \
	   output/noise3.jpg			  \
		output/ico_rotate_0.5.jpg \
     	output/ico_scale_1_2_3.jpg \
		 output/fun1.jpg				\
		 		 output/fun2.jpg					

all: clean $(IMGS)

output/rev1.jpg: input/out.off $(EXE)
	$(EXE) input/out.off output/rev1.off -revolution input/out.off 0 0 0 0 0 1 .1 
	src/meshview output/rev1.off -output_image output/rev1.jpg -exit_immediately

output/split1.jpg: input/cube.off $(EXE)
	$(EXE) input/cube.off output/split1.off
	src/meshview output/split1.off -output_image output/split1.jpg -exit_immediately

output/split2.jpg: input/cube.off $(EXE)
	$(EXE) input/cube.off output/split2.off -loop
	src/meshview output/split2.off -output_image output/split2.jpg -exit_immediately
	
output/split3.jpg: input/cube.off $(EXE)
	$(EXE) input/cube.off output/split3.off -loop -loop
	src/meshview output/split3.off -output_image output/split3.jpg -exit_immediately
	
output/split4.jpg: input/cube.off $(EXE)
	$(EXE) input/cube.off output/split4.off -loop -loop -loop
	src/meshview output/split4.off -output_image output/split4.jpg -exit_immediately


output/truncate1.jpg: input/cube.off $(EXE)
	$(EXE) input/cube.off output/truncate1.off
	src/meshview output/truncate1.off -output_image output/truncate1.jpg -exit_immediately

output/truncate2.jpg: input/cube.off $(EXE)
	$(EXE) input/cube.off output/truncate2.off -truncate .5
	src/meshview output/truncate2.off -output_image output/truncate2.jpg -exit_immediately

output/truncate3.jpg: input/out.off $(EXE)
	$(EXE) input/out.off output/truncate3.off -revolution input/out.off 0 0 0 0 0 1 .1 -rotate 1.5 0 0 0 1 0 0
	src/meshview output/truncate3.off -output_image output/truncate3.jpg -exit_immediately

output/truncate4.jpg: input/out.off $(EXE)
	$(EXE) input/out.off output/truncate4.off -revolution input/out.off 0 0 0 0 0 1 .1 -rotate 1.5 0 0 0 1 0 0 -truncate .3
	src/meshview output/truncate4.off -output_image output/truncate4.jpg -exit_immediately

output/sharpen1.jpg: input/80.off $(EXE)
	$(EXE) input/80.off output/sharpen1.off
	src/meshview output/sharpen1.off -output_image output/sharpen1.jpg -exit_immediately

output/sharpen2.jpg: input/80.off $(EXE)
	$(EXE) input/80.off output/sharpen2.off -sharpen -sharpen
	src/meshview output/sharpen2.off -output_image output/sharpen2.jpg -exit_immediately

output/sharpen3.jpg: input/80.off $(EXE)
	$(EXE) input/80.off output/sharpen3.off -sharpen -sharpen -sharpen -sharpen
	src/meshview output/sharpen3.off -output_image output/sharpen3.jpg -exit_immediately

output/smooth1.jpg: input/cheetah.off $(EXE)
	$(EXE) input/cheetah.off output/smooth1.off -rotate  -1.5  0 0 0  0 1 0 
	src/meshview output/smooth1.off -output_image output/smooth1.jpg -exit_immediately

output/smooth2.jpg: input/cheetah.off $(EXE)
	$(EXE) input/cheetah.off output/smooth2.off -rotate  -1.5  0 0 0  0 1 0 -smooth -smooth
	src/meshview output/smooth2.off -output_image output/smooth2.jpg -exit_immediately

output/smooth3.jpg: input/cheetah.off $(EXE)
	$(EXE) input/cheetah.off output/smooth3.off -rotate  -1.5  0 0 0  0 1 0 -smooth -smooth -smooth -smooth -smooth
	src/meshview output/smooth3.off -output_image output/smooth3.jpg -exit_immediately

output/fun1.jpg: input/80.off $(EXE)
	$(EXE) input/80.off output/fun1.off
	src/meshview output/fun1.off -output_image output/fun1.jpg -exit_immediately

output/fun2.jpg: input/80.off $(EXE)
	$(EXE) input/80.off output/fun2.off -fun -rotate  1.5  0 0 0  0 1 0
	src/meshview output/fun2.off -output_image output/fun2.jpg -exit_immediately

art/TwilightZone.jpg: input/80.off $(EXE)
	$(EXE) input/80.off art/TwilightZone.off -fun -rotate  1.5  0 0 0  0 1 0
	src/meshview art/TwilightZone.off -output_image art/TwilightZone.jpg -exit_immediately

output/inflate1.jpg: input/cheetah.off $(EXE)
	$(EXE) input/cheetah.off output/inflate1.off -inflate .5 -rotate  -1.5  0 0 0  0 1 0
	src/meshview output/inflate1.off -output_image output/inflate1.jpg -exit_immediately

output/inflate2.jpg: input/cheetah.off $(EXE)
	$(EXE) input/cheetah.off output/inflate2.off -inflate 1 -rotate  -1.5  0 0 0  0 1 0
	src/meshview output/inflate2.off -output_image output/inflate2.jpg -exit_immediately

output/inflate3.jpg: input/cheetah.off $(EXE)
	$(EXE) input/cheetah.off output/inflate3.off -inflate 3 -rotate  -1.5  0 0 0  0 1 0
	src/meshview output/inflate3.off -output_image output/inflate3.jpg -exit_immediately

output/inflate4.jpg: input/cheetah.off $(EXE)
	$(EXE) input/cheetah.off output/inflate4.off -inflate 0 -rotate  -1.5  0 0 0  0 1 0
	src/meshview output/inflate4.off -output_image output/inflate4.jpg -exit_immediately

output/noise1.jpg: input/hand.off $(EXE)
	$(EXE) input/hand.off output/noise1.off -noise 0
	src/meshview output/noise1.off -output_image output/noise1.jpg -exit_immediately

output/noise2.jpg: input/hand.off $(EXE)
	$(EXE) input/hand.off output/noise2.off -noise .4
	src/meshview output/noise2.off -output_image output/noise2.jpg -exit_immediately

output/noise3.jpg: input/hand.off $(EXE)
	$(EXE) input/hand.off output/noise3.off -noise .8
	src/meshview output/noise3.off -output_image output/noise3.jpg -exit_immediately

output/ico_scale_1_2_3.jpg: input/ico.off $(EXE)
	$(EXE) input/ico.off output/ico_scale_1_2_3.off -scale 1 2 3
	src/meshview output/ico_scale_1_2_3.off -output_image output/ico_scale_1_2_3.jpg -exit_immediately

output/ico_scale_1_2_3.jpg: input/ico.off $(EXE)
	$(EXE) input/ico.off output/ico_scale_1_2_3.off -scale 1 2 3
	src/meshview output/ico_scale_1_2_3.off -output_image output/ico_scale_1_2_3.jpg -exit_immediately



output/ico_rotate_0.5.jpg: input/ico.off $(EXE)
	$(EXE) input/ico.off output/ico_rotate_0.5.off -rotate  0.5  0 0 0  0 0 1
	src/meshview output/ico_rotate_0.5.off -output_image output/ico_rotate_0.5.jpg -exit_immediately

output/ico_scale_1_2_3.jpg: input/ico.off $(EXE)
	$(EXE) input/ico.off output/ico_scale_1_2_3.off -scale 1 2 3
	src/meshview output/ico_scale_1_2_3.off -output_image output/ico_scale_1_2_3.jpg -exit_immediately

clean:
	rm -f $(IMGS)
	rm -f $(IMGS:.jpg=.off)
