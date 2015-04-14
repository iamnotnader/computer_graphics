 // Source file for image class



// Include files 

#include "R2/R2.h"
#include "R2Pixel.h"
#include "R2Image.h"
#include <cstdlib>
#include <iostream>




////////////////////////////////////////////////////////////////////////
// Constructors/Destructors
////////////////////////////////////////////////////////////////////////
   //static double gaussian(int xsrc, int ysrc, int xdest, int ydest, int width, R2Image& from_image);

   R2Image::
   R2Image(void)
   : pixels(NULL),
    npixels(0),
    width(0), 
    height(0)
   {
   }



   R2Image::
   R2Image(const char *filename)
   : pixels(NULL),
    npixels(0),
    width(0), 
    height(0)
   {
   // Read image
      Read(filename);
   }



   R2Image::
   R2Image(int width, int height)
   : pixels(NULL),
    npixels(width * height),
    width(width), 
    height(height)
   {
   // Allocate pixels
      pixels = new R2Pixel [ npixels ];
      assert(pixels);
   }



   R2Image::
   R2Image(int width, int height, const R2Pixel *p)
   : pixels(NULL),
    npixels(width * height),
    width(width), 
    height(height)
   {
   // Allocate pixels
      pixels = new R2Pixel [ npixels ];
      assert(pixels);
   
   // Copy pixels 
      for (int i = 0; i < npixels; i++) 
         pixels[i] = p[i];
   }



   R2Image::
   R2Image(const R2Image& image)
   : pixels(NULL),
    npixels(image.npixels),
    width(image.width), 
    height(image.height)
    
   {
   // Allocate pixels
      pixels = new R2Pixel [ npixels ];
      assert(pixels);
   
   // Copy pixels 
      for (int i = 0; i < npixels; i++) 
         pixels[i] = image.pixels[i];
   }



   R2Image::
   ~R2Image(void)
   {
   // Free image pixels
      if (pixels) delete [] pixels;
   }



   R2Image& R2Image::
   operator=(const R2Image& image)
   {
   // Delete previous pixels
      if (pixels) { delete [] pixels; pixels = NULL; }
   
   // Reset width and height
      npixels = image.npixels;
      width = image.width;
      height = image.height;
   
   // Allocate new pixels
      pixels = new R2Pixel [ npixels ];
      assert(pixels);
   
   // Copy pixels 
      for (int i = 0; i < npixels; i++) 
         pixels[i] = image.pixels[i];
   
   // Return image
      return *this;
   }



////////////////////////////////////////////////////////////////////////
// Image processing functions
// YOU IMPLEMENT THE FUNCTIONS IN THIS SECTION
////////////////////////////////////////////////////////////////////////

// Per-pixel Operations ////////////////////////////////////////////////

   void R2Image::
   Brighten(double factor)
   {
   // Brighten the image by multiplying each pixel component by the factor,
   // then clamping the result to a valid range.
   
   // MAY ADD CODE HERE FROM ASSIGNMENT 0. NO CREDIT FOR THIS FEATURE
      for (int r = 0; r < Width(); r++)
      {
         for (int c = 0; c < Height(); c++)
         {
            Pixels(r)[c].SetRed(Pixels(r)[c].Red()*factor);
            Pixels(r)[c].SetBlue(Pixels(r)[c].Blue()*factor);
            Pixels(r)[c].SetGreen(Pixels(r)[c].Green()*factor);
            Pixels(r)[c].Clamp();
         }
      }
   }

   void R2Image::
   AddNoise(double factor)
   {
   // Add noise to an image.  The amount of noise is given by the factor
   // in the range [0.0..1.0].  0.0 adds no noise.  1.0 adds a lot of noise.
   
   // MAY FILL IN IMPLEMENTATION HERE (REMOVE PRINT STATEMENT WHEN DONE)
      fprintf(stderr, "AddNoise(%g) not implemented\n", factor);
   }

   void R2Image::
   ChangeContrast(double factor)
   {
   // Change the contrast of an image by interpolating between the image
   // and a constant gray image with the average luminance.
   // Interpolation reduces constrast, extrapolation boosts constrast,
   // and negative factors generate inverted images.
   
   	//Apply Gamma correction
      ApplyGamma(2.2);
   
   	//Calculate the average luminance
      double ltotal = 0;
      for (int i = 0; i < Width(); i++)
         for (int j = 0; j < Height(); j++)
            ltotal += Pixels(i)[j].Luminance();
      ltotal /= Height()*Width();
   	
   	//Interpolate between average luminance for every pixel
      for (int i = 0; i < Width(); i++)
         for (int j = 0; j < Height(); j++)
         {	
            Pixels(i)[j].SetRed((1 - factor)*ltotal + factor*Pixels(i)[j].Red());
            Pixels(i)[j].SetBlue((1 - factor)*ltotal + factor*Pixels(i)[j].Blue());
            Pixels(i)[j].SetGreen((1 - factor)*ltotal + factor*Pixels(i)[j].Green());
            Pixels(i)[j].Clamp();
         }
      
   	//Apply Gamma correction
      ApplyGamma(1/2.2);
   }

   void R2Image::
   ChangeSaturation(double factor)
   {
   // Changes the saturation of an image by interpolating between the
   // image and a gray level version of the image.  Interpolation
   // decreases saturation, extrapolation increases it, negative factors
   // preserve luminance  but invert the hue of the input image.
   
   	//Apply Gamma correction
      ApplyGamma(2.2);
   	
   	//Interpolate between pixel luminance and color value.
      for (int i = 0; i < Width(); i++)
         for (int j = 0; j < Height(); j++)
         {	
            Pixels(i)[j].SetRed((1 - factor)*Pixels(i)[j].Luminance() + factor*Pixels(i)[j].Red());
            Pixels(i)[j].SetBlue((1 - factor)*Pixels(i)[j].Luminance() + factor*Pixels(i)[j].Blue());
            Pixels(i)[j].SetGreen((1 - factor)*Pixels(i)[j].Luminance() + factor*Pixels(i)[j].Green());
            Pixels(i)[j].Clamp();
         }
         //Apply Gamma correction
      ApplyGamma(1/2.2);
   }

   void R2Image::
   ApplyGamma(double exponent)
   {
   // Apply a gamma correction with exponent to each pixel
   	
   	//Exponentiate every pixel.
      for (int r = 0; r < Width(); r++)
      {
         for (int c = 0; c < Height(); c++)
         {
            Pixels(r)[c].SetRed(pow(Pixels(r)[c].Red(), exponent));
            Pixels(r)[c].SetBlue(pow(Pixels(r)[c].Blue(), exponent));
            Pixels(r)[c].SetGreen(pow(Pixels(r)[c].Green(), exponent));
            Pixels(r)[c].Clamp();
         }
      }
   }

   void R2Image::
   BlackAndWhite(void)
   {
   // Replace each pixel with its luminance value
   // Put this in each channel,  so the result is grayscale
   
   	//Luminance declaration
      double luminance = 0;
      
   	//Apply Gamma correction
      ApplyGamma(2.2);
   	
   	//Replace each pixel with its luminance.
      for (int r = 0; r < Width(); r++)
      {
         for (int c = 0; c < Height(); c++)
         {
            luminance = Pixels(r)[c].Luminance();
            Pixels(r)[c].SetRed(luminance);
            Pixels(r)[c].SetBlue(luminance);
            Pixels(r)[c].SetGreen(luminance);
            Pixels(r)[c].Clamp();
         }
      }
      
      //Apply Gamma correction
      ApplyGamma(1/2.2);
   }

   void R2Image::
   ExtractChannel(int channel)
   {
   // Extracts a channel of an image (e.g., R2_IMAGE_RED_CHANNEL).  
   // Leaves the specified channel intact, 
   // and sets all the other ones to zero.
   
   // FILL IN IMPLEMENTATION HERE (REMOVE PRINT STATEMENT WHEN DONE)
   	
   	//Extract channel on every pixel.
      for (int r = 0; r < Width(); r++)
      {
         for (int c = 0; c < Height(); c++)
         {
            Pixels(r)[c][channel] = 0;
         }
      }   
   }

// Linear filtering ////////////////////////////////////////////////


   void R2Image::
   Blur(double sigma)
   {
   // Blur an image with a Gaussian filter with a given sigma.
   
      //Size of filter kernel.
      int size = ceil(3*sigma)*2+1;
   	//Radius of filter kernel.
      int radius = (size - 1)/2;
      
   	//Temporary image.
      R2Image* tempimage = new R2Image(Width(),Height());//-2*radius, Height()-2*radius);
   	
   	//Declarations.
      int dist = 0;
      int xpos = 0;
      int ypos = 0;
      double sum = 0;
      double redsum = 0;
      double bluesum = 0;
      double greensum = 0;
      double kernel[size][size];
   	
   	//Initialize Gaussian kernel.
      for (int i = 0; i < size; i++)
         for (int j = 0; j < size; j++)
            kernel[i][j] = 0;
   	
      for (int i = 0; i < size; i++)
      {
         for (int j = 0; j < size; j++)
         {
            dist = abs(radius-i);
            kernel[i][j] += exp((-.5)*pow(dist,2.0)/(2*pow((sigma),2)));
            dist = abs(radius-j);
            kernel[i][j] += exp((-.5)*pow(dist,2.0)/(2*pow((sigma),2)));
            sum += kernel[i][j];
         }
      }
      
      //Normalize filter kernel.
      for (int i = 0; i < size; i++)
      {
         for (int j = 0; j < size; j++)
         {
            kernel[i][j] = kernel[i][j] / sum;
         }
      }
   	
   	//Edge case if size is 1.
      if (size == 1)
         kernel[0][0] = 1;
   		
   	//Adjust for Gamma correction.
      ApplyGamma(2.2);	
   	
   	//Convolute on every pixel
      for (int r = 0; r < Width(); r++)
      {
         for (int c = 0; c < Height(); c++)
         {				
            redsum = bluesum = greensum = 0;
            for (int i = 0; i < size; i++)
            {
            	//Compute x position and y position of surrounding pixels;
            	//Use a mirror of the pixels when pixel is off the edge. 
               ypos = c + (i - radius);
               if (ypos < 0)
                  ypos = (-1)*ypos;
               else if (ypos >= (Height()))
                  ypos = Height() - (ypos - Height()) - 1;
               for (int j = 0; j < size; j++)
               {
                  xpos = r + (j - radius);
                  if (xpos < 0)
                     xpos = (-1)*xpos;
                  else if (xpos >= (Width()))
                     xpos = Width() - (xpos - Width()) - 1;
               		
                  redsum += Pixels(xpos)[ypos].Red()*kernel[i][j]; 
                  bluesum += Pixels(xpos)[ypos].Blue()*kernel[i][j]; 
                  greensum += Pixels(xpos)[ypos].Green()*kernel[i][j]; 
               }
            }
         	//Set the values of our temporary image.
            tempimage->Pixels(r)[c].SetRed(redsum);
            tempimage->Pixels(r)[c].SetBlue(bluesum);
            tempimage->Pixels(r)[c].SetGreen(greensum);
         }
      }	
      
   	//Make this image point to our temporary image.
      delete this->pixels;	
      this -> pixels = tempimage->pixels;
      this->npixels = tempimage->npixels;
      this->height = tempimage->height;
      this->width = tempimage->width;
   	
   	//Readjust for Gamma correction
      ApplyGamma(1/2.2);
   }

   void R2Image::
   Sharpen()
   {
   // Sharpen an image using a linear filter
   
   	//Radius of filter kernel.
      int size = 3;
      int radius = (size - 1)/2;
      
   	//Temporary image.
      R2Image* tempimage = new R2Image(Width(), Height());
   	
   	//Declarations.
      int xpos = 0;
      int ypos = 0;
      double redsum = 0;
      double bluesum = 0;
      double greensum = 0;
      double kernel[3][3] = {{0, -2/3., 0},{-2/3.,11/3.,-2/3.},{0,-2/3.,0}};
   		
   	//Adjust for Gamma correction.
      ApplyGamma(2.2);	
   	
   	//Convolute on every pixel
      for (int r = 0; r < (Width()); r++)
      {
         for (int c = 0; c < (Height()); c++)
         {				
         	//If out of filter bounds, just leave the pixel alone.
            if ((r < radius) || (r > Width() - radius - 1) 
            || (c < radius) || (c > Height() - radius - 1))
            {
               tempimage->Pixels(r)[c].SetRed(Pixels(r)[c].Red());
               tempimage->Pixels(r)[c].SetBlue(Pixels(r)[c].Blue());
               tempimage->Pixels(r)[c].SetGreen(Pixels(r)[c].Green());
               continue;
            }
            
         	//Otherwise convolute.
            redsum = bluesum = greensum = 0;
            for (int i = 0; i < size; i++)
            {
               ypos = c + (i - radius);
               for (int j = 0; j < size; j++)
               {
                  xpos = r + (j - radius);
                  redsum += Pixels(xpos)[ypos].Red()*kernel[i][j]; 
                  bluesum += Pixels(xpos)[ypos].Blue()*kernel[i][j]; 
                  greensum += Pixels(xpos)[ypos].Green()*kernel[i][j]; 
               }
            }
            
         	//Set values in the temporary image.
            tempimage->Pixels(r)[c].SetRed(redsum);
            tempimage->Pixels(r)[c].SetBlue(bluesum);
            tempimage->Pixels(r)[c].SetGreen(greensum);
         }
      }	
      
   	//Make this image point to our temporary image.
      delete this->pixels;	
      this -> pixels = tempimage->pixels;
      this->npixels = tempimage->npixels;
      this->height = tempimage->height;
      this->width = tempimage->width;
   	
   	//Readjust for Gamma correction
      ApplyGamma(1/2.2);
   }


   void R2Image::
   EdgeDetect(void)
   {
   // Detect edges in an image.
   
      //Size and radius of filter kernel.
      int size = 3;
      int radius = (size - 1)/2;
      
   	//Temporary image.
      R2Image* tempimage = new R2Image(Width(), Height());
   	
   	//Declarations.
      int xpos = 0;
      int ypos = 0;
      double redsum = 0;
      double bluesum = 0;
      double greensum = 0;
      double kernel[3][3] = {{-1,-1,-1},{-1,8,-1},{-1,-1,-1}};
   		
   	//Adjust for Gamma correction.
      ApplyGamma(2.2);	
   	
   	//Convolute on every pixel
      for (int r = 0; r < (Width()); r++)
      {
         for (int c = 0; c < (Height()); c++)
         {				
         	//If out of filter bounds, just leave the pixel alone.
            if ((r < radius) || (r > Width() - radius - 1) 
            || (c < radius) || (c > Height() - radius - 1))
            {
               tempimage->Pixels(r)[c].SetRed(Pixels(r)[c].Red());
               tempimage->Pixels(r)[c].SetBlue(Pixels(r)[c].Blue());
               tempimage->Pixels(r)[c].SetGreen(Pixels(r)[c].Green());
               continue;
            }
            
         	//Otherwise, convolute using filter.
            redsum = bluesum = greensum = 0;
            for (int i = 0; i < size; i++)
            {
               ypos = c + (i - radius);
               for (int j = 0; j < size; j++)
               {
                  xpos = r + (j - radius);
                  redsum += Pixels(xpos)[ypos].Red()*kernel[i][j]; 
                  bluesum += Pixels(xpos)[ypos].Blue()*kernel[i][j]; 
                  greensum += Pixels(xpos)[ypos].Green()*kernel[i][j]; 
               }
            }
            
         	//Set values in the temporary image.
            tempimage->Pixels(r)[c].SetRed(redsum);
            tempimage->Pixels(r)[c].SetBlue(bluesum);
            tempimage->Pixels(r)[c].SetGreen(greensum);
         }
      }	
      
   	//Make this image point to our temporary image.
      delete this->pixels;	
      this -> pixels = tempimage->pixels;
      this->npixels = tempimage->npixels;
      this->height = tempimage->height;
      this->width = tempimage->width;
   	
   	//Readjust for Gamma correction
      ApplyGamma(1/2.2);
   }



   void R2Image::
   MotionBlur(int amount)
   {
   // Perform horizontal motion blur
   
   // convolve in X direction with a linear ramp of amount non-zero pixels
   // the image should be strongest on the right hand side (see example)
   
   	//If amount is zero, don't blur.
      if (amount == 0)
         return;
   
   	//Declarations.
      double sum, redsum, bluesum, greensum;
      double filter[amount];
   	
   	//Temporary image.
      R2Image* tempimage = new R2Image(Width(), Height());
   	
   	//Initialize the 1D filter.
      for (int k = 0; k < amount; k++)
         filter[k] = amount - k;
   	
   	//Adjust for Gamma correction.
      ApplyGamma(2.2);	
   	
   	//Convolute on every pixel
      for (int i = 0; i < Width(); i++)
      {
         for (int j = 0; j < Height(); j++)
         {
            redsum = bluesum = greensum = sum = 0;
            for (int k = 0; k < amount; k++)
            {	
            	//Apply filter to pixels inside bounds.
               if ((i + k) < Width())
               {
                  redsum += Pixels(i+k)[j].Red()*filter[k]; 
                  bluesum += Pixels(i+k)[j].Blue()*filter[k]; 
                  greensum += Pixels(i+k)[j].Green()*filter[k];
                  sum += filter[k];
               }
            }
            
         	//Set the pixel values in the temporary image.
            tempimage->Pixels(i)[j].SetRed(redsum/sum);
            tempimage->Pixels(i)[j].SetBlue(bluesum/sum);
            tempimage->Pixels(i)[j].SetGreen(greensum/sum);
         }
      }
   	
   	//Set this image to point to our temporary image.
      delete this->pixels;	
      this -> pixels = tempimage->pixels;
      this->npixels = tempimage->npixels;
      this->height = tempimage->height;
      this->width = tempimage->width;
      
   	//Adjust for Gamma correction.
      ApplyGamma(1/2.2);	
   }


// Non-Linear filtering ////////////////////////////////////////////////

   void R2Image::
   MedianFilter(double sigma)
   {
   // Perform median filtering with a given width
   
   	//Size and radius of filter kernel.	
      int size = ceil(3*sigma)*2+1;
      int radius = (size - 1)/2;
      
   	//Declarations.
      int xpos, ypos;
      int min;
      R2Pixel* temp; 	
      R2Pixel* medianlum[size*size];
   						
   	//Temporary image.
      R2Image* tempimage = new R2Image(Width(), Height());
   			
   	//Adjust for Gamma correction.
      ApplyGamma(2.2);		
   	
   	//Convolute on every pixel.
      for (int r = 0; r < (Width()); r++)
      {
         for (int c = 0; c < (Height()); c++)
         {
            for (int i = 0; i < size; i++)
            {
            	//Compute x position and y position of surrounding pixels;
            	//Use a mirror of the pixels when pixel is off the edge. 
               ypos = c + (i - radius);
               if (ypos < 0)
                  ypos = (-1)*ypos;
               else if (ypos >= (Height()))
                  ypos = Height() - (ypos - Height()) - 1;
               for (int j = 0; j < size; j++)
               {
                  xpos = r + (j - radius);
                  if (xpos < 0)
                     xpos = (-1)*xpos;
                  else if (xpos >= (Width()))
                     xpos = Width() - (xpos - Width()) - 1;
               	
               	//Store pointers to pixels in an array.
                  medianlum[i*size + j] = &(Pixels(xpos)[ypos]);
               }
            }
         	
         	//Sort first size/2 elements of the pixel array
            for (int i = 0; i < size*size/2+1; i++)
            {
               min = i;
               for (int j = i; j < size*size; j++)
               {
                  if ((*medianlum[j]).Luminance() < (*medianlum[min]).Luminance())
                     min = j;
               }
               temp = medianlum[i];
               medianlum[i] = medianlum[min];
               medianlum[min] = temp;										
            }
            
         	//Pixel at index size/2 is the median.
            tempimage->Pixels(r)[c].Reset(medianlum[size/2]->Red(), 
               medianlum[size/2]->Blue(), medianlum[size/2]->Green(), 0);
         }
      }
      
      //Make this image point to our temporary image.
      delete this->pixels;	
      this -> pixels = tempimage->pixels;
      this->npixels = tempimage->npixels;
      this->height = tempimage->height;
      this->width = tempimage->width;
   
   	//Reapply gamma correction.		
      ApplyGamma(1/2.2);
   }

   void R2Image::
   BilateralFilter(double rangesigma, double domainsigma)
   {
   // Perform bilater val filtering with a given range and domain width
     
     	//Edge case when sigmas are zero.
      if (rangesigma == 0 || domainsigma == 0)
         return;
   
   	//Size of filter kernel.	  
      int size = ceil(3*domainsigma)*2+1;		
   	//Radius of filter kernel.
      int radius = (size - 1)/2;
      
   	//Temporary image.
      R2Image* tempimage = new R2Image(Width(), Height());
   	
   	//Declarations.
      double dist = 0;
      int xpos = 0;
      int ypos = 0;
      double sum = 0;
      double diff = 0;
      double redsum = 0;
      double bluesum = 0;
      double greensum = 0;
      double kernel[size][size];
   	
   	//Initialize Gaussian kernel.
      for (int i = 0; i < size; i++)
         for (int j = 0; j < size; j++)
            kernel[i][j] = 0;
      for (int i = 0; i < size; i++)
      {
         for (int j = 0; j < size; j++)
         {
            dist = abs(radius-i);
            kernel[i][j] += exp((-.5)*pow(dist,2.0)/(2*pow((domainsigma),2)));
            dist = abs(radius-j);
            kernel[i][j] += exp((-.5)*pow(dist,2.0)/(2*pow((domainsigma),2)));
            sum += kernel[i][j];
         }
      }
   	
   	//Edge case if size is 1.
      if (size == 1)
         kernel[0][0] = 1;
   		
   	//Adjust for Gamma correction.
      ApplyGamma(2.2);	
   	
   	//Convolute on every pixel
      for (int r = 0; r < (Width()); r++)
      {
         for (int c = 0; c < (Height()); c++)
         {				
            redsum = bluesum = greensum = 0;
            sum = 0;
            for (int i = 0; i < size; i++)
            {
               //Compute x position and y position of surrounding pixels;
            	//Use a mirror of the pixels when pixel is off the edge. 
               ypos = c + (i - radius);
               if (ypos < 0)
                  ypos = (-1)*ypos;
               else if (ypos >= (Height()))
                  ypos = Height() - (ypos - Height()) - 1;
               for (int j = 0; j < size; j++)
               {
                  xpos = r + (j - radius);
                  if (xpos < 0)
                     xpos = (-1)*xpos;
                  else if (xpos >= (Width()))
                     xpos = Width() - (xpos - Width()) - 1;
               	
               	//Compute color distance for range function.
                  dist = abs(Pixels(r)[c].Luminance() - Pixels(xpos)[ypos].Luminance());
                  diff = exp((-.5)*pow(dist,2.0)/(2*pow((rangesigma),2)));
               	
                  redsum += Pixels(xpos)[ypos].Red()*kernel[i][j]*diff; 
                  bluesum += Pixels(xpos)[ypos].Blue()*kernel[i][j]*diff; 
                  greensum += Pixels(xpos)[ypos].Green()*kernel[i][j]*diff;
               	
                  sum+=kernel[i][j]*diff;
               }
            }
            
            //Set the value of the temporary image.
            tempimage->Pixels(r)[c].SetRed(redsum/sum);
            tempimage->Pixels(r)[c].SetBlue(bluesum/sum);
            tempimage->Pixels(r)[c].SetGreen(greensum/sum);
         }
      }	
      
   	//Make this image point to our temporary image.
      delete this->pixels;	
      this -> pixels = tempimage->pixels;
      this->npixels = tempimage->npixels;
      this->height = tempimage->height;
      this->width = tempimage->width;
   	
   	//Readjust for Gamma correction
      ApplyGamma(1/2.2);
   }


// Resampling operations  ////////////////////////////////////////////////


   void R2Image::
   Scale(double sx, double sy, int sampling_method)
   {
   // Scale an image in x by sx, and y by sy.
   
      //Create a temporary image with scaled dimensions.
      int newxmax = (int)(Width()*sx);
      int newymax = (int)(Height()*sy);
      R2Image* tempimage = new R2Image(newxmax, newymax); 	
   	
   	//Adjust for Gamma correction
      ApplyGamma(2.2);
   	
   	//Use point sampling.
      if (sampling_method == 1)
      {  
         int ixpos, iypos;
      
      	//Sample every pixel in the new image.
         for (int i = 0; i < newxmax; i++)
            for (int j = 0; j < newymax; j++)
            {
               ixpos = i/sx;
               iypos = j/sy;
            	
            	//Set the values in the temporary image.	   
               tempimage->Pixels(i)[j].SetRed(Pixels(ixpos)[iypos].Red());
               tempimage->Pixels(i)[j].SetBlue(Pixels(ixpos)[iypos].Blue());
               tempimage->Pixels(i)[j].SetGreen(Pixels(ixpos)[iypos].Green());
            }
      }
      
      
      /****************************************************/
      /****************************************************/
      //Use Gaussian sampling.
      else if (sampling_method == 2)
      {
      	//Size and radius of filter.
         int size = max(1/sx, 1/sy);
         if (size%2 == 0)
            size++;
         if (size < 3)
            size = 3;
         int radius = (size - 1)/2;
         
       	//Declarations.  
         double redsum, bluesum, greensum, dist, ksum;
         int ixpos, iypos, xlo, xhi, ylo, yhi, filterposx, filterposy;
         double kernel[size][size];
      
      
      	//Initialize gaussian kernel w/ size
         for (int i = 0; i < size; i++)
            for (int j = 0; j < size; j++)
               kernel[i][j] = 0;
      	
         ksum = 0;
         for (int i = 0; i < size; i++)
         {
            for (int j = 0; j < size; j++)
            {
               dist = abs(radius-i);
               kernel[i][j] += exp((-.5)*pow(dist,2.0)/(2*pow((radius/3.),2)));
               dist = abs(radius-j);
               kernel[i][j] += exp((-.5)*pow(dist,2.0)/(2*pow((radius/3.),2)));
               ksum += kernel[i][j];
            }
         }
      
      	//Normalize filter kernel.
         for (int i = 0; i < size; i++)
         {
            for (int j = 0; j < size; j++)
            {
               kernel[i][j] = kernel[i][j] / ksum;
            }
         }
      
      	//Sample every every pixel.
         for (int i = 0; i < newxmax; i++)
            for (int j = 0; j < newymax; j++)
            {
            	//X and Y positions in the original image.
               ixpos = i/sx;
               iypos = j/sy;
            	
            	//Check bounds and adjust-- not important when scaling.
               if (ixpos >= Width() - radius)
                  ixpos = Width() - radius - 1;
               if (iypos >= Height() - radius)
                  iypos = Height() - radius - 1;
            	
            	//Apply Filter.     
               redsum = bluesum = greensum = ksum = 0;
               xlo = ixpos - radius;
               xhi = ixpos + radius;
               ylo = iypos - radius;
               yhi = iypos + radius;
               for (int ix = xlo; ix <= xhi; ix++) 
               {
                  for (int iy = ylo; iy <= yhi; iy++) 
                  {
                     filterposx = ix - xlo;
                     filterposy = iy - ylo;
                     
                     redsum += kernel[filterposx][filterposy]*Pixels(ix)[iy].Red();
                     bluesum += kernel[filterposx][filterposy]*Pixels(ix)[iy].Blue();
                     greensum += kernel[filterposx][filterposy]*Pixels(ix)[iy].Green();
                     ksum += kernel[filterposx][filterposy];
                  }
               }
               //Normalize color values.
               redsum /= ksum;
               bluesum /= ksum;
               greensum /= ksum;           
                  
               //Set values in the temporary image.
               tempimage->Pixels(i)[j].SetRed(redsum);
               tempimage->Pixels(i)[j].SetBlue(bluesum);
               tempimage->Pixels(i)[j].SetGreen(greensum);
            }
      }
      
      /**************************************************/
      /**************************************************/
      //Use Bilerp.
      else if (sampling_method == 3)
      {
      	//Declarations.
         double reda, bluea, greena, redb, blueb, greenb;
         double factorx, factory;
         double red, green, blue;
         double dxpos, dypos;
         int ixpos, iypos;
      	
      	//Sample every pixel.
         for (int i = 0; i < newxmax; i++)
            for (int j = 0; j < newymax; j++)
            {
            	//New pixel locations.
               dxpos = i/sx;
               dypos = j/sy;
               ixpos = (int)dxpos;
               iypos = (int)dypos;
            	
            	//Check bounds.
               if (ixpos >= Width() - 1)
                  ixpos = Width() - 2;
               if (iypos >= Height() - 1)	
                  iypos = Height() - 2;
            	
            	//Set factor to difference pixel positions.
               factorx = dxpos - ixpos;
               factory = dypos - iypos;					
               
            	//Bilerp each color.
               reda = factorx*Pixels(ixpos)[iypos].Red() 
                    + (1 - factorx)*Pixels(ixpos+1)[iypos].Red();
               redb = factorx*Pixels(ixpos)[iypos+1].Red() 
                    + (1 - factorx)*Pixels(ixpos+1)[iypos+1].Red();
               red = factory*reda + (1 - factory)*redb;
            	
            	
               bluea = factorx*Pixels(ixpos)[iypos].Blue() 
                     + (1 - factorx)*Pixels(ixpos+1)[iypos].Blue();
               blueb = factorx*Pixels(ixpos)[iypos+1].Blue() 
                     + (1 - factorx)*Pixels(ixpos+1)[iypos+1].Blue();
               blue = factory*bluea + (1 - factory)*blueb;
            	
            	
               greena = factorx*Pixels(ixpos)[iypos].Green() 
                      + (1 - factorx)*Pixels(ixpos+1)[iypos].Green();
               greenb = factorx*Pixels(ixpos)[iypos+1].Green() 
                  	 + (1 - factorx)*Pixels(ixpos+1)[iypos+1].Green();
               green = factory*greena + (1 - factory)*greenb;
            	
            	//Set the values in the temporary image.    
               tempimage->Pixels(i)[j].SetRed(red);
               tempimage->Pixels(i)[j].SetBlue(blue);
               tempimage->Pixels(i)[j].SetGreen(green);
            }
      }
      
   	//Make this image point to the temporary image.
      delete this->pixels;	
      this -> pixels = tempimage->pixels;
      this->npixels = tempimage->npixels;
      this->height = tempimage->height;
      this->width = tempimage->width;
      
   	//Adjust for Gamma correction
      ApplyGamma(1/2.2);
   }
	
   void R2Image::
   Rotate(double angle, int sampling_method)
   {
   // Rotate an image by the given angle.
      if (sampling_method == 1)
      {
         int ixpos, iypos;
         //Sine and Cosine of angle.
         double sine = sin(-1*angle);
         double cosine = cos(-1*angle);
           //Temporary image.
         R2Image* tempimage = new R2Image(Width(), Height());
       	//Apply Gamma correction.
         ApplyGamma(2.2);  
      	
         //Sample every pixel in the new image.
         for (int i = 0; i < Width(); i++)
            for (int j = 0; j < Height(); j++)
            {
               //Find new positions (shift to center)
               ixpos = ((i-Width()/2.)*cosine + (j-Height()/2.)*sine);
               iypos = ((j-Height()/2.)*cosine - (i-Width()/2.)*sine);							
               
               ixpos += Width()/2.;
               iypos += Height()/2.;
               					
               //Check bounds-- if out of bounds, make black.
               if ((ixpos >= Width() - 1)
                  || (iypos >= Height() - 1)
                  || (ixpos < 0)
                  || (iypos < 0))
               {
                  tempimage->Pixels(i)[j].SetRed(0);
                  tempimage->Pixels(i)[j].SetRed(0);
                  tempimage->Pixels(i)[j].SetRed(0);
                  continue;
               }
               
               //Set the values in the temporary image.	   
               tempimage->Pixels(i)[j].SetRed(Pixels(ixpos)[iypos].Red());
               tempimage->Pixels(i)[j].SetBlue(Pixels(ixpos)[iypos].Blue());
               tempimage->Pixels(i)[j].SetGreen(Pixels(ixpos)[iypos].Green());
            }
         
         //Make this image point to the temporary image.
         delete this->pixels;	
         this -> pixels = tempimage->pixels;
         this->npixels = tempimage->npixels;
         this->height = tempimage->height;
         this->width = tempimage->width;
         ApplyGamma(1/2.2);
      }
      
      
      
      
      else if (sampling_method == 2)
      {
      //Size and radius of Gaussian kernel (for sampling)
         int size = 5;	
         int radius = (size - 1)/2;
      
      //Sine and Cosine of angle.
         double sine = sin(-1*angle);
         double cosine = cos(-1*angle);
      
      //Temporary image.
         R2Image* tempimage = new R2Image(Width(), Height());
      
      //Declarations.
         double redsum, bluesum, greensum, dist, ksum;
         int ixpos, iypos, xlo, xhi, ylo, yhi, filterposx, filterposy;
         double kernel[size][size];
      
      //Apply Gamma correction.
         ApplyGamma(2.2);
      
      
      //Initialize gaussian kernel w/ size.
         for (int i = 0; i < size; i++)
            for (int j = 0; j < size; j++)
               kernel[i][j] = 0;
      	
         ksum = 0;
         for (int i = 0; i < size; i++)
         {
            for (int j = 0; j < size; j++)
            {
               dist = abs(radius-i);
               kernel[i][j] += exp((-.5)*pow(dist,2.0)/(2*pow((radius/3.),2)));
               dist = abs(radius-j);
               kernel[i][j] += exp((-.5)*pow(dist,2.0)/(2*pow((radius/3.),2)));
               ksum += kernel[i][j];
            }
         }
      
      //Normalize filter kernel.
         for (int i = 0; i < size; i++)
         {
            for (int j = 0; j < size; j++)
            {
               kernel[i][j] = kernel[i][j] / ksum;
            }
         }
      
      //Resample every pixel.
         for (int i = 0; i < Width(); i++)
            for (int j = 0; j < Height(); j++)
            {
            //Find new positions (shift to center)
               ixpos = ((i-Width()/2.)*cosine + (j-Height()/2.)*sine);
               iypos = ((j-Height()/2.)*cosine - (i-Width()/2.)*sine);							
            
               ixpos += Width()/2.;
               iypos += Height()/2.;
            						
            //Check bounds-- if out of bounds, make black.
               if ((ixpos >= Width() - 2*radius)
               || (iypos >= Height() - 2*radius)
               || (ixpos <= 2*radius)
               || (iypos <= 2*radius))
               {
                  tempimage->Pixels(i)[j].SetRed(0);
                  tempimage->Pixels(i)[j].SetRed(0);
                  tempimage->Pixels(i)[j].SetRed(0);
                  continue;
               }         	    		
                           
            //Apply Gaussian filter to resample.
               redsum = bluesum = greensum = ksum = 0;
               xlo = ixpos - radius;
               xhi = ixpos + radius;
               ylo = iypos - radius;
               yhi = iypos + radius;
               for (int ix = xlo; ix <= xhi; ix++) 
               {
                  for (int iy = ylo; iy <= yhi; iy++) 
                  {
                     filterposx = ix - xlo;
                     filterposy = iy - ylo;
                  
                     redsum += kernel[filterposx][filterposy]*Pixels(ix)[iy].Red();
                     bluesum += kernel[filterposx][filterposy]*Pixels(ix)[iy].Blue();
                     greensum += kernel[filterposx][filterposy]*Pixels(ix)[iy].Green();
                     ksum += kernel[filterposx][filterposy];
                  }
               }
            //Normalize color values.
               redsum /= ksum;
               bluesum /= ksum;
               greensum /= ksum;           
                  
            //Set values in the temporary image.  
               tempimage->Pixels(i)[j].SetRed(redsum);
               tempimage->Pixels(i)[j].SetBlue(bluesum);
               tempimage->Pixels(i)[j].SetGreen(greensum);
            }
       	//Make this image point to the temporary image.
         delete this->pixels;	
         this -> pixels = tempimage->pixels;
         this->npixels = tempimage->npixels;
         this->height = tempimage->height;
         this->width = tempimage->width;
         
         ApplyGamma(1/2.2);
      } 
      
      else if (sampling_method == 3)
      {
           //Temporary image.
         R2Image* tempimage = new R2Image(Width(), Height());
       	//Apply Gamma correction.
         ApplyGamma(2.2);  
      	
         //Declarations.
         double reda, bluea, greena, redb, blueb, greenb;
         double factorx, factory;
         double red, green, blue;
         double dxpos, dypos;
         int ixpos, iypos;
         //Sine and Cosine of angle.
         double sine = sin(-1*angle);
         double cosine = cos(-1*angle);
      	
      	//Sample every pixel.
         for (int i = 0; i < Width(); i++)
            for (int j = 0; j < Height(); j++)
            {
            	//Find new positions (shift to center)
               dxpos = ((i-Width()/2.)*cosine + (j-Height()/2.)*sine);
               dypos = ((j-Height()/2.)*cosine - (i-Width()/2.)*sine);							
               dxpos += Width()/2.;
               dypos += Height()/2.;
               
            	ixpos = (int)dxpos;
            	iypos = (int)dypos;
            	
            	//Check bounds.
            	if ((ixpos >= Width() - 1)
                  || (iypos >= Height() - 1)
                  || (ixpos < 0)
                  || (iypos < 0))
               {
                  tempimage->Pixels(i)[j].SetRed(0);
                  tempimage->Pixels(i)[j].SetRed(0);
                  tempimage->Pixels(i)[j].SetRed(0);
                  continue;
               }
            	
            	//Set factor to difference pixel positions.
               factorx = dxpos - ixpos;
               factory = dypos - iypos;					
               
            	//Bilerp each color.
               reda = factorx*Pixels(ixpos)[iypos].Red() 
                    + (1 - factorx)*Pixels(ixpos+1)[iypos].Red();
               redb = factorx*Pixels(ixpos)[iypos+1].Red() 
                    + (1 - factorx)*Pixels(ixpos+1)[iypos+1].Red();
               red = factory*reda + (1 - factory)*redb;
            	
            	
               bluea = factorx*Pixels(ixpos)[iypos].Blue() 
                     + (1 - factorx)*Pixels(ixpos+1)[iypos].Blue();
               blueb = factorx*Pixels(ixpos)[iypos+1].Blue() 
                     + (1 - factorx)*Pixels(ixpos+1)[iypos+1].Blue();
               blue = factory*bluea + (1 - factory)*blueb;
            	
            	
               greena = factorx*Pixels(ixpos)[iypos].Green() 
                      + (1 - factorx)*Pixels(ixpos+1)[iypos].Green();
               greenb = factorx*Pixels(ixpos)[iypos+1].Green() 
                  	 + (1 - factorx)*Pixels(ixpos+1)[iypos+1].Green();
               green = factory*greena + (1 - factory)*greenb;
            	
            	//Set the values in the temporary image.    
               tempimage->Pixels(i)[j].SetRed(red);
               tempimage->Pixels(i)[j].SetBlue(blue);
               tempimage->Pixels(i)[j].SetGreen(green);
            }
      
      	//Make this image point to the temporary image.
         delete this->pixels;	
         this -> pixels = tempimage->pixels;
         this->npixels = tempimage->npixels;
         this->height = tempimage->height;
         this->width = tempimage->width;
         ApplyGamma(1/2.2);
      }
   }


   void R2Image::
   Fun(int sampling_method)
   {
   // Warp an image using a creative filter of your choice.
   // Using swirl.
   
   	if (sampling_method == 1)
   	{
   		int ixpos, iypos;
   		   
      	//Declarations.
         double rotamount, rotsize, angle;	
      	
      	//Parameters.
         rotamount = 10;
         rotsize = 200;
      	
           //Temporary image.
         R2Image* tempimage = new R2Image(Width(), Height());
       	//Apply Gamma correction.
         ApplyGamma(2.2);  
      	
         //Sample every pixel in the new image.
         for (int i = 0; i < Width(); i++)
            for (int j = 0; j < Height(); j++)
            {
               //Find the angle based on position.
               angle = rotamount*exp(-((i-Width()/2.)*(i-Width()/2)
                  + (j-Height()/2)*(j-Height()/2))/(rotsize*rotsize));
            //Find resampling position.
               ixpos = cos(angle)*(i-Width()/2) + sin(angle)*(j-Height()/2);
               iypos = -sin(angle)*(i-Width()/2) + cos(angle)*(j-Height()/2);
               ixpos += Width()/2;
               iypos += Height()/2;
            
            //Check bounds-- if out of bounds, set to black.
               if ((ixpos >= Width() - 1)
               || (iypos >= Height() - 1)
               || (ixpos < 0)
               || (iypos < 0))
               {
                  tempimage->Pixels(i)[j].SetRed(0);
                  tempimage->Pixels(i)[j].SetRed(0);
                  tempimage->Pixels(i)[j].SetRed(0);
                  continue;
               }  
               
               //Set the values in the temporary image.	   
               tempimage->Pixels(i)[j].SetRed(Pixels(ixpos)[iypos].Red());
               tempimage->Pixels(i)[j].SetBlue(Pixels(ixpos)[iypos].Blue());
               tempimage->Pixels(i)[j].SetGreen(Pixels(ixpos)[iypos].Green());
            }
         
         //Make this image point to the temporary image.
         delete this->pixels;	
         this -> pixels = tempimage->pixels;
         this->npixels = tempimage->npixels;
         this->height = tempimage->height;
         this->width = tempimage->width;
         ApplyGamma(1/2.2);
   	}
   	
      else if (sampling_method == 2)
      {
      //Size and radius of Gaussian kernel (for sampling).
         int size = 5;	
         int radius = (size - 1)/2;
      
      //Temporary image.	
         R2Image* tempimage = new R2Image(Width(), Height());
      
      //Declarations.
         double redsum, bluesum, greensum, dist, 
            ksum, rotamount, rotsize, angle;
         int ixpos, iypos, xlo, xhi, ylo, yhi, filterposx, filterposy;
         double kernel[size][size];
      
      //Parameters.
         rotamount = 10;
         rotsize = 200;
      
      //Apply Gamma correction.
         ApplyGamma(2.2);
      
      //Initialize gaussian kernel w/ size
         for (int i = 0; i < size; i++)
            for (int j = 0; j < size; j++)
               kernel[i][j] = 0;
      
      
         ksum = 0;
         for (int i = 0; i < size; i++)
         {
            for (int j = 0; j < size; j++)
            {
               dist = abs(radius-i);
               kernel[i][j] += exp((-.5)*pow(dist,2.0)/(2*pow((radius/3.),2)));
               dist = abs(radius-j);
               kernel[i][j] += exp((-.5)*pow(dist,2.0)/(2*pow((radius/3.),2)));
               ksum += kernel[i][j];
            }
         }
      
      //Normalize filter kernel.
         for (int i = 0; i < size; i++)
         {
            for (int j = 0; j < size; j++)
            {
               kernel[i][j] = kernel[i][j] / ksum;
            }
         }
      
      //Resample every pixel.
         for (int i = 0; i < Width(); i++)
            for (int j = 0; j < Height(); j++)
            {
            //Find the angle based on position.
               angle = rotamount*exp(-((i-Width()/2.)*(i-Width()/2)
                  + (j-Height()/2)*(j-Height()/2))/(rotsize*rotsize));
            //Find resampling position.
               ixpos = cos(angle)*(i-Width()/2) + sin(angle)*(j-Height()/2);
               iypos = -sin(angle)*(i-Width()/2) + cos(angle)*(j-Height()/2);
               ixpos += Width()/2;
               iypos += Height()/2;
            
            //Check bounds-- if out of bounds, set to black.
               if ((ixpos >= Width() - 2*radius)
               || (iypos >= Height() - 2*radius)
               || (ixpos <= 2*radius)
               || (iypos <= 2*radius))
               {
                  tempimage->Pixels(i)[j].SetRed(0);
                  tempimage->Pixels(i)[j].SetRed(0);
                  tempimage->Pixels(i)[j].SetRed(0);
                  continue;
               }     
            
            //Apply Gaussian filter.
               redsum = bluesum = greensum = ksum = 0;
               xlo = ixpos - radius;
               xhi = ixpos + radius;
               ylo = iypos - radius;
               yhi = iypos + radius;
               for (int ix = xlo; ix <= xhi; ix++) 
               {
                  for (int iy = ylo; iy <= yhi; iy++) 
                  {
                     filterposx = ix - xlo;
                     filterposy = iy - ylo;
                  
                     redsum += kernel[filterposx][filterposy]*Pixels(ix)[iy].Red();
                     bluesum += kernel[filterposx][filterposy]*Pixels(ix)[iy].Blue();
                     greensum += kernel[filterposx][filterposy]*Pixels(ix)[iy].Green();
                     ksum += kernel[filterposx][filterposy];
                  }
               }
            //Normalize color values.
               redsum /= ksum;
               bluesum /= ksum;
               greensum /= ksum;           
                  
            //Set the values in the temporary image.
               tempimage->Pixels(i)[j].SetRed(redsum);
               tempimage->Pixels(i)[j].SetBlue(bluesum);
               tempimage->Pixels(i)[j].SetGreen(greensum);
            }
         
      //Make this image point to the temporary image.
         delete this->pixels;	
         this -> pixels = tempimage->pixels;
         this->npixels = tempimage->npixels;
         this->height = tempimage->height;
         this->width = tempimage->width;
      
      //Apply Gamma correction.
         ApplyGamma(1/2.2);
      }
      
   	
   	else if (sampling_method == 3)
   	{		   
      	//Declarations.
         double rotamount, rotsize, angle;
         double reda, bluea, greena, redb, blueb, greenb;
         double factorx, factory;
         double red, green, blue;
         double dxpos, dypos;
         int ixpos, iypos;
      	
      	//Parameters.
         rotamount = 10;
         rotsize = 200;
      	
           //Temporary image.
         R2Image* tempimage = new R2Image(Width(), Height());
       	//Apply Gamma correction.
         ApplyGamma(2.2);  
      	
         //Sample every pixel in the new image.
         for (int i = 0; i < Width(); i++)
            for (int j = 0; j < Height(); j++)
            {
               //Find the angle based on position.
               angle = rotamount*exp(-((i-Width()/2.)*(i-Width()/2)
                  + (j-Height()/2)*(j-Height()/2))/(rotsize*rotsize));
            //Find resampling position.
               dxpos = cos(angle)*(i-Width()/2) + sin(angle)*(j-Height()/2);
               dypos = -sin(angle)*(i-Width()/2) + cos(angle)*(j-Height()/2);
               dxpos += Width()/2;
               dypos += Height()/2;
            
            	ixpos = (int)dxpos;
            	iypos = (int)dypos;
            
            //Check bounds-- if out of bounds, set to black.
               if ((ixpos >= Width() - 1)
               || (iypos >= Height() - 1)
               || (ixpos < 0)
               || (iypos < 0))
               {
                  tempimage->Pixels(i)[j].SetRed(0);
                  tempimage->Pixels(i)[j].SetRed(0);
                  tempimage->Pixels(i)[j].SetRed(0);
                  continue;
               }  
               
               //Bilerp each color.
               reda = factorx*Pixels(ixpos)[iypos].Red() 
                    + (1 - factorx)*Pixels(ixpos+1)[iypos].Red();
               redb = factorx*Pixels(ixpos)[iypos+1].Red() 
                    + (1 - factorx)*Pixels(ixpos+1)[iypos+1].Red();
               red = factory*reda + (1 - factory)*redb;
            	
            	
               bluea = factorx*Pixels(ixpos)[iypos].Blue() 
                     + (1 - factorx)*Pixels(ixpos+1)[iypos].Blue();
               blueb = factorx*Pixels(ixpos)[iypos+1].Blue() 
                     + (1 - factorx)*Pixels(ixpos+1)[iypos+1].Blue();
               blue = factory*bluea + (1 - factory)*blueb;
            	
            	
               greena = factorx*Pixels(ixpos)[iypos].Green() 
                      + (1 - factorx)*Pixels(ixpos+1)[iypos].Green();
               greenb = factorx*Pixels(ixpos)[iypos+1].Green() 
                  	 + (1 - factorx)*Pixels(ixpos+1)[iypos+1].Green();
               green = factory*greena + (1 - factory)*greenb;
            	
            	//Set the values in the temporary image.    
               tempimage->Pixels(i)[j].SetRed(red);
               tempimage->Pixels(i)[j].SetBlue(blue);
               tempimage->Pixels(i)[j].SetGreen(green);
            }
         
         //Make this image point to the temporary image.
         delete this->pixels;	
         this -> pixels = tempimage->pixels;
         this->npixels = tempimage->npixels;
         this->height = tempimage->height;
         this->width = tempimage->width;
         ApplyGamma(1/2.2);
   	}
   	
   }

// Dither operations ////////////////////////////////////////////////

   void R2Image::
   Quantize (int nbits)
   {
   // Quantizes an image with "nbits" bits per channel.
   
   	//Declarations.
      double floor;
      double mid;
      double next;
   	
   	//Snap all pixels to range.
      for (int r = 0; r < Width(); r++)
      {
         for (int c = 0; c < Height(); c++)
         {
            for (double i = 0; i < 1; i += 1/(pow(2., (double)nbits)-1))	
            {
            	//Compute midpoint, floor, and ceiling(next) for each pixel.
               floor = i;
               mid = i + 1/((pow(2., (double)nbits)-1)*2);
               next = i + 1/(pow(2., (double)nbits)-1);
            	
            	//Snap pixels to the floor if < mid and ceiling if > mid
               if ((Pixels(r)[c].Red() >= floor) && (Pixels(r)[c].Red() < mid))
                  Pixels(r)[c].SetRed(floor);
               else if ((Pixels(r)[c].Red() >= mid) && (Pixels(r)[c].Red() < next))
                  Pixels(r)[c].SetRed(next);
            		
               if ((Pixels(r)[c].Green() >= floor) && (Pixels(r)[c].Green() < mid))
                  Pixels(r)[c].SetGreen(floor);
               else if ((Pixels(r)[c].Green() >= mid) && (Pixels(r)[c].Green() < next))
                  Pixels(r)[c].SetGreen(next);
            		
               if ((Pixels(r)[c].Blue() >= floor) && (Pixels(r)[c].Blue() < mid))
                  Pixels(r)[c].SetBlue(floor);
               else if ((Pixels(r)[c].Blue() >= mid) && (Pixels(r)[c].Blue() < next))
                  Pixels(r)[c].SetBlue(next);
            }
         }
      }
   }



   void R2Image::
   RandomDither(int nbits)
   {
   // Converts and image to nbits per channel using random dither.
   
      //Random number declaration
      double random; 
   	
   	//Add random noise to every pixel.
      for (int r = 0; r < Width(); r++)
      {
         for (int c = 0; c < Height(); c++)
         {
         	//Generate a random number (noise).
            random = (((double)rand()/(RAND_MAX))/(2*nbits) - 1/(2.*2*nbits));
            
         	//Adjust pixel values by adding noise.
            Pixels(r)[c].SetRed(Pixels(r)[c].Red() + random);
            Pixels(r)[c].SetBlue(Pixels(r)[c].Blue() + random);
            Pixels(r)[c].SetGreen(Pixels(r)[c].Green() + random);
            Pixels(r)[c].Clamp();
         }
      }
      
      //Quantize noisy pixels.
      Quantize(nbits);
   }



   void R2Image::
   OrderedDither(int nbits)
   {
   // Converts an image to nbits per channel using ordered dither, 
   // with a 4x4 Bayer's pattern matrix.
   
   	//Positions in Bayer matrix.
      int i,j;
      //Initialize Bayer matrix.
      double Bayer4[4][4] = {{15,7,13,5},{3,11,1,9},{12,4,14,6},{0,8,2,10}};
   
   	//Apply filter to all pixels.
      for (int r = 0; r < Width(); r++)
      {
         for (int c = 0; c < Height(); c++)
         {
            i = r%4;
            j = c%4;
            Pixels(r)[c].SetRed(Pixels(r)[c].Red() + Bayer4[i][j]/110.);
            Pixels(r)[c].SetBlue(Pixels(r)[c].Blue() + Bayer4[i][j]/110.);
            Pixels(r)[c].SetGreen(Pixels(r)[c].Green() + Bayer4[i][j]/110.);
            Pixels(r)[c].Clamp();
         }
      }
      
      //Quantize after filtering.
      Quantize(nbits);
   }



   void R2Image::
   FloydSteinbergDither(int nbits)
   {
   // Converts an image to nbits per channel using Floyd-Steinberg dither.
   // with error diffusion.
   
   	//Declarations.
      double floor, mid, next;
      double oldred, oldblue, oldgreen;
      double quanterror = 0;
   	
   	//Apply error diffusion to each pixel.
      for (int r = 0; r < Width()-1; r++)
      {
         for (int c = 0; c < Height()-1; c++)
         {
            oldred = Pixels(r)[c].Red();
            oldblue = Pixels(r)[c].Blue();
            oldgreen = Pixels(r)[c].Green();
         	
         	//Quantize each pixel and compute quantization error.
            for (double i = 0; i < 1; i += 1/(pow(2., (double)nbits)-1))	
            {
               floor = i;
               mid = i + 1/((pow(2., (double)nbits)-1)*2);
               next = i + 1/(pow(2., (double)nbits)-1);
            	
               if ((Pixels(r)[c].Red() >= floor) && (Pixels(r)[c].Red() < mid))
                  Pixels(r)[c].SetRed(floor);
               else if ((Pixels(r)[c].Red() >= mid) && (Pixels(r)[c].Red() < next))
                  Pixels(r)[c].SetRed(next);
            		
               if ((Pixels(r)[c].Green() >= floor) && (Pixels(r)[c].Green() < mid))
                  Pixels(r)[c].SetGreen(floor);
               else if ((Pixels(r)[c].Green() >= mid) && (Pixels(r)[c].Green() < next))
                  Pixels(r)[c].SetGreen(next);
            		
               if ((Pixels(r)[c].Blue() >= floor) && (Pixels(r)[c].Blue() < mid))
                  Pixels(r)[c].SetBlue(floor);
               else if ((Pixels(r)[c].Blue() >= mid) && (Pixels(r)[c].Blue() < next))
                  Pixels(r)[c].SetBlue(next);
            }
            Pixels(r)[c].Clamp();
         	
         	//Interpolate between quantization error and actual pixel value
         	//on neighboring pixels.
            quanterror = oldred - Pixels(r)[c].Red();
            Pixels(r+1)[c].SetRed(Pixels(r+1)[c].Red() + 7./16 * quanterror);
            Pixels(r-1)[c+1].SetRed(Pixels(r-1)[c+1].Red() + 3./16 * quanterror);
            Pixels(r)[c+1].SetRed(Pixels(r)[c+1].Red() + 5./16 * quanterror);
            Pixels(r+1)[c+1].SetRed(Pixels(r+1)[c+1].Red() + 1./16 * quanterror);
         	
            quanterror = oldblue - Pixels(r)[c].Blue();
            Pixels(r+1)[c].SetBlue(Pixels(r+1)[c].Blue() + 7./16 * quanterror);
            Pixels(r-1)[c+1].SetBlue(Pixels(r-1)[c+1].Blue() + 3./16 * quanterror);
            Pixels(r)[c+1].SetBlue(Pixels(r)[c+1].Blue() + 5./16 * quanterror);
            Pixels(r+1)[c+1].SetBlue(Pixels(r+1)[c+1].Blue() + 1./16 * quanterror);
         	
            quanterror = oldgreen - Pixels(r)[c].Green();
            Pixels(r+1)[c].SetGreen(Pixels(r+1)[c].Green() + 7./16 * quanterror);
            Pixels(r-1)[c+1].SetGreen(Pixels(r-1)[c+1].Green() + 3./16 * quanterror);
            Pixels(r)[c+1].SetGreen(Pixels(r)[c+1].Green() + 5./16 * quanterror);
            Pixels(r+1)[c+1].SetGreen(Pixels(r+1)[c+1].Green() + 1./16 * quanterror);
         }
      }
   /*********************************************************************/
   	//Fix the last row of pixels.
      int r = Width() - 1;
      for (int c = 0; c < Height(); c++)
      {
         oldred = Pixels(r)[c].Red();
         oldblue = Pixels(r)[c].Blue();
         oldgreen = Pixels(r)[c].Green();
         	
         	//Quantize each pixel and compute quantization error.
         for (double i = 0; i < 1; i += 1/(pow(2., (double)nbits)-1))	
         {
            floor = i;
            mid = i + 1/((pow(2., (double)nbits)-1)*2);
            next = i + 1/(pow(2., (double)nbits)-1);
            	
            if ((Pixels(r)[c].Red() >= floor) && (Pixels(r)[c].Red() < mid))
               Pixels(r)[c].SetRed(floor);
            else if ((Pixels(r)[c].Red() >= mid) && (Pixels(r)[c].Red() < next))
               Pixels(r)[c].SetRed(next);
            		
            if ((Pixels(r)[c].Green() >= floor) && (Pixels(r)[c].Green() < mid))
               Pixels(r)[c].SetGreen(floor);
            else if ((Pixels(r)[c].Green() >= mid) && (Pixels(r)[c].Green() < next))
               Pixels(r)[c].SetGreen(next);
            		
            if ((Pixels(r)[c].Blue() >= floor) && (Pixels(r)[c].Blue() < mid))
               Pixels(r)[c].SetBlue(floor);
            else if ((Pixels(r)[c].Blue() >= mid) && (Pixels(r)[c].Blue() < next))
               Pixels(r)[c].SetBlue(next);
         }
         Pixels(r)[c].Clamp();	
      }
      int c = Height() - 1;
      for (int r = 0; r < Width(); r++)
      {
         oldred = Pixels(r)[c].Red();
         oldblue = Pixels(r)[c].Blue();
         oldgreen = Pixels(r)[c].Green();
         	
         	//Quantize each pixel and compute quantization error.
         for (double i = 0; i < 1; i += 1/(pow(2., (double)nbits)-1))	
         {
            floor = i;
            mid = i + 1/((pow(2., (double)nbits)-1)*2);
            next = i + 1/(pow(2., (double)nbits)-1);
            	
            if ((Pixels(r)[c].Red() >= floor) && (Pixels(r)[c].Red() < mid))
               Pixels(r)[c].SetRed(floor);
            else if ((Pixels(r)[c].Red() >= mid) && (Pixels(r)[c].Red() < next))
               Pixels(r)[c].SetRed(next);
            		
            if ((Pixels(r)[c].Green() >= floor) && (Pixels(r)[c].Green() < mid))
               Pixels(r)[c].SetGreen(floor);
            else if ((Pixels(r)[c].Green() >= mid) && (Pixels(r)[c].Green() < next))
               Pixels(r)[c].SetGreen(next);
            		
            if ((Pixels(r)[c].Blue() >= floor) && (Pixels(r)[c].Blue() < mid))
               Pixels(r)[c].SetBlue(floor);
            else if ((Pixels(r)[c].Blue() >= mid) && (Pixels(r)[c].Blue() < next))
               Pixels(r)[c].SetBlue(next);
         }
         Pixels(r)[c].Clamp();	
      }
      /***********************************************************************/
   }



// Miscellaneous operations ////////////////////////////////////////////////

   void R2Image::
   CopyChannel (const R2Image& from_image, int from_channel, int to_channel)
   {
   // Copies one channel of an image (e.g., R2_IMAGE_RED_CHANNEL).  
   // to another channel
   
   	//Cast away constness of parameter (avoid compiler warning).
      R2Image nonconst = (R2Image)from_image;
      
   	//Copy channel for each pixel.
      for (int i = 0; i < Width(); i++)
      {
         for (int j = 0; j < Height(); j++)
         {
            if ((i < from_image.Width()) && (j < from_image.Height()))
               Pixels(i)[j][to_channel] = nonconst.Pixels(i)[j][from_channel];	
         }
      }
   }

   void R2Image::
   Composite(const R2Image& top, int operation)
   {
   // Composite passed image on top of this one using operation OVER
   
   	//Declarations.
      double topalpha, topred, topblue, topgreen;    
      R2Image nonconst = (R2Image)top;
   	
   	//Apply the over operation on each pixel.
      for (int i = 0; i < Width(); i++)
      {
         for (int j = 0; j < Height(); j++)
         {
            if ((i < top.Width()) && (j < top.Height()))
            {
            	//Store alpha value of over image in variables (easier).
               topalpha = nonconst.Pixels(i)[j].Alpha();
               topred = nonconst.Pixels(i)[j].Red();
               topblue = nonconst.Pixels(i)[j].Blue();
               topgreen = nonconst.Pixels(i)[j].Green();
            
            	//Apply the over operator.
               Pixels(i)[j].SetRed((1-topalpha)*topred 
                  + (topalpha)*Pixels(i)[j].Red());
               Pixels(i)[j].SetBlue((1-topalpha)*topblue 
                  + (topalpha)*Pixels(i)[j].Blue());
               Pixels(i)[j].SetGreen((1-topalpha)*topgreen 
                  + (topalpha)*Pixels(i)[j].Green());
               Pixels(i)[j].SetAlpha((1-topalpha)*topalpha 
                  + (topalpha)*Pixels(i)[j].Alpha());
            }
         }
      }
   }

   void R2Image::
   Morph(const R2Image& target, 
   R2Segment *source_segments, R2Segment *target_segments, int nsegments, 
   double t, int sampling_method)
   {
   // Morph this source image towards a passed target image by t using pairwise line segment correspondences
   
   	//Declarations.
      R2Point X, Xp, DSum;
      R2Segment PQ, PpQp, PX;
      R2Vector Disp;
      double u, v, a, b, p, dist, weight, weightsum;
      int xpos, ypos;
      //Cast away constness of parameter (avoid compiler warnings.)
      R2Image targ = (R2Image)target;
      //New arrays for interpolated line segments
      R2Segment interpsource[nsegments];
      R2Segment interptarg[nsegments];
   	
   	//Temporary images.
      R2Image* tempimage = new R2Image(Width(), Height());
      R2Image* tempimage2 = new R2Image(targ.Width(), targ.Height());
      
   	//Tweakable parameters.
      a = .8;
      b = .9;
      p = .1;
      
   	//Apply Gamma correction.
      ApplyGamma(2.2);
   	
   	//Interpolate source lines and target lines based on time step.
   	//Store lines in a new array.
      for (int k = 0; k < nsegments; k++)
      {
         interpsource[k].SetStart(t*source_segments[k].Start() 
            + (1-t)*target_segments[k].Start());
         interpsource[k].SetEnd(t*source_segments[k].End() 
            + (1-t)*target_segments[k].End());
      }
   	
   	//Apply Beier algoritm to every pixel in destination.
      for (int r = 0; r < targ.Width(); r++)
      {
         for (int c = 0; c < targ.Height(); c++)
         {
         	//Reset variables
            X.Reset(r,c);
            DSum.Reset(0,0);
            weightsum = 0;
            
            //Iterate over every line.
            for (int i = 0; i < nsegments; i++)
            {
            	//Initialize segments.
               PQ = target_segments[i];
               PpQp = interpsource[i];
               PX.SetStart(PQ.Start());
               PX.SetEnd(X);
            
            	//Find u and v values
               u = (PX.Vector()*PX.Length()).Dot(PQ.Vector()
                  * PQ.Length())/pow(PQ.Length(),2);
               v = (PX.Vector()*PX.Length()).Dot(PQ.Normal()
                  * PQ.Length())/PQ.Length();
            
               Xp = PpQp.Start() + u*(PpQp.Vector()*PpQp.Length())
                  + v*(PpQp.Normal()*PpQp.Length())/PpQp.Length();
                  
            	//Compute displacement vector and distance to line.
               Disp = Xp - X;
               dist = abs(PQ.T(X));	
            	
            	//Compute weight of line.
               weight = (pow(PQ.Length(),p)/pow(a + dist, b));
            	
            	//Add displacement vector to sum.
               DSum += Disp*weight;
             
             	//Add line weight to sum.  
               weightsum += weight;
            }
            //Adjust Xp.
            Xp = X + DSum/weightsum;
            xpos = Xp.X();
            ypos = Xp.Y();
            
         	//Make sure still in bounds-- black if not.
            if (xpos < 0 || ypos < 0 || xpos >= (Width()-1) || ypos >= (Height()-1))
            {
               tempimage->Pixels(r)[c].SetRed(0);
               tempimage->Pixels(r)[c].SetBlue(0);
               tempimage->Pixels(r)[c].SetGreen(0);
               continue;
            } 
         	
         	//Set value in temporary image.
            tempimage->Pixels(r)[c].SetRed(Pixels(xpos)[ypos].Red());
            tempimage->Pixels(r)[c].SetBlue(Pixels(xpos)[ypos].Blue());
            tempimage->Pixels(r)[c].SetGreen(Pixels(xpos)[ypos].Green());
         }
      }
      
   	//Interpolate lines for the destination image.
      for (int k = 0; k < nsegments; k++)
      {
         interptarg[k].SetStart(t*source_segments[k].Start() + (1-t)*target_segments[k].Start());
         interptarg[k].SetEnd(t*source_segments[k].End() + (1-t)*target_segments[k].End());
      }
   	
   	//Apply Beier algorithm to every pixel.
      for (int r = 0; r < Width(); r++)
      {
         for (int c = 0; c < Height(); c++)
         {
         	//Same process.
            X.Reset(r,c);
            DSum.Reset(0,0);
            weightsum = 0;
            for (int i = 0; i < nsegments; i++)
            {
               PQ = source_segments[i];
               PpQp = interptarg[i];
               PX.SetStart(PQ.Start());
               PX.SetEnd(X);
            
               u = (PX.Vector()*PX.Length()).Dot(PQ.Vector()*PQ.Length())/pow(PQ.Length(),2);
               v = (PX.Vector()*PX.Length()).Dot(PQ.Normal()*PQ.Length())/PQ.Length();
            
               Xp = PpQp.Start() + u*(PpQp.Vector()*PpQp.Length())
                  + v*(PpQp.Normal()*PpQp.Length())/PpQp.Length();
                  
               Disp = Xp - X;
               dist = abs(PQ.T(X));	
            	
               weight = (pow(PQ.Length(),p)/pow(a + dist, b));
            	
               DSum += Disp*weight;
            	
               weightsum += weight;
            }
            Xp = X + DSum/weightsum;
            xpos = Xp.X();
            ypos = Xp.Y();
            
         	//Check bounts again-- black if not inside.
            if (xpos < 0 || ypos < 0 || xpos >= 
            	(targ.Width()-1) || ypos >= (targ.Height()-1))
            {
               tempimage->Pixels(r)[c].SetRed(0);
               tempimage->Pixels(r)[c].SetBlue(0);
               tempimage->Pixels(r)[c].SetGreen(0);
               continue;
            }
            
         	//Set values in temporary image.
            tempimage2->Pixels(r)[c].SetRed(targ.Pixels(xpos)[ypos].Red());
            tempimage2->Pixels(r)[c].SetBlue(targ.Pixels(xpos)[ypos].Blue());
            tempimage2->Pixels(r)[c].SetGreen(targ.Pixels(xpos)[ypos].Green());
         }
      }
   	
   	//Interpolate pixels between source morph and target morph.
      for (int r = 0; r < Width(); r++)
      {
         for (int c = 0; c < Height(); c++)
         {
            tempimage->Pixels(r)[c].SetRed((1-t)*tempimage->Pixels(r)[c].Red()
               + t*tempimage2->Pixels(r)[c].Red());
            tempimage->Pixels(r)[c].SetBlue((1-t)*tempimage->Pixels(r)[c].Blue()
               + t*tempimage2->Pixels(r)[c].Blue());
            tempimage->Pixels(r)[c].SetGreen((1-t)*tempimage->Pixels(r)[c].Green()
               + t*tempimage2->Pixels(r)[c].Green());
         }
      }
      
   	//Use bilateral filter to smooth out artifacts.
      tempimage->BilateralFilter(.5,.1); 
   	
   	//Set this image to point to the temporary image.
      delete this->pixels;	
      this -> pixels = tempimage->pixels;
      this->npixels = tempimage->npixels;
      this->height = tempimage->height;
      this->width = tempimage->width;
      
   	//Apply Gamma correction.
      ApplyGamma(1/2.2);
   }

   void R2Image::
   Crop(int x, int y, int w, int h)
   {
   // Extracts a sub image from the image, 
   // at position (x, y), width w, and height h.
   
   // MAY FILL IN IMPLEMENTATION HERE (REMOVE PRINT STATEMENT WHEN DONE) (NO CREDIT FOR ASSIGNMENT)
      fprintf(stderr, "Crop(%d %d  %d %d) not implemented\n", x, y, w, h);
   }

   void R2Image::
   Add(const R2Image& image)
   {
   // Add passed image pixel-by-pixel.
   
   // MAY FILL IN IMPLEMENTATION HERE (REMOVE PRINT STATEMENT WHEN DONE) (NO CREDIT FOR ASSIGNMENT)
      fprintf(stderr, "Add not implemented\n");
   }

   void R2Image::
   Subtract(const R2Image& image)
   {
   // Subtract passed image from this image.
   
   // MAY FILL IN IMPLEMENTATION HERE (REMOVE PRINT STATEMENT WHEN DONE) (NO CREDIT FOR ASSIGNMENT)
      fprintf(stderr, "Subtract not implemented\n");
   }


////////////////////////////////////////////////////////////////////////
// I/O Functions
////////////////////////////////////////////////////////////////////////

   int R2Image::
   Read(const char *filename)
   {
   // Initialize everything
      if (pixels) { delete [] pixels; pixels = NULL; }
      npixels = width = height = 0;
   
   // Parse input filename extension
      char *input_extension;
      if (!(input_extension = (char*)strrchr(filename, '.'))) {
         fprintf(stderr, "Input file has no extension (e.g., .jpg).\n");
         return 0;
      }
   
   // Read file of appropriate type
      if (!strncmp(input_extension, ".bmp", 4)) 
         return ReadBMP(filename);
      else if (!strncmp(input_extension, ".ppm", 4)) 
         return ReadPPM(filename);
      else if (!strncmp(input_extension, ".jpg", 4)) 
         return ReadJPEG(filename);
      else if (!strncmp(input_extension, ".jpeg", 5)) 
         return ReadJPEG(filename);
   
   // Should never get here
      fprintf(stderr, "Unrecognized image file extension");
      return 0;
   }



   int R2Image::
   Write(const char *filename) const
   {
   // Parse input filename extension
      char *input_extension;
      if (!(input_extension = (char*)strrchr(filename, '.'))) {
         fprintf(stderr, "Input file has no extension (e.g., .jpg).\n");
         return 0;
      }
   
   // Write file of appropriate type
      if (!strncmp(input_extension, ".bmp", 4)) 
         return WriteBMP(filename);
      else if (!strncmp(input_extension, ".ppm", 4)) 
         return WritePPM(filename, 1);
      else if (!strncmp(input_extension, ".jpg", 5)) 
         return WriteJPEG(filename);
      else if (!strncmp(input_extension, ".jpeg", 5)) 
         return WriteJPEG(filename);
   
   // Should never get here
      fprintf(stderr, "Unrecognized image file extension");
      return 0;
   }



////////////////////////////////////////////////////////////////////////
// BMP I/O
////////////////////////////////////////////////////////////////////////

#if !defined(_WIN32)

   typedef 
      struct tagBITMAPFILEHEADER {
         unsigned short int bfType;
         unsigned int bfSize;
         unsigned short int bfReserved1;
         unsigned short int bfReserved2;
         unsigned int bfOffBits;
      } BITMAPFILEHEADER;

   typedef 
      struct tagBITMAPINFOHEADER {
         unsigned int biSize;
         int biWidth;
         int biHeight;
         unsigned short int biPlanes;
         unsigned short int biBitCount;
         unsigned int biCompression;
         unsigned int biSizeImage;
         int biXPelsPerMeter;
         int biYPelsPerMeter;
         unsigned int biClrUsed;
         unsigned int biClrImportant;
      } BITMAPINFOHEADER;

   typedef 
      struct tagRGBTRIPLE {
         unsigned char rgbtBlue;
         unsigned char rgbtGreen;
         unsigned char rgbtRed;
      } RGBTRIPLE;

   typedef 
      struct tagRGBQUAD {
         unsigned char rgbBlue;
         unsigned char rgbGreen;
         unsigned char rgbRed;
         unsigned char rgbReserved;
      } RGBQUAD;

#endif

#define BI_RGB        0L
#define BI_RLE8       1L
#define BI_RLE4       2L
#define BI_BITFIELDS  3L

#define BMP_BF_TYPE 0x4D42 /* word BM */
#define BMP_BF_OFF_BITS 54 /* 14 for file header + 40 for info header (not sizeof(), but packed size) */
#define BMP_BI_SIZE 40 /* packed size of info header */


   static unsigned short int WordReadLE(FILE *fp)
   {
   // Read a unsigned short int from a file in little endian format 
      unsigned short int lsb, msb;
      lsb = getc(fp);
      msb = getc(fp);
      return (msb << 8) | lsb;
   }



   static void WordWriteLE(unsigned short int x, FILE *fp)
   {
   // Write a unsigned short int to a file in little endian format
      unsigned char lsb = (unsigned char) (x & 0x00FF); putc(lsb, fp); 
      unsigned char msb = (unsigned char) (x >> 8); putc(msb, fp);
   }



   static unsigned int DWordReadLE(FILE *fp)
   {
   // Read a unsigned int word from a file in little endian format 
      unsigned int b1 = getc(fp);
      unsigned int b2 = getc(fp);
      unsigned int b3 = getc(fp);
      unsigned int b4 = getc(fp);
      return (b4 << 24) | (b3 << 16) | (b2 << 8) | b1;
   }



   static void DWordWriteLE(unsigned int x, FILE *fp)
   {
   // Write a unsigned int to a file in little endian format 
      unsigned char b1 = (x & 0x000000FF); putc(b1, fp);
      unsigned char b2 = ((x >> 8) & 0x000000FF); putc(b2, fp);
      unsigned char b3 = ((x >> 16) & 0x000000FF); putc(b3, fp);
      unsigned char b4 = ((x >> 24) & 0x000000FF); putc(b4, fp);
   }



   static int LongReadLE(FILE *fp)
   {
   // Read a int word from a file in little endian format 
      int b1 = getc(fp);
      int b2 = getc(fp);
      int b3 = getc(fp);
      int b4 = getc(fp);
      return (b4 << 24) | (b3 << 16) | (b2 << 8) | b1;
   }



   static void LongWriteLE(int x, FILE *fp)
   {
   // Write a int to a file in little endian format 
      char b1 = (x & 0x000000FF); putc(b1, fp);
      char b2 = ((x >> 8) & 0x000000FF); putc(b2, fp);
      char b3 = ((x >> 16) & 0x000000FF); putc(b3, fp);
      char b4 = ((x >> 24) & 0x000000FF); putc(b4, fp);
   }



   int R2Image::
   ReadBMP(const char *filename)
   {
   // Open file
      FILE *fp = fopen(filename, "rb");
      if (!fp) {
         fprintf(stderr, "Unable to open image file: %s\n", filename);
         return 0;
      }
   
   /* Read file header */
      BITMAPFILEHEADER bmfh;
      bmfh.bfType = WordReadLE(fp);
      bmfh.bfSize = DWordReadLE(fp);
      bmfh.bfReserved1 = WordReadLE(fp);
      bmfh.bfReserved2 = WordReadLE(fp);
      bmfh.bfOffBits = DWordReadLE(fp);
   
   /* Check file header */
      assert(bmfh.bfType == BMP_BF_TYPE);
   /* ignore bmfh.bfSize */
   /* ignore bmfh.bfReserved1 */
   /* ignore bmfh.bfReserved2 */
      assert(bmfh.bfOffBits == BMP_BF_OFF_BITS);
   
   /* Read info header */
      BITMAPINFOHEADER bmih;
      bmih.biSize = DWordReadLE(fp);
      bmih.biWidth = LongReadLE(fp);
      bmih.biHeight = LongReadLE(fp);
      bmih.biPlanes = WordReadLE(fp);
      bmih.biBitCount = WordReadLE(fp);
      bmih.biCompression = DWordReadLE(fp);
      bmih.biSizeImage = DWordReadLE(fp);
      bmih.biXPelsPerMeter = LongReadLE(fp);
      bmih.biYPelsPerMeter = LongReadLE(fp);
      bmih.biClrUsed = DWordReadLE(fp);
      bmih.biClrImportant = DWordReadLE(fp);
   
   // Check info header 
      assert(bmih.biSize == BMP_BI_SIZE);
      assert(bmih.biWidth > 0);
      assert(bmih.biHeight > 0);
      assert(bmih.biPlanes == 1);
      assert(bmih.biBitCount == 24);  /* RGB */
      assert(bmih.biCompression == BI_RGB);   /* RGB */
      int lineLength = bmih.biWidth * 3;  /* RGB */
      if ((lineLength % 4) != 0) lineLength = (lineLength / 4 + 1) * 4;
      assert(bmih.biSizeImage == (unsigned int) lineLength * (unsigned int) bmih.biHeight);
   
   // Assign width, height, and number of pixels
      width = bmih.biWidth;
      height = bmih.biHeight;
      npixels = width * height;
   
   // Allocate unsigned char buffer for reading pixels
      int rowsize = 3 * width;
      if ((rowsize % 4) != 0) rowsize = (rowsize / 4 + 1) * 4;
      int nbytes = bmih.biSizeImage;
      unsigned char *buffer = new unsigned char [nbytes];
      if (!buffer) {
         fprintf(stderr, "Unable to allocate temporary memory for BMP file");
         fclose(fp);
         return 0;
      }
   
   // Read buffer 
      fseek(fp, (long) bmfh.bfOffBits, SEEK_SET);
      if (fread(buffer, 1, bmih.biSizeImage, fp) != bmih.biSizeImage) {
         fprintf(stderr, "Error while reading BMP file %s", filename);
         return 0;
      }
   
   // Close file
      fclose(fp);
   
   // Allocate pixels for image
      pixels = new R2Pixel [ width * height ];
      if (!pixels) {
         fprintf(stderr, "Unable to allocate memory for BMP file");
         fclose(fp);
         return 0;
      }
   
   // Assign pixels
      for (int j = 0; j < height; j++) {
         unsigned char *p = &buffer[j * rowsize];
         for (int i = 0; i < width; i++) {
            double b = (double) *(p++) / 255;
            double g = (double) *(p++) / 255;
            double r = (double) *(p++) / 255;
            R2Pixel pixel(r, g, b, 1);
            SetPixel(i, j, pixel);
         }
      }
   
   // Free unsigned char buffer for reading pixels
      delete [] buffer;
   
   // Return success
      return 1;
   }



   int R2Image::
   WriteBMP(const char *filename) const
   {
   // Open file
      FILE *fp = fopen(filename, "wb");
      if (!fp) {
         fprintf(stderr, "Unable to open image file: %s\n", filename);
         return 0;
      }
   
   // Compute number of bytes in row
      int rowsize = 3 * width;
      if ((rowsize % 4) != 0) rowsize = (rowsize / 4 + 1) * 4;
   
   // Write file header 
      BITMAPFILEHEADER bmfh;
      bmfh.bfType = BMP_BF_TYPE;
      bmfh.bfSize = BMP_BF_OFF_BITS + rowsize * height;
      bmfh.bfReserved1 = 0;
      bmfh.bfReserved2 = 0;
      bmfh.bfOffBits = BMP_BF_OFF_BITS;
      WordWriteLE(bmfh.bfType, fp);
      DWordWriteLE(bmfh.bfSize, fp);
      WordWriteLE(bmfh.bfReserved1, fp);
      WordWriteLE(bmfh.bfReserved2, fp);
      DWordWriteLE(bmfh.bfOffBits, fp);
   
   // Write info header 
      BITMAPINFOHEADER bmih;
      bmih.biSize = BMP_BI_SIZE;
      bmih.biWidth = width;
      bmih.biHeight = height;
      bmih.biPlanes = 1;
      bmih.biBitCount = 24;       /* RGB */
      bmih.biCompression = BI_RGB;    /* RGB */
      bmih.biSizeImage = rowsize * (unsigned int) bmih.biHeight;  /* RGB */
      bmih.biXPelsPerMeter = 2925;
      bmih.biYPelsPerMeter = 2925;
      bmih.biClrUsed = 0;
      bmih.biClrImportant = 0;
      DWordWriteLE(bmih.biSize, fp);
      LongWriteLE(bmih.biWidth, fp);
      LongWriteLE(bmih.biHeight, fp);
      WordWriteLE(bmih.biPlanes, fp);
      WordWriteLE(bmih.biBitCount, fp);
      DWordWriteLE(bmih.biCompression, fp);
      DWordWriteLE(bmih.biSizeImage, fp);
      LongWriteLE(bmih.biXPelsPerMeter, fp);
      LongWriteLE(bmih.biYPelsPerMeter, fp);
      DWordWriteLE(bmih.biClrUsed, fp);
      DWordWriteLE(bmih.biClrImportant, fp);
   
   // Write image, swapping blue and red in each pixel
      int pad = rowsize - width * 3;
      for (int j = 0; j < height; j++) {
         for (int i = 0; i < width; i++) {
            const R2Pixel& pixel = (*this)[i][j];
            double r = 255.0 * pixel.Red();
            double g = 255.0 * pixel.Green();
            double b = 255.0 * pixel.Blue();
            if (r >= 255) r = 255;
            if (g >= 255) g = 255;
            if (b >= 255) b = 255;
            fputc((unsigned char) b, fp);
            fputc((unsigned char) g, fp);
            fputc((unsigned char) r, fp);
         }
      
      // Pad row
         for (int i = 0; i < pad; i++) fputc(0, fp);
      }
   
   // Close file
      fclose(fp);
   
   // Return success
      return 1;  
   }



////////////////////////////////////////////////////////////////////////
// PPM I/O
////////////////////////////////////////////////////////////////////////

   int R2Image::
   ReadPPM(const char *filename)
   {
   // Open file
      FILE *fp = fopen(filename, "rb");
      if (!fp) {
         fprintf(stderr, "Unable to open image file: %s\n", filename);
         return 0;
      }
   
   // Read PPM file magic identifier
      char buffer[128];
      if (!fgets(buffer, 128, fp)) {
         fprintf(stderr, "Unable to read magic id in PPM file");
         fclose(fp);
         return 0;
      }
   
   // skip comments
      int c = getc(fp);
      while (c == '#') {
         while (c != '\n') c = getc(fp);
         c = getc(fp);
      }
      ungetc(c, fp);
   
   // Read width and height
      if (fscanf(fp, "%d%d", &width, &height) != 2) {
         fprintf(stderr, "Unable to read width and height in PPM file");
         fclose(fp);
         return 0;
      }
   
      npixels = width * height;
   
   // Read max value
      double max_value;
      if (fscanf(fp, "%lf", &max_value) != 1) {
         fprintf(stderr, "Unable to read max_value in PPM file");
         fclose(fp);
         return 0;
      }
   
   // Allocate image pixels
      pixels = new R2Pixel [ width * height ];
      if (!pixels) {
         fprintf(stderr, "Unable to allocate memory for PPM file");
         fclose(fp);
         return 0;
      }
   
   // Check if raw or ascii file
      if (!strcmp(buffer, "P6\n")) {
      // Read up to one character of whitespace (\n) after max_value
         int c = getc(fp);
         if (!isspace(c)) putc(c, fp);
      
      // Read raw image data 
      // First ppm pixel is top-left, so read in opposite scan-line order
         for (int j = height-1; j >= 0; j--) {
            for (int i = 0; i < width; i++) {
               double r = (double) getc(fp) / max_value;
               double g = (double) getc(fp) / max_value;
               double b = (double) getc(fp) / max_value;
               R2Pixel pixel(r, g, b, 1);
               SetPixel(i, j, pixel);
            }
         }
      }
      else {
      // Read asci image data 
      // First ppm pixel is top-left, so read in opposite scan-line order
         for (int j = height-1; j >= 0; j--) {
            for (int i = 0; i < width; i++) {
            // Read pixel values
               int red, green, blue;
               if (fscanf(fp, "%d%d%d", &red, &green, &blue) != 3) {
                  fprintf(stderr, "Unable to read data at (%d,%d) in PPM file", i, j);
                  fclose(fp);
                  return 0;
               }
            
            // Assign pixel values
               double r = (double) red / max_value;
               double g = (double) green / max_value;
               double b = (double) blue / max_value;
               R2Pixel pixel(r, g, b, 1);
               SetPixel(i, j, pixel);
            }
         }
      }
   
   // Close file
      fclose(fp);
   
   // Return success
      return 1;
   }



   int R2Image::
   WritePPM(const char *filename, int ascii) const
   {
   // Check type
      if (ascii) {
      // Open file
         FILE *fp = fopen(filename, "w");
         if (!fp) {
            fprintf(stderr, "Unable to open image file: %s\n", filename);
            return 0;
         }
      
      // Print PPM image file 
      // First ppm pixel is top-left, so write in opposite scan-line order
         fprintf(fp, "P3\n");
         fprintf(fp, "%d %d\n", width, height);
         fprintf(fp, "255\n");
         for (int j = height-1; j >= 0 ; j--) {
            for (int i = 0; i < width; i++) {
               const R2Pixel& p = (*this)[i][j];
               int r = (int) (255 * p.Red());
               int g = (int) (255 * p.Green());
               int b = (int) (255 * p.Blue());
               fprintf(fp, "%-3d %-3d %-3d  ", r, g, b);
               if (((i+1) % 4) == 0) fprintf(fp, "\n");
            }
            if ((width % 4) != 0) fprintf(fp, "\n");
         }
         fprintf(fp, "\n");
      
      // Close file
         fclose(fp);
      }
      else {
      // Open file
         FILE *fp = fopen(filename, "wb");
         if (!fp) {
            fprintf(stderr, "Unable to open image file: %s\n", filename);
            return 0;
         }
      
      // Print PPM image file 
      // First ppm pixel is top-left, so write in opposite scan-line order
         fprintf(fp, "P6\n");
         fprintf(fp, "%d %d\n", width, height);
         fprintf(fp, "255\n");
         for (int j = height-1; j >= 0 ; j--) {
            for (int i = 0; i < width; i++) {
               const R2Pixel& p = (*this)[i][j];
               int r = (int) (255 * p.Red());
               int g = (int) (255 * p.Green());
               int b = (int) (255 * p.Blue());
               fprintf(fp, "%c%c%c", r, g, b);
            }
         }
      
      // Close file
         fclose(fp);
      }
   
   // Return success
      return 1;  
   }



////////////////////////////////////////////////////////////////////////
// JPEG I/O
////////////////////////////////////////////////////////////////////////


extern "C" { 
#   define XMD_H // Otherwise, a conflict with INT32
#   undef FAR // Otherwise, a conflict with windows.h
#   include "jpeg/jpeglib.h"
};



   int R2Image::
   ReadJPEG(const char *filename)
   {
   // Open file
      FILE *fp = fopen(filename, "rb");
      if (!fp) {
         fprintf(stderr, "Unable to open image file: %s\n", filename);
         return 0;
      }
   
   // Initialize decompression info
      struct jpeg_decompress_struct cinfo;
      struct jpeg_error_mgr jerr;
      cinfo.err = jpeg_std_error(&jerr);
      jpeg_create_decompress(&cinfo);
      jpeg_stdio_src(&cinfo, fp);
      jpeg_read_header(&cinfo, TRUE);
      jpeg_start_decompress(&cinfo);
   
   // Remember image attributes
      width = cinfo.output_width;
      height = cinfo.output_height;
      npixels = width * height;
      int ncomponents = cinfo.output_components;
   
   // Allocate pixels for image
      pixels = new R2Pixel [ npixels ];
      if (!pixels) {
         fprintf(stderr, "Unable to allocate memory for BMP file");
         fclose(fp);
         return 0;
      }
   
   // Allocate unsigned char buffer for reading image
      int rowsize = ncomponents * width;
      if ((rowsize % 4) != 0) rowsize = (rowsize / 4 + 1) * 4;
      int nbytes = rowsize * height;
      unsigned char *buffer = new unsigned char [nbytes];
      if (!buffer) {
         fprintf(stderr, "Unable to allocate temporary memory for JPEG file");
         fclose(fp);
         return 0;
      }
   
   // Read scan lines 
   // First jpeg pixel is top-left, so read pixels in opposite scan-line order
      while (cinfo.output_scanline < cinfo.output_height) {
         int scanline = cinfo.output_height - cinfo.output_scanline - 1;
         unsigned char *row_pointer = &buffer[scanline * rowsize];
         jpeg_read_scanlines(&cinfo, &row_pointer, 1);
      }
   
   // Free everything
      jpeg_finish_decompress(&cinfo);
      jpeg_destroy_decompress(&cinfo);
   
   // Close file
      fclose(fp);
   
   // Assign pixels
      for (int j = 0; j < height; j++) {
         unsigned char *p = &buffer[j * rowsize];
         for (int i = 0; i < width; i++) {
            double r, g, b, a;
            if (ncomponents == 1) {
               r = g = b = (double) *(p++) / 255;
               a = 1;
            }
            else if (ncomponents == 1) {
               r = g = b = (double) *(p++) / 255;
               a = 1;
               p++;
            }
            else if (ncomponents == 3) {
               r = (double) *(p++) / 255;
               g = (double) *(p++) / 255;
               b = (double) *(p++) / 255;
               a = 1;
            }
            else if (ncomponents == 4) {
               r = (double) *(p++) / 255;
               g = (double) *(p++) / 255;
               b = (double) *(p++) / 255;
               a = (double) *(p++) / 255;
            }
            else {
               fprintf(stderr, "Unrecognized number of components in jpeg image: %d\n", ncomponents);
               return 0;
            }
            R2Pixel pixel(r, g, b, a);
            SetPixel(i, j, pixel);
         }
      }
   
   // Free unsigned char buffer for reading pixels
      delete [] buffer;
   
   // Return success
      return 1;
   }


	

   int R2Image::
   WriteJPEG(const char *filename) const
   {
   // Open file
      FILE *fp = fopen(filename, "wb");
      if (!fp) {
         fprintf(stderr, "Unable to open image file: %s\n", filename);
         return 0;
      }
   
   // Initialize compression info
      struct jpeg_compress_struct cinfo;
      struct jpeg_error_mgr jerr;
      cinfo.err = jpeg_std_error(&jerr);
      jpeg_create_compress(&cinfo);
      jpeg_stdio_dest(&cinfo, fp);
      cinfo.image_width = width; 	/* image width and height, in pixels */
      cinfo.image_height = height;
      cinfo.input_components = 3;		/* # of color components per pixel */
      cinfo.in_color_space = JCS_RGB; 	/* colorspace of input image */
      cinfo.dct_method = JDCT_ISLOW;
      jpeg_set_defaults(&cinfo);
      cinfo.optimize_coding = TRUE;
      jpeg_set_quality(&cinfo, 75, TRUE);
      jpeg_start_compress(&cinfo, TRUE);
   
   // Allocate unsigned char buffer for reading image
      int rowsize = 3 * width;
      if ((rowsize % 4) != 0) rowsize = (rowsize / 4 + 1) * 4;
      int nbytes = rowsize * height;
      unsigned char *buffer = new unsigned char [nbytes];
      if (!buffer) {
         fprintf(stderr, "Unable to allocate temporary memory for JPEG file");
         fclose(fp);
         return 0;
      }
   
   // Fill buffer with pixels
      for (int j = 0; j < height; j++) {
         unsigned char *p = &buffer[j * rowsize];
         for (int i = 0; i < width; i++) {
            const R2Pixel& pixel = (*this)[i][j];
            int r = (int) (255 * pixel.Red());
            int g = (int) (255 * pixel.Green());
            int b = (int) (255 * pixel.Blue());
            if (r > 255) r = 255;
            if (g > 255) g = 255;
            if (b > 255) b = 255;
            *(p++) = r;
            *(p++) = g;
            *(p++) = b;
         }
      }
   
   
   
   // Output scan lines
   // First jpeg pixel is top-left, so write in opposite scan-line order
      while (cinfo.next_scanline < cinfo.image_height) {
         int scanline = cinfo.image_height - cinfo.next_scanline - 1;
         unsigned char *row_pointer = &buffer[scanline * rowsize];
         jpeg_write_scanlines(&cinfo, &row_pointer, 1);
      }
   
   // Free everything
      jpeg_finish_compress(&cinfo);
      jpeg_destroy_compress(&cinfo);
   
   // Close file
      fclose(fp);
   
   // Free unsigned char buffer for reading pixels
      delete [] buffer;
   
   // Return number of bytes written
      return 1;
   }






