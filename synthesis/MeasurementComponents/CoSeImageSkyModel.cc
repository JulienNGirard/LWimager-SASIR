//# CoSeImageSkyModel.cc: Implementation of CoSeImageSkyModel class
//# Copyright (C) 1996,1997,1998,1999,2000,2001,2002,2003
//# Associated Universities, Inc. Washington DC, USA.
//#
//# This library is free software; you can redistribute it and/or modify it
//# under the terms of the GNU Library General Public License as published by
//# the Free Software Foundation; either version 2 of the License, or (at your
//# option) any later version.
//#
//# This library is distributed in the hope that it will be useful, but WITHOUT
//# ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
//# FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Library General Public
//# License for more details.
//#
//# You should have received a copy of the GNU Library General Public License
//# along with this library; if not, write to the Free Software Foundation,
//# Inc., 675 Massachusetts Ave, Cambridge, MA 02139, USA.
//#
//# Correspondence concerning AIPS++ should be addressed as follows:
//#        Internet email: aips2-request@nrao.edu.
//#        Postal address: AIPS++ Project Office
//#                        National Radio Astronomy Observatory
//#                        520 Edgemont Road
//#                        Charlottesville, VA 22903-2475 USA
//#
//# $Id$

#include <casa/Arrays/ArrayMath.h>
#include <casa/Arrays/Matrix.h>
#include <synthesis/MeasurementComponents/CoSeImageSkyModel.h>
#include <coordinates/Coordinates/DirectionCoordinate.h>
#include <images/Images/PagedImage.h>
#include <casa/OS/File.h>
#include <lattices/Lattices/LatticeExpr.h>
#include <lattices/Lattices/LatticeExprNode.h>
#include <lattices/Lattices/LatticeStepper.h>
#include <lattices/Lattices/LatticeIterator.h>
#include <synthesis/MeasurementEquations/SkyEquation.h>
#include <casa/Exceptions/Error.h>
#include <casa/BasicSL/String.h>
#include <casa/Utilities/Assert.h>

#include <casa/sstream.h>

#include <casa/Logging/LogMessage.h>
#include <casa/Logging/LogSink.h>
#include <casa/Logging/LogIO.h>

#include <msvis/MSVis/StokesVector.h>
#include <synthesis/MeasurementEquations/LatConvEquation.h>
#include <synthesis/MeasurementEquations/ClarkCleanLatModel.h>
#include <lattices/Lattices/SubLattice.h>
#include <lattices/Lattices/LCBox.h>
#include <CEA_comp_sens.h>

namespace casa {

Int CoSeImageSkyModel::add(ImageInterface<Float>& image,
			      const Int maxNumXfr)
{
  return CleanImageSkyModel::add(image, maxNumXfr);
};

Bool CoSeImageSkyModel::addMask(Int image,
				   ImageInterface<Float>& mask)
{
  return CleanImageSkyModel::addMask(image, mask);
};

Bool CoSeImageSkyModel::addResidual(Int image,
				       ImageInterface<Float>& residual) 
{
  return ImageSkyModel::addResidual(image, residual);
};

#define STORE_FLOAT
#ifdef STORE_FLOAT
static void store (const Matrix<Complex> &data, const string &name)
  {
    Matrix<Double> xform(2, 2);
    char buf[1024];
    Matrix<Float> dataf(data.shape()[0],data.shape()[1]);

    
    for (int i=0; i<data.shape()[0]; ++i)
      for (int j=0; j<data.shape()[1]; ++j) {
	dataf(i,j) = data(i,j).real();
      }
    
    xform = 0.0;
    xform.diagonal() = 1.0;
    Quantum<Double> incLon((8.0 / data.shape()(0)) * C::pi / 180.0, "rad");
    Quantum<Double> incLat((8.0 / data.shape()(1)) * C::pi / 180.0, "rad");
    Quantum<Double> refLatLon(45.0 * C::pi / 180.0, "rad");
    DirectionCoordinate dir(MDirection::J2000, Projection(Projection::SIN),
                            refLatLon, refLatLon, incLon, incLat,
                            xform, data.shape()(0) / 2, data.shape()(1) / 2);
    
    //cout<<"Saving... "<<name<<endl;
    Vector<Int> stokes(1);
    stokes(0) = Stokes::I;
    CoordinateSystem csys;
    csys.addCoordinate(dir);
    csys.addCoordinate(StokesCoordinate(stokes));
    csys.addCoordinate(SpectralCoordinate(casa::MFrequency::TOPO, 60e6, 0.0, 0.0, 60e6));
    PagedImage<Float> im(TiledShape(IPosition(4, data.shape()(0), data.shape()(1), 1, 1)), csys, name);
    im.putSlice(dataf, IPosition(4, 0, 0, 0, 0));
    
    // save map
    sprintf(buf,"%s.map",name.c_str());
    FILE *f = fopen(buf,"w");
    fprintf(f,"%d\n",(int)dataf.shape()[0]);
    for (int i=0; i<dataf.shape()[0]; ++i)
      for (int j=0; j<dataf.shape()[1]; ++j) {
	fprintf(f,"%f\n",dataf(i,j));
      }
     fclose(f);
}
#else
static void store (const Matrix<Complex> &data, const string &name)
  {
    Matrix<Double> xform(2, 2);
    char buf[1024];
 
    xform = 0.0;
    xform.diagonal() = 1.0;
    Quantum<Double> incLon((8.0 / data.shape()(0)) * C::pi / 180.0, "rad");
    Quantum<Double> incLat((8.0 / data.shape()(1)) * C::pi / 180.0, "rad");
    Quantum<Double> refLatLon(45.0 * C::pi / 180.0, "rad");
    DirectionCoordinate dir(MDirection::J2000, Projection(Projection::SIN),
                            refLatLon, refLatLon, incLon, incLat,
                            xform, data.shape()(0) / 2, data.shape()(1) / 2);
    
    //cout<<"Saving... "<<name<<endl;
    Vector<Int> stokes(1);
    stokes(0) = Stokes::I;
    CoordinateSystem csys;
    csys.addCoordinate(dir);
    csys.addCoordinate(StokesCoordinate(stokes));
    csys.addCoordinate(SpectralCoordinate(casa::MFrequency::TOPO, 60e6, 0.0, 0.0, 60e6));
    PagedImage<Complex> im(TiledShape(IPosition(4, data.shape()(0), data.shape()(1), 1, 1)), csys, name);
    im.putSlice(data, IPosition(4, 0, 0, 0, 0));
    
    // save map
    sprintf(buf,"%s.map",name.c_str());
    FILE *f = fopen(buf,"w");
    fprintf(f,"%d\n",(int)data.shape()[0]);
    for (int i=0; i<data.shape()[0]; ++i)
      for (int j=0; j<data.shape()[1]; ++j) {
	fprintf(f,"%f\n",data(i,j).real());
      }
     fclose(f);
}
#endif

static void store(SubLattice<Float> sl, const string& filename) {
  ofstream output;
  IPosition p(4, 0, 0, 0, 0);
  Matrix<Complex> M(sl.shape()[0], sl.shape()[1]);
  
  for (int xxx=0; xxx<sl.shape()[0]; ++xxx)
    for (int yyy=0; yyy<sl.shape()[1]; ++yyy) {
      p[0] = xxx; p[1] = yyy; 
      M(xxx,yyy) = sl.getAt(p); 
    }
  store(M, filename);
}

static void store(LatticeExpr<Float> sl, const string& filename) {
  ofstream output;
  IPosition p(4, 0, 0, 0, 0);
  Matrix<Complex> M(sl.shape()[0], sl.shape()[1]);
  
  for (int xxx=0; xxx<sl.shape()[0]; ++xxx)
    for (int yyy=0; yyy<sl.shape()[1]; ++yyy) {
      p[0] = xxx; p[1] = yyy; 
      M(xxx,yyy) = sl.getAt(p); 
    }
  store(M, filename);
}


static const char *name(const char *s, int i) {
  static char buf[1024];
  sprintf(buf,"%s_%d.img",s,i);
  return buf;
}


// Find the RMS of the image pixels and chop everything below that to
// hopefully leave more of the source
static void  chop_background(SubLattice<Float>& res) {
  float rms;
  IPosition pos(4,0,0,0,0);

  rms = 0;
  for (int i=0; i<res.shape()[0]; ++i)
    for (int j=0; j<res.shape()[1]; ++j) {
      pos[0] = i; pos[1] = j; 
      if ( res.getAt(pos) > 0 ) rms += res.getAt(pos)*res.getAt(pos);
    } 
  rms = sqrt(rms/res.shape()[0]/res.shape()[1]);
  for (int i=0; i<res.shape()[0]; ++i)
    for (int j=0; j<res.shape()[1]; ++j) {
      pos[0] = i; pos[1] = j;  
      if ( res.getAt(pos) < rms ) res.putAt(0,pos);
    } 

}     

 
Bool CoSeImageSkyModel::solve(SkyEquation& se) {
  cout << "Compressed Sensing\n";
	
  LogIO os(LogOrigin("CoSeImageSkyModel","compressed sensing"));

  
  if( numberIterations() < 1)
    return True;
  //Make the PSFs, one per field

  //os << LogIO::NORMAL    // Loglevel PROGRESS
  //   << "Making approximate Point Spread Functions" << LogIO::POST;
  if(!donePSF_p)
    makeApproxPSFs(se);
 
  // Loop over major cycles
  Int cycle=0;
  Bool stop_cs=False;

  if (displayProgress_p) {
    progress_p = new ClarkCleanProgress( pgplotter_p );
  } else {
    progress_p = 0;
  }
  
  float_array dirty_image, next_image, difference;
  double scale_factor=1; 
  while( !stop_cs ) {				
   
    ++cycle;
    
    makeNewtonRaphsonStep(se);

   
    if ( numberOfModels() > 1 ) {
    	cout << "Only one model allowed for compressed sensing\n"; exit(1);
    }
    
    Int model = 0;
	
    Int nx=image(model).shape()(0);
    Int ny=image(model).shape()(1);
    Int npol=image(model).shape()(2);
    Int nchan=image(model).shape()(3);
    AlwaysAssert((npol==1)||(npol==2)||(npol==4), AipsError);
	
    if ( nchan > 1 ) {
    	cerr << "Only 1 channel allowed for compressed sensing\n"; exit(1);
    }   
    
    Int chan = 0;
	    
    IPosition onePlane(4, nx, ny, 1, 1);
	    
    IPosition oneCube(4, nx, ny, npol, 1);
	    
    LCBox psfbox(IPosition(4, 0, 0, 0, chan), 
     IPosition(4, nx-1, ny-1, 0, chan),
     PSF(model).shape());
    SubLattice<Float>  psf_sl (PSF(model), psfbox, True);
		
    LCBox imagebox(IPosition(4, 0, 0, 0, chan), 
       IPosition(4, nx-1, ny-1, npol-1, chan), 
       residual(model).shape());
		
	
    SubLattice<Float> residual_sl (residual(model), imagebox, True);
    SubLattice<Float> image_sl (image(model), imagebox, True);
    SubLattice<Float> deltaimage_sl (deltaImage(model), imagebox, True);
    	  
    LatConvEquation eqn(psf_sl, residual_sl);
    ClarkCleanLatModel cleaner(deltaimage_sl);
 	 		 
    if ( cycle == 1 ) {
        store(residual_sl,"compressed_sensing_dirty.img");		
        store(psf_sl,"compressed_sensing_psf.img"); 
    }		
    // For major cycle compressed sensing
    IPosition pos(4, 0, 0, 0, 0);
    int algorithm; int transform; 
   
    AlwaysAssert(getMinimization() == "fista" || getMinimization() == "mca",  AipsError);
    if ( getMinimization() == "fista" ) algorithm = 2;    // FISTA 
    else if ( getMinimization() == "mca" ) algorithm = 3;    // MCA
    
    AlwaysAssert(getTransform() == "wavelets" || getTransform() == "curvelets", AipsError);
    if ( getTransform() == "wavelets" ) transform = 3;   // wavelets
    else if ( getTransform() == "curvelets" ) transform = 4;   // curvelets 
    
    if ( cycle == 1 ) {
        Matrix<Float> psf_patch;
        cleaner.getPsfPatch(psf_patch, eqn);   // Need this to scale down model
        scale_factor = 0;
        for (int i=0; i<psf_patch.shape()[0]; ++i)
            for (int j=0; j<psf_patch.shape()[1]; ++j) if ( psf_patch(i,j) > 0 ) scale_factor += psf_patch(i,j);
        cout << "[CS] PSF Scale Factor " << scale_factor << endl; fflush(stdout); 

	dirty_image.alloc(residual_sl.shape()[0]);  
	next_image.alloc(residual_sl.shape()[0]); difference.alloc(residual_sl.shape()[0]);
	chop_background(residual_sl);
	
	for (int i=0; i<residual_sl.shape()[0]; ++i)
	   for (int j=0; j<residual_sl.shape()[1]; ++j) {
	       pos[0] = i; pos[1] = j;
	       if (algorithm == 2 ) dirty_image[i][j] = 0;
	       else dirty_image[i][j] = residual_sl.getAt(pos); 
           }
	cout << "*****************" << endl;
 	cout << " CS parameters: " << endl;
	cout << "*****************" << endl;
	cout << " - Verbose=" << getVerbose() << endl;
	cout << " - AllOutput=" << getAllOutput() << endl;
	cout << " - doConvolve=" << getdoConvolve() << endl;
	cout << " - Number of iterations=" << numberIterations() << endl;
        cout << " - Positivity constraint=" << getPositivity() << endl;
        cout << " - Mu=" << gain() << " (Mu=0 --> estimation of Mu)" <<endl;
        cout << " - Threshold=" << threshold() << endl;
        cout << " - TolVar=" << getTolVar() << endl;
        cout << " - Num Scales=" << getNbScales() << endl;
        cout << " - Algorithm=" << algorithm << endl;
        cout << " - Transform=" << transform << endl;
        cout << " - No Coarse=" << getNoCoarse() << endl;
        cout << " - kSigma=" << getkSigma() << " (ksigma >0 --> MAD)" << endl;
        cout << "**************" << endl;

	compressed_sensing_begin(dirty_image, next_image, getVerbose(), numberIterations(), getPositivity(), /*mu*/gain(), threshold(), getTolVar(), getNbScales(), algorithm, transform, getNoCoarse(), 
getkSigma(),getAllOutput());
    } else {
	// Check the residual for stopping based on RMS
        // Also copy it
	double rms = 0;
       for (int i=0; i<residual_sl.shape()[0]; ++i)
	    for (int j=0; j<residual_sl.shape()[1]; ++j) {
		 pos[0] = i; pos[1] = j;
		 difference[i][j] = residual_sl.getAt(pos);
		 rms += (difference[i][j])*(difference[i][j]);
	    }
	rms = sqrt(rms/residual_sl.shape()[0]/residual_sl.shape()[1]); cout << " Residual RMS " << rms << endl; 
//	if ( false && rms <= sqrt(2)*threshold() )  { stop_cs = True; os << LogIO::NORMAL << "Compressed sensing reached threshold" << LogIO::POST; }
//	else stop_cs = !compressed_sensing_next(difference, next_image, algorithm);
stop_cs = !compressed_sensing_next(difference, next_image, algorithm);
  }
if (getAllOutput()) {
    store(image_sl,name("compressed_sensing_output",cycle)); store(residual_sl,name("compressed_sensing_residual",cycle));
}
    if ( !stop_cs ) { 		  
      for (int i=0; i<image_sl.shape()[0]; ++i)
      	for (int j=0; j<image_sl.shape()[1]; ++j) {
          pos[0] = i; pos[1] = j; 
          image_sl.putAt(next_image[i][j]/scale_factor,pos); 
       }
		
       store(image_sl,"compressed_sensing_output_temp.img"); store(residual_sl,"compressed_sensing_residual_temp.img");
    }		    
    modified_p = True; 		 

  }  //   while( ! stop_cs );
  if (progress_p) delete progress_p;
  
  os << LogIO::NORMAL    // Loglevel INFO
       << LatticeExprNode(sum(image(0))).getFloat() 
       << " Jy is the sum of the model " << LogIO::POST;

  os << LogIO::POST;

  return(True);
};


// Find maximum residual
Float CoSeImageSkyModel::maxField(Block<Vector<Float> >& imagemax,
				     Block<Vector<Float> >& imagemin) {

  LogIO os(LogOrigin("ImageSkyModel","maxField"));
  
  Float absmax=0.0;

  // Loop over all models
  for (Int model=0;model<numberOfModels();model++) {

    imagemax[model].resize(image(model).shape()(3));
    imagemin[model].resize(image(model).shape()(3));
    // Remember that the residual image can be either as specified
    // or created specially.
    ImageInterface<Float>* imagePtr=0;
    if(residual_p[model]) {
      imagePtr=residual_p[model];
    }
    else {
      imagePtr=(ImageInterface<Float> *)residualImage_p[model];
    }
    AlwaysAssert(imagePtr, AipsError);
    AlwaysAssert(imagePtr->shape().nelements()==4, AipsError);
    Int nx=imagePtr->shape()(0);
    Int ny=imagePtr->shape()(1);
    Int npol=imagePtr->shape()(2);
    
    AlwaysAssert((npol==1)||(npol==2)||(npol==4), AipsError);
    
    // Loop over all channels
    IPosition onePlane(4, nx, ny, 1, 1);
    LatticeStepper ls(imagePtr->shape(), onePlane, IPosition(4, 0, 1, 2, 3));
    RO_LatticeIterator<Float> imageli(*imagePtr, ls);
    
    // If we are using a mask then reset the region to be
    // cleaned
    Array<Float> maskArray;
    RO_LatticeIterator<Float> maskli;
    if(hasMask(model)) {
      Int mx=mask(model).shape()(0);
      Int my=mask(model).shape()(1);
      Int mpol=mask(model).shape()(2);
      //AlwaysAssert(mx==nx, AipsError);
      //AlwaysAssert(my==ny, AipsError);
      //AlwaysAssert(mpol==npol, AipsError);
      if((mx != nx) || (my != ny) || (mpol != npol)){
	throw(AipsError("Mask image shape is not the same as dirty image"));
      }
      LatticeStepper mls(mask(model).shape(), onePlane,
			 IPosition(4, 0, 1, 2, 3));
      
      maskli=RO_LatticeIterator<Float> (mask(model), mls);
      maskli.reset();
      if (maskli.cursor().shape().nelements() > 1) maskArray=maskli.cursor();
    }
    
    Int chan=0;
    Int polpl=0;
    Float imax, imin;
    imax=-1E20; imagemax[model]=imax;
    imin=+1E20; imagemin[model]=imin;
    imageli.reset();

    for (imageli.reset();!imageli.atEnd();imageli++) {
      imax=-1E20;
      imin=+1E20;
      IPosition imageposmax(imageli.cursor().ndim());
      IPosition imageposmin(imageli.cursor().ndim());
      
      // If we are using a mask then multiply by it
      if (hasMask(model)) {
	Array<Float> limage=imageli.cursor();
	//limage*=maskArray;
	minMaxMasked(imin, imax, imageposmin, imageposmax, limage, maskArray);
	maskli++;
	if (maskli.cursor().shape().nelements() > 1) maskArray=maskli.cursor();
      
      }
      
      else {
	minMax(imin, imax, imageposmin, imageposmax, imageli.cursor());
      }
      if(abs(imax)>absmax) absmax=abs(imax);
      if(abs(imin)>absmax) absmax=abs(imin);
      if(imin<imagemin[model][chan]) imagemin[model][chan]=imin;
      if(imax>imagemax[model][chan]) imagemax[model][chan]=imax;
      ++polpl;
      if(polpl==npol){
	++chan;
	polpl=0;	  
      }
    }
  }
  return absmax;
};
    

Vector<Float> CoSeImageSkyModel::outerMinMax(Lattice<Float> & lat, const uInt nCenter ) 
{
  Array<Float> arr = lat.get();
  IPosition pos( arr.shape() );
  uInt nx = lat.shape()(0);
  uInt ny = lat.shape()(1);
  uInt innerx = lat.shape()(0)/4;
  uInt innery = lat.shape()(1)/4;
  uInt nxc = 0;
  uInt nyc = 0;
  Float amax = 0.0;
  Vector<Float> amax2,minMax(2);
  //
  // First locate the location of the peak
  //
  for (uInt ix = 0; ix < nx; ix++) 
    for (uInt iy = 0; iy < ny; iy++) 
      if (arr(IPosition(4, ix, iy, 0, 0)) > amax) 
	{
	  nxc = ix;
	  nyc = iy;
	  amax = arr(IPosition(4, ix, iy, 0, 0));
	}
  //
  // Now exclude the mainlobe center on the location of the peak to
  // get the max. outer sidelobe.
  //
  Float myMax=0.0;
  Float myMin=0.0;

  uInt nxL = nxc - nCenter;
  uInt nxH = nxc + nCenter;
  uInt nyL = nyc - nCenter;
  uInt nyH = nyc + nCenter;
  uInt nx0 = nxc - innerx/2, nx1 = nxc + innerx/2;
  uInt ny0 = nyc - innery/2, ny1 = nyc + innery/2;
  
  //
  // Search only in the square with innerx and innery pixels on each side.
  //
  for (uInt ix = nx0; ix < nx1; ix++) {
    for (uInt iy = ny0; iy < ny1; iy++) {
      if ( !(ix >= nxL && ix <= nxH &&  iy >= nyL && iy <= nyH) ) {
	if (arr(IPosition(4, ix, iy, 0, 0)) > myMax) 
	  myMax = arr(IPosition(4, ix, iy, 0, 0));
	if (arr(IPosition(4, ix, iy, 0, 0)) < myMin)
	  myMin = arr(IPosition(4, ix, iy, 0, 0));
      }
    }
  }

  // Float absMax = max( abs(myMin), myMax );
  // return absMax;
  minMax(0) = myMin;
  minMax(1) = max( abs(myMin), myMax );
  return minMax;
};

} //#End casa namespace
