//
//  main.m
//  MomentTensorModelTest
//
//  Created by Jeffrey J. Early on 5/22/14.
//  Copyright (c) 2014 Jeffrey J. Early. All rights reserved.
//

#import <Foundation/Foundation.h>
#import <GLNumericalModelingKit/GLNumericalModelingKit.h>
#import <GLNumericalModelingKit/GLOperationOptimizer.h>
#import "MomentTensorModels.h"
#import "MomentEllipseFileReader.h"
#import "EllipseErrorOperation.h"

int main(int argc, const char * argv[])
{

	@autoreleasepool {
	    
        GLEquation *equation = [[GLEquation alloc] init];
        
		GLDimension *tDim = [[GLDimension alloc] initDimensionWithGrid: kGLEndpointGrid nPoints: 10 domainMin: 0 length: 9*1800];
		GLFunction *t = [GLFunction functionOfRealTypeFromDimension: tDim withDimensions: @[tDim] forEquation: equation];
		GLScalar *kappa = [GLScalar scalarWithValue: 0.2014 forEquation: equation];
		GLScalar *sigma = [GLScalar scalarWithValue: 3.4342e-06 forEquation: equation];
		GLScalar *theta = [GLScalar scalarWithValue: -32.2493*M_PI/180 forEquation: equation];
		GLScalar *zeta = [GLScalar scalarWithValue: -3.0000e-06 forEquation: equation];
		//GLScalar *zeta = [GLScalar scalarWithValue: 1e-9 forEquation: equation];
		GLFloat Mxx0 = 1.8499e+05;
		GLFloat Myy0 = 2.7761e+05;
		GLFloat Mxy0 = -5.0412e+04;
		
		NSArray *tensorComps = strainVorticityDiffusivityModel( Mxx0, Myy0, Mxy0, t, kappa, sigma, theta, zeta);
		//NSArray *tensorComps = strainDiffusivityModel( Mxx0, Myy0, Mxy0, t, kappa, sigma, theta);
		GLFunction *Mxx=tensorComps[0]; GLFunction *Myy=tensorComps[1]; GLFunction *Mxy=tensorComps[2];
		for (int i=0; i<tDim.nPoints; i++) {
			printf("(Mxx,Myy,Mxy)=(%f, %f, %f)\n", Mxx.pointerValue[i], Myy.pointerValue[i], Mxy.pointerValue[i]);
		}
		printf("\n\n");
		NSArray *ellipseComps = tensorCompsToEllipseComps( tensorComps );
		GLFunction *semiMajor=ellipseComps[0]; GLFunction *semiMinor=ellipseComps[1]; GLFunction *angle=ellipseComps[2];
		for (int i=0; i<tDim.nPoints; i++) {
			printf("(a,b,theta)=(%f, %f, %f)\n", semiMajor.pointerValue[i], semiMinor.pointerValue[i], angle.pointerValue[i]*180/M_PI);
		}
	    
	}
    return 0;
}

